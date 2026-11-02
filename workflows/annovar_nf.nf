/*
================================================================================
    ANNOVAR_NF workflow

    input file ──► PREPARE_VCF ──► NORMALIZE_VCF ──► ANNOVAR_TABLE ──► GATHER_CHUNKS
                   (build, bgzip,   (per chunk:       (per chunk)        (per input file)
                    sort, index,     split multi-          │          GATHER_VCF (--save_vcf)
                    chunk plan)      allelics, drop        ▼
                        │            '*' / <SYM>)   ANNOTATION_SUMMARY
                        ├──► SAMPLE_QC ───────────► (per input file: QC,
                        │    (sex, relatedness, PED)  inheritance, P/LP, VUS)
                        └──► DB_INSTALL: protocol databases of the build not installed yet
                             (--auto_install_db, once per build); its chunks wait for it
    DB_GENES: the report's gene table (GenCC, gnomAD constraint) when missing;
              ANNOTATION_SUMMARY waits for it
================================================================================
*/

include { PREPARE_VCF        } from '../modules/local/prepare_vcf'
include { NORMALIZE_VCF      } from '../modules/local/normalize_vcf'
include { ANNOVAR_TABLE      } from '../modules/local/annovar_table'
include { GATHER_CHUNKS      } from '../modules/local/gather_chunks'
include { GATHER_VCF         } from '../modules/local/gather_vcf'
include { ANNOTATION_SUMMARY } from '../modules/local/annotation_summary'
include { SAMPLE_QC          } from '../modules/local/sample_qc'
include { DB_INSTALL         } from '../modules/local/annovar_db'
include { DB_GENES           } from '../modules/local/annovar_db'

include { annovarProtocol; annovarVersion; autoInstall; flag; chunkCpus; chunkParts; codeChecksum; checkAnnovarResources; databaseVersions; databaseTable; referenceFasta; referenceFai; fileDetectedLog; fileDoneLog; chunkProgressLog; geneTable; genePanels; installLog; liveDisplay; missingDatabases; protocolEntries } from '../subworkflows/local/utils_annovarnf'

workflow ANNOVAR_NF {

    take:
    ch_input    // channel: [ meta, vcf, index_or_[] ]

    main:
    ch_versions = channel.empty()

    //
    // 1. Build detection, BGZF/sort/index, chunk plan
    //
    PREPARE_VCF(ch_input, codeChecksum(['detect_genome_build.py', 'plan_chunks.py']), chunkParts())
    ch_versions = ch_versions.mix(PREPARE_VCF.out.versions)

    //
    // 2. ANNOVAR databases of every detected build: the protocol databases that
    //    are not installed yet are installed first (--auto_install_db, e.g. on a
    //    new humandb), once per build; files of that build wait for it.
    //
    ch_db = PREPARE_VCF.out.vcf
        .map { _meta, _vcf, _index, build, _source, _plan -> build }
        .unique()
        .branch { build ->
            install: autoInstall() && missingDatabases(build)
            ready: true
        }
    ch_install = ch_db.install.map { build ->
        log.info installLog(build, missingDatabases(build))
        [build, protocolEntries(build)]
    }
    DB_INSTALL(ch_install, file(params.annovar_db), file(params.annovar_src))
    ch_ready = ch_db.ready.mix(DB_INSTALL.out.build)

    //
    // 3. One item per chunk. The detected build is stored in meta.genome and
    //    the ANNOVAR databases for it are verified before any annotation runs.
    //    The chunks of a file share its --threads CPUs (meta.cpus: one each
    //    when the file has a chunk per CPU) and run at the same time.
    //
    ch_chunks = PREPARE_VCF.out.vcf
        .map { meta, vcf, index, build, source, chunk_file -> [build, meta, vcf, index, source, chunk_file] }
        .combine(ch_ready, by: 0)
        .flatMap { build, meta, vcf, index, source, chunk_file ->
            def rows = chunk_file.readLines().findAll { line -> line.trim() }.collect { line -> line.split('\t') }
            if (!rows) {
                log.warn "[annovar-nf] ${meta.id}: no variant records, skipping"
                return []
            }
            checkAnnovarResources(build)
            def counts = rows.collect { r -> r[3].toLong() }
            log.info fileDetectedLog(meta, build, source, rows.size(), counts.sum() as long)
            def file_meta = meta + [genome: build, n_chunks: rows.size()]
            def cpus = chunkCpus(counts)
            rows.withIndex().collect { r, i ->
                def chunk_meta = file_meta + [chunk_idx: r[0].toInteger(), chunk: r[1], regions: r[2], cpus: cpus[i]]
                [chunk_meta, vcf, index, referenceFasta(build), referenceFai(build)]
            }
        }

    //
    // 4. Normalize each chunk
    //
    NORMALIZE_VCF(ch_chunks)
    ch_versions = ch_versions.mix(NORMALIZE_VCF.out.versions.first())

    //
    // 5. ANNOVAR per chunk
    //
    ch_annovar_in = NORMALIZE_VCF.out.vcf
        .map { meta, vcf ->
            def protocol = annovarProtocol(meta.genome)
            [meta, vcf, protocol[0], protocol[1]]
        }

    ANNOVAR_TABLE(
        ch_annovar_in,
        file(params.annovar_db, checkIfExists: true),
        file(params.annovar_src, checkIfExists: true),
        codeChecksum(['format_annovar_table.py'])
    )
    ch_versions = ch_versions.mix(ANNOVAR_TABLE.out.versions.first())

    // Progress of files split into several chunks (about every 10%), for plain
    // logs only: Nextflow's live display already counts the ANNOVAR_TABLE tasks
    // on one line, and messages there cannot be updated in place.
    // subscribe() handles items one at a time, so the counter needs no locking.
    def chunks_done = [:]
    def plain_log = !liveDisplay()
    ANNOVAR_TABLE.out.chunk
        .filter { meta, _has, _h, _b -> plain_log && meta.n_chunks > 1 }
        .subscribe { meta, _has, _h, _b ->
            def n = (chunks_done[meta.id] ?: 0) + 1
            chunks_done[meta.id] = n
            def step = [1, (meta.n_chunks / 10).intValue()].max()
            if (n % step == 0 || n == meta.n_chunks) log.info chunkProgressLog(meta, n)
        }


    //
    // 6. Gather chunks per input file, in genomic order. groupKey lets each
    //    file proceed as soon as its last chunk finishes.
    //
    ch_gather = ANNOVAR_TABLE.out.chunk
        .map { meta, _has_variants, header, body ->
            def file_meta = meta.subMap(['id', 'genome', 'input', 'n_chunks'])
            [groupKey(file_meta, meta.n_chunks), meta.chunk_idx, header, body]
        }
        .groupTuple()
        .map { key, idx, headers, bodies ->
            def order = (0..<idx.size()).toList().sort { i -> idx[i] }
            [key.getGroupTarget(), order.collect { i -> headers[i] }, order.collect { i -> bodies[i] }]
        }

    GATHER_CHUNKS(ch_gather)
    ch_versions = ch_versions.mix(GATHER_CHUNKS.out.versions.first())
    GATHER_CHUNKS.out.txt.subscribe { meta, txt -> log.info fileDoneLog(meta, txt) }

    // Annotated VCF (optional): ANNOVAR_TABLE emits per-chunk VCFs only with --save_vcf.
    ch_vcf_out = channel.empty()
    if (flag('save_vcf')) {
        ch_gather_vcf = ANNOVAR_TABLE.out.vcf
            .map { meta, has_variants, vcf ->
                def file_meta = meta.subMap(['id', 'genome', 'input', 'n_chunks'])
                [groupKey(file_meta, meta.n_chunks), meta.chunk_idx, has_variants == 'true', vcf]
            }
            .groupTuple()
            .map { key, idx, has_variants, vcfs ->
                def order = (0..<idx.size()).toList().sort { i -> idx[i] }
                // Empty chunks carry a header-only VCF; leave them out of the
                // concat unless every chunk is empty.
                def vcf_order = order.findAll { i -> has_variants[i] } ?: [order[0]]
                [key.getGroupTarget(), vcf_order.collect { i -> vcfs[i] }]
            }
        GATHER_VCF(ch_gather_vcf)
        ch_versions = ch_versions.mix(GATHER_VCF.out.versions.first())
        ch_vcf_out = GATHER_VCF.out.vcf
    }

    ch_run_builds = PREPARE_VCF.out.vcf
        .map { _meta, _vcf, _index, build, _source, _chunks -> build }
        .unique()

    //
    // 7. Per-sample QC of every input file: sex from chrX, relatedness of the
    //    samples of a family, checks against the PED file
    //
    ch_sample_qc = channel.empty()
    if (!flag('skip_sample_qc')) {
        SAMPLE_QC(
            PREPARE_VCF.out.vcf.map { meta, vcf, index, build, _source, _chunks -> [meta, vcf, index, build] },
            params.ped ? file(params.ped, checkIfExists: true) : [],
            codeChecksum(['sample_qc.py'])
        )
        ch_versions = ch_versions.mix(SAMPLE_QC.out.versions.first())
        ch_sample_qc = SAMPLE_QC.out.json
    }

    //
    // 8. Summary report per input file (+ combined one for several files),
    //    including how the genome build of each file was determined, the
    //    sample QC, the inheritance analysis and the ANNOVAR databases used
    //
    ch_db_table = ch_run_builds
        .map { build -> databaseTable(build) }
        .collectFile(name: 'databases.tsv', sort: true, newLine: true)

    // The gene table (diseases, modes of inheritance, LOEUF): built in humandb
    // while the files are annotated when it is missing (--auto_install_db).
    ch_genes = channel.value([])
    if (flag('gene_info')) {
        if (geneTable().exists()) {
            ch_genes = channel.value(geneTable())
        } else if (autoInstall()) {
            DB_GENES(file(params.annovar_db), '')
            ch_genes = DB_GENES.out.done.map { _done -> geneTable().exists() ? geneTable() : [] }
        } else {
            log.warn "[annovar-nf] no gene table in ${params.annovar_db}: the report shows no gene information " +
                "(--mode update_db builds it)"
        }
    }
    ANNOTATION_SUMMARY(ANNOVAR_TABLE.out.stats.collect(), PREPARE_VCF.out.build_report.collect(),
                       ch_sample_qc.collect().ifEmpty([]), ch_db_table, ch_genes, genePanels(), annovarVersion(),
                       codeChecksum(['summarize_annotations.py', 'inheritance.py', 'gene_info.py']))
    ch_versions = ch_versions.mix(ANNOTATION_SUMMARY.out.versions)

    //
    // Software versions (one entry per process) and the ANNOVAR databases
    // used for every genome build in the run
    //
    ch_db_versions = ch_run_builds.map { build -> databaseVersions(build) }

    ch_versions
        .map { f -> f.text.trim() }
        .unique()
        .mix(ch_db_versions)
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true, newLine: true)

    emit:
    txt     = GATHER_CHUNKS.out.txt
    vcf     = ch_vcf_out
    summary = ANNOTATION_SUMMARY.out.html
}
