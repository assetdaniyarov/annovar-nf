/*
================================================================================
    ANNOVAR_DB workflow (--mode check_db | update_db)

    check_db : audit every database of the protocols (installed? intact?
               newer version? newer ClinVar at NCBI?) - read-only.
    update_db: audit ─► download missing / damaged official databases
                     ─► build the newest ClinVar from NCBI (for 'clinvar_latest')
                     ─► audit again and report;
               the report's gene table (GenCC, gnomAD constraint) is rebuilt
               when missing or older than 30 days (--gene_info).
    The protocols in params.annovar_protocols are the single source of truth
    for what has to be installed.
================================================================================
*/

include { DB_AUDIT                } from '../modules/local/annovar_db'
include { DB_AUDIT as DB_REAUDIT  } from '../modules/local/annovar_db'
include { DB_DOWNLOAD             } from '../modules/local/annovar_db'
include { DB_CLINVAR              } from '../modules/local/annovar_db'
include { DB_GENES                } from '../modules/local/annovar_db'

include { flag; protocolEntries; dbAuditLog } from '../subworkflows/local/utils_annovarnf'

workflow ANNOVAR_DB {

    take:
    builds      // list: ['hg19', 'hg38']
    update      // boolean

    main:
    ch_builds = channel.fromList(builds.collect { b -> [b, protocolEntries(b)] })

    DB_AUDIT(ch_builds, file(params.annovar_db), 'before')

    // The full audit tables are printed by the completion summary (from the
    // published TSVs): blocks taller than the terminal, printed while Nextflow's
    // live progress is on screen, get cut or duplicated when it redraws.
    if (update) {
        ch_rows = DB_AUDIT.out.tsv
            .flatMap { build, tsv ->
                log.info dbAuditLog(build, tsv, 'plan')
                tsv.splitCsv(header: true, sep: '\t')
            }

        // Official databases: install missing ones, replace damaged ones. Tables of one
        // package (ALL/EUR/... .sites from a 1000 Genomes .zip) are downloaded once.
        ch_download = ch_rows
            .filter { row -> row.action in ['download', 'redownload'] }
            .map { row -> [row.build, row.database, row.action == 'redownload', row.package ?: row.database] }
            .unique { row -> "${row[0]}:${row[3]}" }
            .map { build, database, force, _package -> [build, database, force] }

        DB_DOWNLOAD(ch_download, file(params.annovar_db), file(params.annovar_src))

        // ClinVar from NCBI: 'clinvar_latest' follows --clinvar_release, a pinned
        // 'clinvar_YYYYMMDD' is (re)built for exactly that date.
        ch_clinvar = ch_rows
            .filter { row -> row.action == 'clinvar' }
            .map { row ->
                def release = row.entry.endsWith('_latest') ? params.clinvar_release : row.entry.replaceFirst('clinvar_', '')
                [row.build, release, row.status == 'damaged']
            }
            .unique { row -> "${row[0]}:${row[1]}" }

        DB_CLINVAR(ch_clinvar, file(params.annovar_db))

        // Re-audit once every download/build has finished.
        ch_after = DB_DOWNLOAD.out.done.map { build, _db -> build }
            .mix(DB_CLINVAR.out.done.map { build, _db -> build })
            .collect()
            .ifEmpty([])
            .flatMap { _done -> builds.collect { b -> [b, protocolEntries(b)] } }

        DB_REAUDIT(ch_after, file(params.annovar_db), 'after')

        if (flag('gene_info')) DB_GENES(file(params.annovar_db), 30)
    }
}
