/*
================================================================================
    annovar-nf helper functions
    Parameter validation, input discovery, sample naming and ANNOVAR protocol
    handling. Kept free of channel logic so they can be reused and tested.
================================================================================
*/

// File extensions recognised as variant files when scanning a directory.
// Index files (.tbi/.csi/.idx) are never matched.
def VARIANT_FILE_RE() { ~/(?i)\.(g\.)?(vcf|gvcf|bcf)(\.b?gz)?$/ }

// Pipeline-specific suffixes stripped from file names to derive sample ids,
// applied repeatedly: HG002.hard-filtered.FINAL.vcf.gz -> HG002.
def SAMPLE_SUFFIX_RE() { ~/(?i)\.(hard-filtered|hard_filtered|filtered|FINAL|PASS|snp_indel|norm|normalized)$/ }

// Gene-based ANNOVAR databases ('g' operation) used to infer --operation.
def GENE_DB_RE() { ~/(?i)^(refGene|knownGene|ensGene|ccdsGene|ncbiRefSeq\w*|wgEncodeGencode\w+|refGeneWithVer|ensGene\d*)$/ }

/*
 * Help ------------------------------------------------------------------------
 */
def helpMessage() {
    return """
    ${workflow.manifest.name} v${workflow.manifest.version}
    ${workflow.manifest.description}

    Usage:
      nextflow run main.nf --input <samplesheet.csv | dir | 'glob' | file.vcf.gz> \\
          --annovar_db /path/to/humandb --annovar_src /path/to/annovar [options]

    Input / output
      --input             Samplesheet (.csv/.tsv: sample,vcf[,genome]), a directory,
                          a quoted glob ('/data/*.vcf.gz') or a single VCF/BCF file
      --input_pattern     File pattern used when --input is a directory [${params.input_pattern}]
      --recursive         Also search sub-directories of an input directory [${params.recursive}]
      --outdir            Output directory [${params.outdir}]

    Genome build
      --genome            auto | hg19 | hg38 [${params.genome}]
                          auto = per-file detection from ##contig lengths, header and file name
      --skip_build_check  Downgrade build contradictions to warnings [${params.skip_build_check}]

    ANNOVAR
      --annovar_db        ANNOVAR humandb directory [required]
      --annovar_src       Directory containing table_annovar.pl [required]
      --protocol          Comma-separated protocol overriding the per-build defaults
      --operation         Matching operations (g,f,r,...); inferred from --protocol if omitted
      --onetranscript     Report one transcript per gene [${params.onetranscript}]
      --annovar_args      Extra arguments appended to table_annovar.pl ['${params.annovar_args}']
      --threads           CPUs for ANNOVAR per input file: one chunk and process each [${params.threads}]
      --skip_db_check     Do not verify that every database file exists [${params.skip_db_check}]
      --auto_install_db   Install missing protocol databases before annotating (ANNOVAR
                          packages, ClinVar from NCBI) [${params.auto_install_db}]

    Database maintenance
      --mode              annotate | check_db | update_db [${params.mode}]
                          check_db : read-only audit of every protocol database (installed,
                                     intact, newer version, newer ClinVar at NCBI)
                          update_db: install missing/damaged databases, build the newest
                                     ClinVar from NCBI for 'clinvar_latest', then re-audit
      --db_builds         Builds to check/update, comma-separated [all builds with a protocol]
      --clinvar_release   ClinVar release built for 'clinvar_latest': latest | YYYYMMDD [${params.clinvar_release}]
      --db_download_forks Parallel database downloads [${params.db_download_forks}]
      --db_offline        Audit without contacting ANNOVAR/NCBI [${params.db_offline}]

    Pre-processing (bcftools)
      --chunk_size        At most about this many records per chunk; 0 = never split [${params.chunk_size}]
      --pass_only         Keep FILTER=PASS records only [${params.pass_only}]
      --drop_genotypes    Annotate sites only (drop FORMAT/sample columns) [${params.drop_genotypes}]
      --min_ac            Drop alleles with AC below this after splitting (0 = keep all) [${params.min_ac}]
      --fasta_hg19        Reference FASTA for left-alignment of hg19 inputs (optional)
      --fasta_hg38        Reference FASTA for left-alignment of hg38 inputs (optional)
      --check_ref         bcftools norm --check-ref mode when a FASTA is given [${params.check_ref}]

    Samples, families and genotype quality
      --ped               PED file (family, sample, father, mother, sex 1/2, phenotype 2 = affected),
                          or its CSV / TSV spreadsheet export: inheritance analysis for the affected
                          samples, sex and relationship checks
      --skip_sample_qc    No per-sample QC (sex from chrX, relatedness) [${params.skip_sample_qc}]
      --family_max_samples  Files with at most this many samples get all genotypes listed and the
                          inheritance analysis [${params.family_max_samples}]
      --gene_info         Diseases and modes of inheritance (GenCC) and LOEUF (gnomAD) of every gene
                          in the report; the gene table is built in humandb [${params.gene_info}]
      --gene_panel        Gene panel(s), comma-separated: gene lists or PanelApp TSV exports; other
                          genes are shown on request only [${params.gene_panel ?: 'none'}]
      --qc_min_dp         Flag LOW_DP below this read depth [${params.qc_min_dp}]
      --qc_min_gq         Flag LOW_GQ below this genotype quality [${params.qc_min_gq}]
      --qc_min_alt_reads  Flag LOW_ALT_READS for carriers with fewer variant reads [${params.qc_min_alt_reads}]
      --qc_het_vaf        Flag VAF_SKEW for heterozygous calls outside MIN,MAX [${params.qc_het_vaf}]
      --qc_min_hom_vaf    Flag VAF_SKEW for homozygous / hemizygous calls below [${params.qc_min_hom_vaf}]

    Outputs
      --compress_output   bgzip the annotated tables [${params.compress_output}]
      --save_vcf          Also write the annotated VCF (+ .tbi) [${params.save_vcf}]
      --save_normalized   Publish the normalized per-chunk VCFs [${params.save_normalized}]
      --publish_dir_mode  Nextflow publish mode [${params.publish_dir_mode}]
      --monochrome_logs   Disable coloured terminal messages [${params.monochrome_logs}]

    Profiles (-profile)
      standard (local, default) | slurm | conda | mamba | docker | singularity | apptainer
      test (end-to-end with a mock ANNOVAR, no licence needed)
    """.stripIndent()
}

/*
 * Terminal output -----------------------------------------------------------------
 * Colours are disabled with --monochrome_logs (e.g. when logging to a file).
 */
def ansi(String code) {
    return flag('monochrome_logs') ? '' : "\u001B[${code}m"
}

def colors() {
    return [reset: ansi('0'), bold: ansi('1'), dim: ansi('2'), red: ansi('0;31'), green: ansi('0;32'),
            yellow: ansi('0;33'), blue: ansi('0;34'), purple: ansi('0;35'), cyan: ansi('0;36')]
}

// 1234567 -> 1,234,567 (locale independent); non-numbers are returned unchanged.
def groupDigits(value) {
    def s = value.toString()
    return s ==~ /\d{4,}/ ? s.reverse().replaceAll(/(\d{3})(?=\d)/, '$1,').reverse() : s
}

def rule(String colour) {
    return "${colour}${'─' * 72}${colors().reset}"
}

// True when Nextflow's live (ANSI) display is on. It redraws every message it shows in
// place; a block taller than the terminal cannot be redrawn and leaves partial copies
// behind, so long tables are shortened there (full versions are in the report files).
// [rows, columns] of the terminal (Nextflow's own jline); 24 x 80 when unknown.
def terminalSize() {
    try {
        def t = jline.TerminalFactory.get()
        return [Math.max(t.getHeight(), 10), Math.max(t.getWidth(), 40)]
    } catch (Exception _e) {
        return [24, 80]
    }
}

def liveDisplay() {
    def console = System.console()
    if (!nextflow.Global.session?.ansiLog || console == null) return false
    try {
        return console.isTerminal()     // Java 22+: a console object exists even when redirected
    } catch (Exception _e) {
        return true
    }
}

// Checksum of bin/ scripts, given to the processes that run them: Nextflow does not hash
// bin/ into the task, so -resume would otherwise reuse results made by older scripts.
def codeChecksum(List names) {
    return names.collect { n -> file("${projectDir}/bin/${n}").text }.join('\n').md5()
}

// ANNOVAR release = newest script date: a release does not update every script
// (2025-03-02 ships annotate_variation.pl 2025-03-02 with table_annovar.pl 2022-08-02).
def annovarVersion() {
    def dates = ['table_annovar.pl', 'annotate_variation.pl', 'coding_change.pl'].collect { name ->
        def script = file("${params.annovar_src}/${name}")
        script.exists() ? script.text.find(/Date: (\d{4}-\d{2}-\d{2})/) { _match, date -> date } : null
    }.findAll { d -> d }
    return dates ? dates.max() : 'unknown version'
}

def protocolBuilds() {
    return params.protocol ? ['hg19', 'hg38'] : ['hg19', 'hg38'].findAll { g -> params.annovar_protocols?.get(g) }
}

def databaseSummary() {
    return params.protocol
        ? "custom protocol (${params.protocol.toString().split(',').size()} databases)"
        : protocolBuilds().collect { g -> "${g}: ${params.annovar_protocols.get(g).size()}" }.join(' · ')
}

/*
 * Banner rows for what '<family>_latest' protocol entries resolve to, one per alias:
 * ['ClinVar', 'clinvar_latest → clinvar_20261025 (hg19, hg38)']; a note goes after the ClinVar row.
 */
def latestAliasRows(String clinvarNote) {
    def c = colors()
    def resolved = [:]
    protocolBuilds().each { g ->
        protocolPairs(g).findAll { p -> p[0].endsWith('_latest') }.each { p ->
            if (!resolved.containsKey(p[0])) resolved[p[0]] = [:]
            resolved[p[0]][g] = resolveDatabase(g, p[0]) ?: "${c.red}not installed${c.reset}"
        }
    }
    return resolved.collect { alias, dbs ->
        def text = dbs.values().toSet().size() == 1
            ? "${alias} → ${dbs.values().first()}  ${c.dim}(${dbs.keySet().join(', ')})${c.reset}"
            : "${alias} → " + dbs.collect { g, db -> "${g} ${db}" }.join(' · ')
        def clinvar = alias == 'clinvar_latest'
        [clinvar ? 'ClinVar' : alias.replace('_latest', ''), text + (clinvar && clinvarNote ? "  ${c.dim}· ${clinvarNote}${c.reset}" : '')]
    }
}

/*
 * Start-up banner with the run configuration and the input files found.
 */
def pipelineBanner(List inputs) {
    def c = colors()
    def ids = inputs.collect { e -> e[0].id }
    def shown = ids.take(4).join(', ') + (ids.size() > 4 ? ", … (+${ids.size() - 4})" : '')
    def genome = normalizeGenome(params.genome) == 'auto' ? 'auto — detected per file' : normalizeGenome(params.genome)
    def rows = [
        ['Input files',  "${c.bold}${ids.size()}${c.reset}  ${c.dim}(${shown})${c.reset}"],
        ['Genome build', genome],
        ['Output',       file(params.outdir).toString()],
        ['ANNOVAR',      "${params.annovar_src}  ${c.dim}(${annovarVersion()})${c.reset}"],
        ['Databases',    databaseSummary()],
    ] + latestAliasRows('') + [
        ['Processing',   "${annovarCpus()} CPU${annovarCpus() == 1 ? '' : 's'} per input file, one chunk each" +
                         (annovarCpus() < intParam(params.threads) ? " ${c.dim}(--threads ${params.threads} capped by --max_cpus)${c.reset}" : '') +
                         " · chunks of at most about ${groupDigits(params.chunk_size)} records"],
        ['Filters',      "PASS only: ${flag('pass_only') ? 'yes' : 'no'} · genotypes: ${flag('drop_genotypes') ? 'dropped' : 'kept'}" +
                         (intParam(params.min_ac) > 0 ? " · min AC: ${params.min_ac}" : '')],
        ['Genes',        geneTableSummary()],
    ] + (genePanels() ? [['Gene panel', genePanels().collect { p -> p.name }.join(', ')]] : []) + [
        ['Run',          "${workflow.runName} · profile ${workflow.profile} · Nextflow ${workflow.nextflow.version}"],
        ['Work dir',     workflow.workDir.toString()],
    ]
    def body = rows.collect { r -> "  ${c.dim}${r[0].padRight(13)}${c.reset}${r[1]}" }.join('\n')
    return """
${rule(c.blue)}
  ${c.bold}${c.cyan}annovar-nf${c.reset} ${c.bold}v${workflow.manifest.version}${c.reset}  ${c.dim}·  ANNOVAR annotation of VCF/BCF files (hg19 / hg38)${c.reset}
${rule(c.blue)}
${body}
${rule(c.blue)}"""
}

/*
 * Start-up banner for --mode check_db / update_db.
 */
def dbBanner() {
    def c = colors()
    def update = params.mode == 'update_db'
    def rows = [
        ['Mode',         update ? "${c.yellow}update_db${c.reset} — installs into humandb" : 'check_db — read-only'],
        ['Builds',       dbBuilds().join(', ')],
        ['humandb',      params.annovar_db.toString()],
        ['Databases',    databaseSummary()],
    ] + latestAliasRows("update_db builds NCBI release: ${params.clinvar_release}") + [
        ['Versions',     flag('db_offline') ? 'offline (integrity checks only)' : 'ANNOVAR download page + NCBI ClinVar'],
        ['Genes',        geneTableSummary() + (update && flag('gene_info') ? ' · rebuilt when older than 30 days' : '')],
    ]
    if (update) rows << ['ANNOVAR', "${params.annovar_src}  ${c.dim}(${annovarVersion()})${c.reset}"]
    def body = rows.collect { r -> "  ${c.dim}${r[0].padRight(13)}${c.reset}${r[1]}" }.join('\n')
    return """
${rule(c.blue)}
  ${c.bold}${c.cyan}annovar-nf${c.reset} ${c.bold}v${workflow.manifest.version}${c.reset}  ${c.dim}·  ANNOVAR database ${update ? 'update' : 'check'}${c.reset}
${rule(c.blue)}
${body}
${rule(c.blue)}"""
}

/*
 * One line per input file once its build is known and its chunks are planned.
 */
def fileDetectedLog(Map meta, String build, String source, int nChunks, long nRecords) {
    def c = colors()
    def sources = [contig_lengths: 'build from ##contig lengths', explicit: 'build set by user',
                   header_hints: 'build from header hints', filename: "${c.yellow}build from file name only — check it${c.reset}"]
    def chunks = nChunks == 1 ? '1 chunk' : "${nChunks} chunks"
    return "  ${c.blue}▸${c.reset} ${c.bold}${meta.id.padRight(24)}${c.reset} ${c.cyan}${build}${c.reset}" +
        "  ${groupDigits(nRecords).padLeft(11)} input records · ${chunks.padRight(9)}  ${c.dim}${sources.get(source) ?: source}${c.reset}"
}

/*
 * Databases installed before the files of a build are annotated (--auto_install_db).
 */
def installLog(String build, List missing) {
    def c = colors()
    return "  ${c.yellow}▾${c.reset} ${c.bold}${build}${c.reset}: installing ${missing.size()} missing database(s) into " +
        "${params.annovar_db} before annotation ${c.dim}(once; large ones take hours; log: " +
        "${params.outdir}/database_check/logs/${build}.db_install.log)${c.reset}\n" +
        "    ${c.dim}${missing.join(', ')}${c.reset}"
}

/*
 * One line per input file once its annotated table has been written.
 */
def fileDoneLog(Map meta, txt) {
    def c = colors()
    return "  ${c.green}✔${c.reset} ${c.bold}${meta.id.padRight(24)}${c.reset} ${c.dim}→${c.reset} ${params.outdir}/annotation/${meta.id}/${txt.name}"
}

/*
 * Genes with carried pathogenic / likely pathogenic variants (ClinVar, plus
 * InterVar-only calls), per input file.
 */
def plpGenesLog(List ids, int maxFiles) {
    def c = colors()
    def rows = summaryFiles(ids, '_plp.tsv').collectMany { f -> f.splitCsv(header: true, sep: '	') }
    if (!rows) return []
    def lines = ["  ${c.bold}Pathogenic / likely pathogenic, carried${c.reset}"]
    // One TSV row per variant and carrier sample: count each variant once per file.
    def bySample = rows.groupBy { r -> r.file }.collectEntries { file_id, rs -> [file_id, rs.unique { r -> r.vcf }] }
    bySample.take(maxFiles).each { sample, vs ->
        def genes = vs.countBy { r -> r.gene }.sort { a, b -> b.value <=> a.value ?: a.key <=> b.key }
        def shown = genes.take(8).collect { g, n -> n > 1 ? "${c.red}${g}${c.reset}(${n})" : "${c.red}${g}${c.reset}" }
        def more = genes.size() > 8 ? " ${c.dim}+${genes.size() - 8} genes${c.reset}" : ''
        lines << "  ${c.dim}${sample.padRight(24)}${c.reset} ${vs.size()} variants · ${shown.join(', ')}${more}"
    }
    if (bySample.size() > maxFiles)
        lines << "  ${c.dim}… ${bySample.size() - maxFiles} more input files: ${params.outdir}/summary/annotation_summary_plp.tsv${c.reset}"
    lines << ''
    return lines
}

/*
 * Progress line for files split into several chunks (large cohorts).
 */
def chunkProgressLog(Map meta, int done) {
    def c = colors()
    def pct = (100 * done / meta.n_chunks).intValue()
    def width = 24
    def filled = (width * done / meta.n_chunks).intValue()
    def bar = "${c.cyan}${'█' * filled}${c.reset}${c.dim}${'░' * (width - filled)}${c.reset}"
    return "  ${c.blue}⋯${c.reset} ${c.bold}${meta.id.padRight(24)}${c.reset} ${bar} ${done}/${meta.n_chunks} chunks annotated (${pct}%)"
}

// Published per-file summary files of this run (<id>.annotation_summary<suffix>);
// empty ones (-stub runs) are skipped.
def summaryFiles(List ids, String suffix) {
    return ids.collect { id -> file("${params.outdir}/summary/${id}.annotation_summary${suffix}") }
        .findAll { f -> f.exists() && f.size() > 0 }
}

/*
 * Final report: status, run statistics and a per-file results table.
 * ids: input-file ids of this run (empty in database modes).
 */
def completionSummary(List ids) {
    def c = colors()
    def stats = workflow.stats
    def lines = []
    if (workflow.success) {
        lines << rule(c.green)
        lines << "  ${c.green}${c.bold}✔ annovar-nf completed successfully${c.reset}"
        lines << rule(c.green)
    } else {
        lines << rule(c.red)
        lines << "  ${c.red}${c.bold}✘ annovar-nf failed${c.reset}"
        lines << rule(c.red)
    }
    lines << "  ${c.dim}Duration     ${c.reset}${workflow.duration}  ${c.dim}(CPU time ${stats.computeTimeFmt.replaceAll(/[()]/, '')})${c.reset}"
    lines << "  ${c.dim}Tasks        ${c.reset}${stats.succeedCount} succeeded · ${stats.cachedCount} cached · " +
        "${stats.failedCount > 0 ? c.red : ''}${stats.failedCount} failed${c.reset}"

    if (!workflow.success) {
        def message = (workflow.errorMessage ?: 'see the log for details').readLines().findAll { l -> l.trim() }.take(6)
        lines << "  ${c.dim}Error        ${c.reset}${c.red}${message.join("\n               ")}${c.reset}"
        lines << "  ${c.dim}Log          ${c.reset}${workflow.launchDir}/.nextflow.log"
        lines << "  ${c.dim}Next         ${c.reset}fix the problem and rerun the same command with ${c.bold}-resume${c.reset}"
        lines << rule(c.red)
        return '\n' + lines.join('\n')
    }

    def live = liveDisplay()
    if (params.mode != 'annotate') {
        // Audit tables, read back from the published reports.
        def update = params.mode == 'update_db'
        def report = update ? 'after' : 'before'
        def stage = update ? 'result' : 'check'
        def reports = dbBuilds().collect { b -> file("${params.outdir}/database_check/${b}.db_audit.${report}.tsv") }
            .findAll { f -> f.exists() }
        def render = { String layout, int width -> reports.collect { f -> dbAuditLog(f.name.tokenize('.')[0], f, stage, layout, width) } }
        def layout = 'full'
        def tables = render.call('full', 80)
        if (live) {
            // Every database when the block (tables + this summary + progress lines) fits the window.
            def term = terminalSize()
            tables = render.call('dense', term[1])
            layout = tables.sum { t -> t.readLines().size() } + lines.size() + 8 <= term[0] ? 'dense' : 'compact'
            if (layout == 'compact') tables = render.call('compact', term[1])
        }
        if (tables) lines.add(0, tables.join('\n') + '\n')
        lines << "  ${c.dim}Reports      ${c.reset}${params.outdir}/database_check/<build>.db_audit.*.{tsv,json,txt}"
        if (layout == 'compact' && tables)
            lines << "  ${c.dim}Full table   ${c.reset}cat ${params.outdir}/database_check/<build>.db_audit.${report}.txt  ${c.dim}(window too small)${c.reset}"
        lines << rule(c.green)
        return '\n' + lines.join('\n')
    }
    lines << "  ${c.dim}Results      ${c.reset}${file(params.outdir)}"
    def totals = summaryFiles(ids, '_totals.tsv')
    if (totals) {
        def per_file = totals.collect { f -> f.readLines().findAll { l -> l.trim() }.collect { l -> l.split('\t').toList() } }
        def rows = [per_file[0][0]] + per_file.collectMany { r -> r.drop(1) }
        def header = rows[0]
        def wanted = [['file', 'file'], ['genome', 'build'], ['n_genotype_columns', 'samples'], ['total', 'records'], ['snv', 'SNV'],
                      ['indel', 'indel'], ['ti_tv', 'Ti/Tv'], ['pass', 'PASS'], ['exonic', 'exonic'],
                      ['coding_rare', 'rare coding'], ['clinvar_plp', 'P/LP'], ['vus_candidates', 'VUS']]
            .findAll { w -> header.contains(w[0]) }
        // The live display cannot redraw a block taller than the terminal: list 10 files there.
        def files = rows.drop(1)
        def shown = live ? files.take(10) : files
        def table = [wanted.collect { w -> w[1] }] +
            shown.collect { r -> wanted.collect { w -> groupDigits(r[header.indexOf(w[0])]) } }
        def widths = (0..<wanted.size()).collect { i -> table.collect { r -> r[i].size() }.max() }
        lines << ''
        table.eachWithIndex { r, n ->
            def cells = r.withIndex().collect { cell, i -> i < 2 ? cell.padRight(widths[i]) : cell.padLeft(widths[i]) }
            lines << (n == 0 ? "  ${c.dim}${cells.join('  ')}${c.reset}" : "  ${cells.join('  ')}")
        }
        if (files.size() > shown.size())
            lines << "  ${c.dim}… ${files.size() - shown.size()} more input files: ${params.outdir}/summary/annotation_summary_totals.tsv${c.reset}"
        lines << ''
    }
    lines.addAll(plpGenesLog(ids, live ? 10 : Integer.MAX_VALUE))
    lines << "  ${c.dim}Tables       ${c.reset}${params.outdir}/annotation/<sample>/<sample>.<build>_multianno.txt.gz"
    if (ids.size() == 1) {
        lines << "  ${c.dim}HTML report  ${c.reset}${params.outdir}/summary/${ids[0]}.annotation_summary.html"
    } else {
        lines << "  ${c.dim}HTML reports ${c.reset}${params.outdir}/summary/<sample>.annotation_summary.html"
        lines << "  ${c.dim}             ${c.reset}${params.outdir}/summary/annotation_summary.html  ${c.dim}(all files)${c.reset}"
    }
    lines << rule(c.green)
    return '\n' + lines.join('\n')
}

/*
 * Database audit table (--mode check_db / update_db) from annovar_db.py TSV.
 * layout: 'full' (two-column grid with headers), or for the live display (see
 * liveDisplay()) 'dense' (every database, as many columns as fit termWidth) or
 * 'compact' (only what needs attention).
 */
def dbAuditLog(String build, tsv, String stage, String layout = 'full', int termWidth = 80) {
    def c = colors()
    def rows = tsv.splitCsv(header: true, sep: '\t')
    def titles = [check: 'database check', plan: 'database update plan', result: 'database status after update']
    def w = [rows.collect { r -> r.database.size() }.max() ?: 8, 9]
    def todo = rows.findAll { r -> r.action != 'none' }
    // The update plan is printed while downloads are running: keep it short (only what
    // will change), since tall blocks get mangled by Nextflow's live progress display.
    // The full table follows once the update has finished.
    if (stage == 'plan') {
        def head = "  ${c.bold}update plan${c.reset} · ${c.cyan}${build}${c.reset} · ${rows.size() - todo.size()} of ${rows.size()} databases ok"
        if (!todo) return "${head} · ${c.green}nothing to do${c.reset}"
        return ([head + ', to do:'] + todo.collect { r ->
            "    ${c.bold}→ ${r.action}${c.reset}  ${r.database}  ${c.dim}(${r.status})${c.reset}"
        }).join('\n')
    }
    // Two columns (read down, then right) with a short status badge per database;
    // anything that is not plainly 'ok' is explained in full below the grid.
    def badge = { String s ->
        if (s == 'ok')                return ['✔ ok', c.green]
        if (s == 'ok (warning)')      return ['✔ ok*', c.yellow]
        if (s in ['missing', 'damaged']) return ["✘ ${s}".toString(), c.red]
        if (s.startsWith('outdated')) return ['↑ outdated', c.yellow]
        if (s.startsWith('newer'))    return ['↑ newer', c.yellow]
        if (s == 'custom build')      return ['● custom', c.cyan]
        return ["? ${s}".toString(), c.yellow]
    }
    def size = { r ->
        if (r.installed != 'True') return '-'
        def gb = r.size_gb.toDouble()
        return String.format(gb >= 100 ? '%.0f' : (gb >= 10 ? '%.1f' : '%.2f'), gb) + 'G'
    }
    def nw = Math.max(8, w[0])
    def cell = { String name, String sz, String status, String colour, boolean pad ->
        "${name.padRight(nw)} ${sz.padLeft(5)}  ${colour}${pad ? status.padRight(10) : status}${c.reset}"
    }
    def width = 4 + 2 * (nw + 18)
    def line = "${c.blue}${'─' * Math.max(72, width)}${c.reset}"
    def half = (rows.size() + 1).intdiv(2)

    def issues = rows.findAll { r -> r.status != 'ok' }
    def issueLine = { r ->
        def b = badge.call(r.status)
        def action = r.action != 'none' ? "  ${c.bold}→ ${r.action}${c.reset}" : ''
        // 'not installed' only repeats the 'missing' status
        def notes = r.notes && r.status != 'missing' ? "${c.dim} · ${r.notes}${c.reset}" : ''
        "  ${b[1]}${b[0].take(1)}${c.reset} ${r.database}: ${b[1]}${r.status}${c.reset}${action}${notes}"
    }
    def summary = todo
        ? "${c.yellow}${todo.size()} of ${rows.size()} need action${c.reset}" +
          (stage == 'check' ? " · run with ${c.bold}--mode update_db${c.reset} to fix" : '')
        : "${c.green}all ${rows.size()} databases ready${c.reset}"

    // Live display: a title line per build, then (dense) every database in as many
    // columns as the terminal is wide, or (compact) only the databases that need attention.
    if (layout != 'full') {
        def gb = rows.findAll { r -> r.installed == 'True' }.sum { r -> r.size_gb.toDouble() } ?: 0.0
        def total = gb >= 1000 ? String.format('%.2f TB', gb / 1000) : String.format('%.1f GB', gb)
        def out = ["  ${c.bold}${titles.get(stage)}${c.reset} · ${c.cyan}${build}${c.reset} · ${summary}  ${c.dim}(${total})${c.reset}"]
        if (layout == 'dense' && rows) {
            def cw = nw + 9
            def ncol = Math.max(1, Math.min(4, (termWidth - 4).intdiv(cw + 2)))
            def nrow = (rows.size() + ncol - 1).intdiv(ncol)
            (0..<nrow).each { i ->
                def cells = (0..<ncol).collect { j -> i + j * nrow }.findAll { k -> k < rows.size() }.collect { k -> rows[k] }
                out << '    ' + cells.withIndex().collect { r, n ->
                    def b = badge.call(r.status)
                    def text = "${r.database.padRight(nw)} ${size.call(r).padLeft(5)} ${b[1]}${b[0].take(1)}${c.reset}"
                    n < cells.size() - 1 ? text + '   ' : text
                }.join('')
            }
        }
        out.addAll(issues.take(10).collect { r -> '  ' + issueLine.call(r) })
        if (issues.size() > 10) out << "      ${c.dim}… ${issues.size() - 10} more${c.reset}"
        return out.join('\n')
    }

    def lines = ['', line, "  ${c.bold}${titles.get(stage)}${c.reset} · ${c.cyan}${build}${c.reset} · ${params.annovar_db}", line]
    if (rows) {
        def head = cell.call('database', 'size', 'status', c.dim, true)
        lines << "  ${c.dim}${head}  ${rows.size() > 1 ? cell.call('database', 'size', 'status', c.dim, false) : ''}${c.reset}"
    }
    (0..<half).each { i ->
        def cells = [rows[i]] + (i + half < rows.size() ? [rows[i + half]] : [])
        lines << '  ' + cells.withIndex().collect { r, n ->
            def b = badge.call(r.status)
            cell.call(r.database, size.call(r), b[0], b[1], n == 0 && cells.size() > 1)
        }.join('  ')
    }

    if (issues) {
        lines << line
        issues.each { r -> lines << issueLine.call(r) }
    }
    lines << line
    lines << "  ${summary}"
    return lines.join('\n')
}

/*
 * Parameter validation ------------------------------------------------------------
 */
def validateParameters() {
    def errors = []
    def mode = params.mode
    if (!(mode in ['annotate', 'check_db', 'update_db']))
        errors << "--mode must be annotate, check_db or update_db (got '${mode}')"
    if (mode == 'annotate' && !params.input) errors << '--input is required'
    if (!params.annovar_db)  errors << '--annovar_db is required (ANNOVAR humandb directory)'
    if (mode != 'check_db' && !params.annovar_src)
        errors << '--annovar_src is required (directory with table_annovar.pl)'
    if (params.clinvar_release != 'latest' && !(params.clinvar_release.toString() ==~ /\d{8}/))
        errors << "--clinvar_release must be 'latest' or YYYYMMDD (got '${params.clinvar_release}')"
    dbBuilds().findAll { b -> !(b in ['hg19', 'hg38']) }.each { b -> errors << "--db_builds: unsupported build '${b}'" }

    // A new humandb is created when databases are going to be installed into it.
    if (params.annovar_db && !file(params.annovar_db).exists() && (mode == 'update_db' || autoInstall()))
        file(params.annovar_db).mkdirs()
    if (params.annovar_db && !file(params.annovar_db).isDirectory())
        errors << "--annovar_db is not a directory: ${params.annovar_db}"
    if (params.annovar_src && !file("${params.annovar_src}/table_annovar.pl").exists())
        errors << "table_annovar.pl not found in --annovar_src: ${params.annovar_src}"
    else if (params.annovar_src) {
        // table_annovar.pl (run with perl) starts these as programs: without the
        // execute bit ANNOVAR fails with 'Permission denied'.
        def scripts = ['annotate_variation.pl', 'convert2annovar.pl', 'coding_change.pl']
        def not_exec = scripts.findAll { s ->
            def f = file("${params.annovar_src}/${s}")
            f.exists() && !java.nio.file.Files.isExecutable(f)
        }
        if (not_exec)
            errors << "ANNOVAR scripts are not executable: ${not_exec.join(', ')}\n    fix: chmod +x ${params.annovar_src}/*.pl"
    }

    if (!(normalizeGenome(params.genome) in ['auto', 'hg19', 'hg38']))
        errors << "--genome must be auto, hg19 or hg38 (got '${params.genome}')"
    if (!(params.check_ref in ['e', 'w', 'x', 's', 'ws', 'xs', 'wx']))
        errors << "--check_ref must be a bcftools norm --check-ref mode (e, w, x, s; got '${params.check_ref}')"
    [chunk_size: params.chunk_size, threads: params.threads, min_ac: params.min_ac,
     db_download_forks: params.db_download_forks, family_max_samples: params.family_max_samples,
     qc_min_dp: params.qc_min_dp, qc_min_gq: params.qc_min_gq, qc_min_alt_reads: params.qc_min_alt_reads].each { name, value ->
        if (!(value.toString() ==~ /\d+/)) errors << "--${name} must be a non-negative integer (got '${value}')"
    }
    if (!(params.qc_het_vaf.toString() ==~ /0?\.\d+,(0?\.\d+|1(\.0+)?)/))
        errors << "--qc_het_vaf must be MIN,MAX fractions, e.g. 0.2,0.8 (got '${params.qc_het_vaf}')"
    if (!(params.qc_min_hom_vaf.toString() ==~ /0?\.\d+|1(\.0+)?/))
        errors << "--qc_min_hom_vaf must be a fraction, e.g. 0.8 (got '${params.qc_min_hom_vaf}')"
    if (params.ped && !file(params.ped).exists())
        errors << "--ped not found: ${params.ped}"
    genePanels().findAll { p -> !p.exists() }.each { p -> errors << "--gene_panel not found: ${p}" }
    if (params.threads.toString() == '0') errors << '--threads must be at least 1'
    if (params.operation && !params.protocol)
        errors << '--operation requires --protocol'
    if (params.protocol && params.operation &&
        params.protocol.toString().split(',').size() != params.operation.toString().split(',').size())
        errors << '--protocol and --operation must have the same number of comma-separated entries'
    ['hg19', 'hg38'].each { g ->
        def fasta = fastaParam(g)
        if (fasta && !file(fasta).exists()) errors << "--fasta_${g} not found: ${fasta}"
    }
    BOOLEAN_PARAMS().findAll { name -> !(params.get(name).toString() in ['true', 'false']) }.each { name ->
        errors << "--${name} must be true or false (got '${params.get(name)}')"
    }

    if (errors) {
        error("Invalid parameters:\n  - " + errors.join('\n  - ') + "\n\nRun with --help for usage.")
    }

    // Fail fast when the build is fixed up front; with 'auto' every detected
    // build is checked as soon as it is known (after installing what is missing).
    def genome = normalizeGenome(params.genome)
    if (mode == 'annotate' && genome != 'auto' && !autoInstall()) checkAnnovarResources(genome)
}

// Builds handled by --mode check_db/update_db: --db_builds, else every build with a protocol.
def dbBuilds() {
    if (params.db_builds) return params.db_builds.toString().split(',').collect { b -> normalizeGenome(b) }
    return ['hg19', 'hg38'].findAll { g -> params.protocol || params.annovar_protocols?.get(g) }
}

// Numeric parameters arrive as strings when given on the command line.
def intParam(value) {
    return value.toString().toInteger()
}

// So do boolean ones: '--onetranscript false' is the string 'false', which Groovy treats as true.
def flag(String name) {
    return params.get(name).toString() == 'true'
}

def BOOLEAN_PARAMS() {
    return ['recursive', 'skip_build_check', 'onetranscript', 'skip_db_check', 'auto_install_db', 'db_offline',
            'pass_only', 'drop_genotypes', 'skip_sample_qc', 'gene_info', 'compress_output', 'save_vcf',
            'save_normalized', 'monochrome_logs']
}

// The report's gene table in humandb (built by 'annovar_db.py genes', DB_GENES).
def geneTable() {
    return file("${params.annovar_db}/annovar_nf/genes.tsv")
}

// --gene_panel files (comma-separated).
def genePanels() {
    return params.gene_panel ? params.gene_panel.toString().split(',').collect { p -> file(p.trim()) } : []
}

// Banner text for the gene table: when it was built, or what happens without it.
def geneTableSummary() {
    if (!flag('gene_info')) return 'off (--gene_info false)'
    def table = geneTable()
    if (!table.exists()) {
        return autoInstall() || params.mode == 'update_db' ? 'GenCC + gnomAD constraint: built in humandb on this run'
            : 'not installed - the report shows no gene information (--mode update_db builds it)'
    }
    def built = table.withReader { r -> r.readLine() }?.find(/built (\d{4}-\d{2}-\d{2})/) { _m, d -> d }
    return "GenCC + gnomAD constraint (annovar_nf/genes.tsv${built ? ", built ${built}" : ''})"
}

// CPUs ANNOVAR gets per input file: --threads capped by resourceLimits (--max_cpus).
def annovarCpus() {
    return Math.min(intParam(params.threads), intParam(params.max_cpus))
}

// CPUs of the chunks of a file: its CPUs shared in proportion to records
// (largest remainder), at least 1 each, so that all chunks run at the same
// time. PREPARE_VCF plans one chunk per CPU, so most chunks get 1 CPU and run
// ANNOVAR as a single process; a file split into fewer chunks than CPUs gives
// the larger chunks ANNOVAR threads, at most one per 1,000 records (every
// thread loads the database indexes, about a minute).
def chunkCpus(List<Long> records) {
    def cpus = annovarCpus()
    def total = Math.max(records.sum() as long, 1L)
    def exact = records.collect { n -> cpus * n / total }
    def alloc = exact.collect { x -> Math.max(1, Math.floor(x) as int) }
    def spare = cpus - (alloc.sum() as int)
    def order = (0..<records.size()).toList().sort { i -> -(exact[i] - Math.floor(exact[i])) }
    order.take(Math.max(spare, 0)).each { i -> alloc[i] += 1 }
    return alloc.withIndex().collect { n, i -> Math.min(n, Math.max(1L, records[i].intdiv(1000)) as int) }
}

// Chunks PREPARE_VCF plans per file: one per CPU (plan_chunks.py --parts).
def chunkParts() {
    return annovarCpus()
}

/*
 * Genome helpers ------------------------------------------------------------------
 */
def normalizeGenome(value) {
    def g = (value ?: 'auto').toString().trim().toLowerCase()
    def aliases = [grch37: 'hg19', b37: 'hg19', hs37d5: 'hg19', grch38: 'hg38', hs38dh: 'hg38', '': 'auto', 'null': 'auto']
    return aliases.containsKey(g) ? aliases[g] : g
}

def fastaParam(String genome) {
    return genome == 'hg19' ? params.fasta_hg19 : (genome == 'hg38' ? params.fasta_hg38 : null)
}

def referenceFasta(String genome) {
    def path = fastaParam(genome)
    return path ? file(path, checkIfExists: true) : []
}

def referenceFai(String genome) {
    def path = fastaParam(genome)
    if (!path) return []
    def fai = file("${path}.fai")
    return fai.exists() ? fai : []
}

/*
 * ANNOVAR protocol --------------------------------------------------------------
 * Entries are '<database>:<operation>'. '<family>_latest' (e.g. clinvar_latest)
 * resolves to the newest installed '<family>_<YYYYMMDD>' database at run time.
 */
def protocolPairs(String genome) {
    if (params.protocol) {
        def dbs = params.protocol.toString().split(',').collect { d -> d.trim() }
        def ops = params.operation
            ? params.operation.toString().split(',').collect { o -> o.trim() }
            : dbs.collect { d -> d ==~ GENE_DB_RE() ? 'g' : 'f' }
        return [dbs, ops].transpose()
    }
    def entries = params.annovar_protocols?.get(genome)
    if (!entries) error("No default ANNOVAR protocol for '${genome}'. Set --protocol or params.annovar_protocols.${genome}.")
    return entries.collect { e ->
        def parts = e.toString().split(':')
        [parts[0].trim(), parts.size() > 1 ? parts[1].trim() : (parts[0].trim() ==~ GENE_DB_RE() ? 'g' : 'f')]
    }
}

// 'db:op,db:op,...' exactly as configured (aliases unresolved), for annovar_db.py.
def protocolEntries(String genome) {
    return protocolPairs(genome).collect { p -> "${p[0]}:${p[1]}" }.join(',')
}

// Newest installed version for '<family>_latest'; the name itself otherwise; null if none.
def resolveDatabase(String genome, String name) {
    if (!name.endsWith('_latest')) return name
    def family = name[0..<(name.size() - '_latest'.size())]
    def dir = file(params.annovar_db)
    if (!dir.isDirectory()) return null
    // <genome>_<family>_<digits>.txt, e.g. hg38_clinvar_20261025.txt
    def prefix = "${genome}_${family}_".toString()
    def versions = dir.listFiles()
        .collect { f -> f.name }
        .findAll { n -> n.startsWith(prefix) && n.endsWith('.txt') && n.substring(prefix.size(), n.size() - 4).isLong() }
        .collect { n -> n.substring(genome.size() + 1, n.size() - 4) }
        .sort { n -> n.substring(family.size() + 1).toLong() }
    return versions ? versions.last() : null
}

/*
 * Returns [protocol, operation] strings for table_annovar.pl (aliases resolved).
 */
def annovarProtocol(String genome) {
    def pairs = protocolPairs(genome).collect { p ->
        def db = resolveDatabase(genome, p[0])
        if (db == null) error("${p[0]}: no installed ${genome} version in ${params.annovar_db}. Run with --mode update_db.")
        [db, p[1]]
    }
    return [pairs.collect { p -> p[0] }.join(','), pairs.collect { p -> p[1] }.join(',')]
}

/*
 * The ANNOVAR databases a run uses for one build: protocol entry, resolved name
 * (version), operation, file, file date and size in bytes (-1: missing).
 * The file date is the ANNOVAR release date for files fetched by
 * annotate_variation.pl, otherwise the date the file was built/installed.
 */
def databaseInfo(String genome) {
    return protocolPairs(genome).collect { p ->
        def db = resolveDatabase(genome, p[0]) ?: p[0]
        def txt = file("${params.annovar_db}/${genome}_${db}.txt")
        [protocol: p[0], database: db, operation: p[1], file: txt.name,
         date: txt.exists() ? new Date(txt.lastModified()).format('yyyy-MM-dd') : 'missing',
         bytes: txt.exists() ? txt.size() : -1]
    }
}

/*
 * software_versions.yml block listing the ANNOVAR databases of one build.
 */
def databaseVersions(String genome) {
    def lines = ["\"ANNOVAR databases (${genome})\":", "    humandb: \"${params.annovar_db}\""]
    databaseInfo(genome).each { d ->
        def fields = []
        if (d.database != d.protocol) fields << "protocol: \"${d.protocol}\""
        fields << "operation: \"${d.operation}\""
        fields << "file: \"${d.file}\""
        fields << "file_date: \"${d.date}\""
        if (d.bytes >= 0) {
            def bytes = d.bytes
            def size = bytes >= 1024 ** 3 ? String.format('%.1f GB', bytes / (1024 ** 3))
                : bytes >= 1024 ** 2 ? String.format('%.0f MB', bytes / (1024 ** 2))
                : String.format('%.0f KB', bytes / 1024)
            fields << "size: \"${size}\""
        }
        lines << "    \"${d.database}\": {${fields.join(', ')}}"
    }
    return lines.join('\n')
}

/*
 * The same as TSV rows for the summary report (columns: see summarize_annotations.py).
 */
def databaseTable(String genome) {
    return databaseInfo(genome).collect { d ->
        [genome, d.protocol, d.database, d.operation, d.file, d.date, d.bytes].join('\t')
    }.join('\n')
}

/*
 * Verify the protocol's databases before any annotation starts (ANNOVAR would
 * otherwise fail - or silently skip an outdated index - hours into a run):
 * table present, *Mrna.fa for gene-based databases, .idx consistent with the
 * table size. Deep checks: --mode check_db.
 */
// Protocol databases of a build without an installed table.
def missingDatabases(String genome) {
    return protocolPairs(genome).findAll { p ->
        def db = resolveDatabase(genome, p[0])
        !db || !file("${params.annovar_db}/${genome}_${db}.txt").exists()
    }.collect { p -> p[0] }.unique()
}

// Annotation installs missing protocol databases before it starts (--auto_install_db).
def autoInstall() {
    return params.mode == 'annotate' && flag('auto_install_db') && !flag('skip_db_check')
}

def checkAnnovarResources(String genome) {
    if (flag('skip_db_check')) return
    def problems = []
    def missing = missingDatabases(genome)
    protocolPairs(genome).findAll { p -> !(p[0] in missing) }.each { p ->
        def db = resolveDatabase(genome, p[0])
        def txt = file("${params.annovar_db}/${genome}_${db}.txt")
        if (p[1].startsWith('g') && !file("${params.annovar_db}/${genome}_${db}Mrna.fa").exists())
            problems << "${genome}_${db}Mrna.fa is missing (gene-based database ${db})"
        def idx = file("${txt}.idx")
        if (idx.exists()) {
            def header = idx.withReader { r -> r.readLine() } ?: ''
            def expected = header.find(/^#BIN\t\d+\t(\d+)/) { _m, size -> size }
            if (expected && expected.toLong() != txt.size())
                problems << "${genome}_${db}.txt.idx does not match the table size (incomplete download or modified file)"
        }
    }
    if (missing) {
        def available = file(params.annovar_db).listFiles()
            .collect { f -> f.name }
            .findAll { n -> n.startsWith("${genome}_") && n.endsWith('.txt') }
            .collect { n -> n.replaceFirst("^${genome}_", '').replaceFirst(/\.txt$/, '') }
            .findAll { n -> missing.any { m -> n.toLowerCase().startsWith(m.toLowerCase().replaceAll(/(_latest|[\d_]+)$/, '')) } }
        problems.add(0, "missing: ${missing.join(', ')}" + (available ? "\n    similar databases installed: ${available.sort().join(', ')}" : ''))
    }
    if (problems) {
        error("""
            ANNOVAR databases for ${genome} in ${params.annovar_db} are not ready:
              - ${problems.join('\n              - ')}

            Fix:  nextflow run main.nf --mode update_db [your profile/params]   (installs / repairs / builds ClinVar)
            Check: nextflow run main.nf --mode check_db  (read-only report)
            Or adjust params.annovar_protocols.${genome} / --protocol, or use --skip_db_check.""".stripIndent())
    }
}

/*
 * Input discovery -----------------------------------------------------------------
 * Returns a list of [meta, vcf, index] where index is [] when absent.
 */
def sampleIdFromName(String name) {
    def base = name.replaceFirst(VARIANT_FILE_RE(), '')
    // Strip stacked suffixes (e.g. .hard-filtered.FINAL); 5 passes is plenty.
    def id = (1..5).inject(base) { acc, _pass -> acc.replaceFirst(SAMPLE_SUFFIX_RE(), '') }
    return id.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

def findIndex(vcf) {
    def idx = ['tbi', 'csi'].collect { ext -> file("${vcf}.${ext}") }.find { f -> f.exists() }
    return idx ?: []
}

def readInputs(String input) {
    def entries = []
    def lower = input.toLowerCase()

    if (lower.endsWith('.csv') || lower.endsWith('.tsv')) {
        def sheet = file(input, checkIfExists: true)
        def rows = sheet.splitCsv(header: true, sep: lower.endsWith('.tsv') ? '\t' : ',', strip: true)
        if (!rows) error("Samplesheet is empty: ${input}")
        def columns = rows[0].keySet()
        if (!columns.contains('vcf')) error("Samplesheet ${input} must have a 'vcf' column (found: ${columns.join(', ')})")
        rows.eachWithIndex { row, i ->
            def raw = row.vcf?.toString()?.trim()
            if (!raw) error("Samplesheet ${input}, row ${i + 2}: empty 'vcf' value")
            // Relative paths are resolved against the samplesheet location.
            def vcf = (raw.startsWith('/') || raw =~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\//) ? file(raw) : file(sheet.parent.resolve(raw).toString())
            if (!vcf.exists()) error("Samplesheet ${input}, row ${i + 2}: file not found: ${vcf}")
            def id = row.sample ? row.sample.toString().replaceAll(/[^A-Za-z0-9._-]/, '_') : sampleIdFromName(vcf.name)
            entries << [[id: id, genome: genomeOrNull(row.genome ?: params.genome), input: vcf.name], vcf, findIndex(vcf)]
        }
    } else {
        def target = file(input)
        def matches
        if (!(target instanceof List) && target.isDirectory()) {
            def prefix = flag('recursive') ? '**/' : ''
            matches = files("${target}/${prefix}${params.input_pattern}")
        } else {
            matches = files(input)
        }
        matches = matches.findAll { f -> f.name ==~ /(?i).*\.(vcf|gvcf|bcf)(\.b?gz)?$/ }.sort { f -> f.toString() }
        if (!matches) error("No VCF/BCF files found for --input '${input}' (pattern: ${params.input_pattern})")
        matches.each { vcf ->
            entries << [[id: sampleIdFromName(vcf.name), genome: genomeOrNull(params.genome), input: vcf.name], vcf, findIndex(vcf)]
        }
    }

    def duplicated = entries.countBy { e -> e[0].id }.findAll { _id, n -> n > 1 }.keySet()
    if (duplicated) {
        error("Duplicate sample ids derived from --input: ${duplicated.join(', ')}. " +
              "Use a samplesheet with an explicit 'sample' column to disambiguate.")
    }
    return entries
}

// 'auto' becomes null so PREPARE_VCF runs auto-detection for that file.
def genomeOrNull(value) {
    def g = normalizeGenome(value)
    if (!(g in ['auto', 'hg19', 'hg38'])) error("Unsupported genome '${value}' (use hg19, hg38 or auto)")
    return g == 'auto' ? null : g
}
