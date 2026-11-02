/*
================================================================================
    ANNOVAR_TABLE
    Annotate one normalized chunk with table_annovar.pl and post-process the
    table in the same task (so the large raw table never leaves the work dir):
      - real VCF column names instead of Otherinfo4..N
      - ';' between genes in Gene.* / GeneDetail.* (table and VCF)
      - header and data rows written separately for cheap concatenation
      - per-chunk statistics for the summary report
================================================================================
*/

process ANNOVAR_TABLE {
    tag "${meta.id}:${meta.chunk}"
    label 'process_annovar'

    input:
    tuple val(meta), path(vcf), val(protocol), val(operation)
    path annovar_db
    path annovar_src
    val  code               // checksum of the bin/ scripts: -resume reruns after they change

    output:
    tuple val(meta), env('HAS_VARIANTS'), path("${meta.id}.${meta.chunk}.header.tsv"), path("${meta.id}.${meta.chunk}.body.tsv.gz"), emit: chunk
    // Annotated VCF only with --save_vcf.
    tuple val(meta), env('HAS_VARIANTS'), path("${meta.id}.${meta.chunk}.multianno.vcf.gz"), optional: true, emit: vcf
    path "${meta.id}.${meta.chunk}.stats.json",               emit: stats
    path "versions.yml",                                      emit: versions

    script:
    def prefix  = "${meta.id}.${meta.chunk}"
    def build   = meta.genome
    def threads = task.cpus > 1 ? "--thread ${task.cpus}" : ''
    def onetx   = params.onetranscript.toString() == 'true' ? '--onetranscript' : ''
    def args    = params.annovar_args ?: ''
    def save_vcf = params.save_vcf.toString() == 'true'
    def vcf_args = save_vcf ? "--annovar-vcf ${prefix}.${build}_multianno.vcf --vcf-out ${prefix}.multianno.vcf" : ''
    """
    # code ${code}
    set -euo pipefail

    # ANNOVAR needs an uncompressed VCF to build the *_multianno.vcf output.
    gzip -dc "${vcf}" > ${prefix}.input.vcf

    # table_annovar.pl runs annotate_variation.pl children that can outlive
    # the task: on Ctrl+C the Nextflow wrapper (.command.run) dies at once,
    # while this script and ANNOVAR ignore SIGINT (they run as background
    # jobs), so nothing ever tells them to stop. Run ANNOVAR in the
    # background, kill its whole process tree on TERM/HUP, and watch the
    # parent wrapper: when it disappears, stop ANNOVAR too.
    kill_tree() {
        local child
        for child in \$(pgrep -P "\$1" || true); do kill_tree "\$child"; done
        kill -TERM "\$1" 2>/dev/null || true
    }
    ANNOVAR_PID=''
    WATCHDOG_PID=''
    trap '[ -n "\$ANNOVAR_PID" ] && kill_tree "\$ANNOVAR_PID"; exit 143' TERM HUP

    if grep -q -v -m1 '^#' ${prefix}.input.vcf; then
        HAS_VARIANTS=true

        # Mitochondrial records get a twin under the other name of the
        # mitochondrion (chrM <-> MT): ANNOVAR matches names exactly, and its
        # databases use both. format_annovar_table.py merges the twins back.
        TWIN_ARG=''
        if [ "\$(format_annovar_table.py mito-twins ${prefix}.input.vcf --genome ${build})" -gt 0 ]; then
            TWIN_ARG='--mito-twins'
        fi

        perl ${annovar_src}/table_annovar.pl \\
            ${prefix}.input.vcf \\
            ${annovar_db} \\
            --buildver  ${build} \\
            --outfile   ${prefix} \\
            --protocol  ${protocol} \\
            --operation ${operation} \\
            --nastring  . \\
            --vcfinput \\
            --otherinfo \\
            --remove \\
            ${onetx} \\
            ${threads} \\
            ${args} \\
            > ${prefix}.annovar.log 2>&1 &
        ANNOVAR_PID=\$!
        (
            while kill -0 "\$PPID" 2>/dev/null && kill -0 "\$ANNOVAR_PID" 2>/dev/null; do sleep 5; done
            kill -0 "\$PPID" 2>/dev/null || kill_tree "\$ANNOVAR_PID"
        ) > /dev/null 2>&1 &
        WATCHDOG_PID=\$!
        ANNOVAR_STATUS=0
        wait "\$ANNOVAR_PID" || ANNOVAR_STATUS=\$?
        kill "\$WATCHDOG_PID" 2>/dev/null || true
        ANNOVAR_PID=''
        if [ "\$ANNOVAR_STATUS" -ne 0 ]; then
            echo "table_annovar.pl failed (exit \$ANNOVAR_STATUS); last lines of the log:" >&2
            tail -n 40 ${prefix}.annovar.log >&2
            exit 1
        fi

        format_annovar_table.py \\
            --txt        ${prefix}.${build}_multianno.txt \\
            --vcf        ${prefix}.input.vcf \\
            --header-out ${prefix}.header.tsv \\
            --stats      ${prefix}.stats.json \\
            --sample     "${meta.id}" \\
            --chunk      "${meta.chunk}" \\
            --genome     ${build} \\
            --qc-min-dp ${params.qc_min_dp} --qc-min-gq ${params.qc_min_gq} \\
            --qc-min-alt-reads ${params.qc_min_alt_reads} --qc-het-vaf ${params.qc_het_vaf} \\
            --qc-min-hom-vaf ${params.qc_min_hom_vaf} --family-max-samples ${params.family_max_samples} \\
            ${vcf_args} \$TWIN_ARG \\
        | bgzip --threads ${task.cpus} > ${prefix}.body.tsv.gz

        if [ "${save_vcf}" = "true" ]; then
            bgzip --threads ${task.cpus} -c ${prefix}.multianno.vcf > ${prefix}.multianno.vcf.gz
        fi
        rm -f ${prefix}.${build}_multianno.txt ${prefix}.${build}_multianno.vcf ${prefix}.multianno.vcf ${prefix}.avinput
    else
        # Empty chunk (e.g. everything filtered out): emit empty placeholders.
        HAS_VARIANTS=false
        : > ${prefix}.header.tsv
        printf '' | bgzip > ${prefix}.body.tsv.gz
        if [ "${save_vcf}" = "true" ]; then
            bgzip -c ${prefix}.input.vcf > ${prefix}.multianno.vcf.gz
        fi
        printf '{"sample": "%s", "chunk": "%s", "genome": "%s", "counts": {"total": 0}, "func": {}, "columns": {}}\\n' \\
            "${meta.id}" "${meta.chunk}" "${build}" > ${prefix}.stats.json
    fi
    rm -f ${prefix}.input.vcf

    # ANNOVAR release = newest script date (not every script changes per release).
    script_date() {
        local s
        for s in "\$@"; do grep -m1 -o 'Date: [0-9-]*' ${annovar_src}/\$s 2> /dev/null | sed 's/Date: //' || true; done
    }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$(script_date annotate_variation.pl table_annovar.pl coding_change.pl | sort | tail -n 1)
        table_annovar.pl: \$(script_date table_annovar.pl)
        annotate_variation.pl: \$(script_date annotate_variation.pl)
        perl: \$(perl -e 'print substr(\$^V, 1)')
        htslib: \$(bgzip --version | sed -n '1s/^bgzip (htslib) //p')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.chunk}"
    """
    HAS_VARIANTS=true
    printf 'Chr\\tStart\\tEnd\\tRef\\tAlt\\n' > ${prefix}.header.tsv
    echo "" | gzip > ${prefix}.body.tsv.gz
    if [ "${params.save_vcf}" = "true" ]; then echo "" | gzip > ${prefix}.multianno.vcf.gz; fi
    echo '{"sample": "${meta.id}", "chunk": "${meta.chunk}", "counts": {"total": 0}, "func": {}}' > ${prefix}.stats.json
    printf '"${task.process}":\\n    annovar: stub\\n' > versions.yml
    """
}
