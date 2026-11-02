/*
================================================================================
    GATHER_CHUNKS
    Merge the per-chunk ANNOVAR tables of one input file (chunks arrive in
    genomic order): header once + data rows; BGZF members are simply
    concatenated, so nothing is decompressed or recompressed.
================================================================================
*/

process GATHER_CHUNKS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(headers, stageAs: 'headers/*'), path(bodies, stageAs: 'bodies/*')

    output:
    tuple val(meta), path("${meta.id}.${meta.genome}_multianno.txt{,.gz}"),  emit: txt
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.${meta.genome}_multianno"
    // Keep the order chosen upstream (genomic chunk order). A single staged
    // file arrives as a Path, not a List, hence the normalisation.
    // Staged values already carry their stageAs directory (headers/..., ...).
    def header_list = [headers].flatten().join(' ')
    def body_list   = [bodies].flatten().join(' ')
    """
    set -euo pipefail

    # ---- annotated table -------------------------------------------------------
    header=""
    for h in ${header_list}; do
        if [ -s "\$h" ]; then header="\$h"; break; fi
    done
    if [ -z "\$header" ]; then
        echo "[GATHER_CHUNKS] WARNING: ${meta.id} has no annotated variants" >&2
        : > empty_header.tsv; header=empty_header.tsv
    fi

    if [ "${params.compress_output}" = "true" ]; then
        bgzip -c "\$header" > ${prefix}.txt.gz
        cat ${body_list} >> ${prefix}.txt.gz
    else
        cat "\$header" > ${prefix}.txt
        for b in ${body_list}; do gzip -dc "\$b" >> ${prefix}.txt; done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(bgzip --version | sed -n '1s/^bgzip (htslib) //p')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.genome}_multianno"
    def txt    = params.compress_output.toString() == 'true' ? "${prefix}.txt.gz" : "${prefix}.txt"
    """
    touch ${txt}
    printf '"${task.process}":\\n    htslib: stub\\n' > versions.yml
    """
}
