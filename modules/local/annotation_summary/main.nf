/*
================================================================================
    ANNOTATION_SUMMARY
    One report per input file (<id>.annotation_summary.{html,json,tsv,...})
    built from the per-chunk statistics, the genome-build reports, the
    per-sample QC (sex, relatedness, PED data), the list of ANNOVAR
    databases used, the gene table (GenCC, gnomAD constraint) and the gene
    panels; with several input files also a combined
    annotation_summary.* across all of them.
================================================================================
*/

process ANNOTATION_SUMMARY {
    label 'process_single'

    input:
    path stats,     stageAs: 'stats/*'
    path builds,    stageAs: 'builds/*'
    path sample_qc, stageAs: 'sample_qc/*'
    path databases
    path genes,     stageAs: 'genes/*'
    path panels,    stageAs: 'panels/?/*'     // one folder each: panels may share a file name
    val  annovar_version
    val  code               // checksum of the bin/ scripts: -resume reruns after they change

    output:
    path "*annotation_summary.tsv",            emit: matrix
    path "*annotation_summary_totals.tsv",     emit: totals
    path "*annotation_summary_samples.tsv",    emit: samples
    path "*annotation_summary_plp.tsv",        emit: plp
    path "*annotation_summary_vus.tsv",        emit: vus
    path "*annotation_summary_candidates.tsv", emit: candidates
    path "*annotation_summary.json",           emit: json
    path "*annotation_summary.html",           emit: html
    path "versions.yml",                       emit: versions

    script:
    def qc_args = sample_qc ? '--sample-qc sample_qc/*' : ''
    def gene_args = (genes ? '--genes genes/* ' : '') + (panels ? '--panel panels/*/*' : '')
    """
    # code ${code}
    summarize_annotations.py \\
        --stats   stats/* \\
        --builds  builds/* \\
        ${qc_args} ${gene_args} \\
        --databases ${databases} \\
        --annovar-version "${annovar_version}" \\
        --family-max-samples ${params.family_max_samples} \\
        --out     annotation_summary \\
        --version "${workflow.manifest.version}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    for f in stats/*; do
        id=\$(basename "\$f" | sed 's/\\.[^.]*\\.stats\\.json\$//')
        for s in .tsv _totals.tsv _samples.tsv _plp.tsv _vus.tsv _candidates.tsv .json .html; do
            touch \$id.annotation_summary\$s
        done
    done
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}
