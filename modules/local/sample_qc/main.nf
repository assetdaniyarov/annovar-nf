/*
================================================================================
    SAMPLE_QC
    Per-sample quality control of one input file (PASS records): counts from
    bcftools stats, sex from chrX outside the pseudo-autosomal regions and,
    for files with few samples, pairwise relatedness and trio Mendelian errors.
    Pedigree: --ped, else the VCF's ##PEDIGREE header lines, else inferred
    from the relatedness.
================================================================================
*/

process SAMPLE_QC {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(index), val(build)
    path ped
    val  code               // checksum of the bin/ scripts: -resume reruns after they change

    output:
    path "${meta.id}.sample_qc.json", emit: json
    path "versions.yml",              emit: versions

    script:
    // chrX outside the pseudo-autosomal regions (PAR1 / PAR2).
    def non_par = build == 'hg19' ? '2699521-154931043' : '2781480-155701382'
    def ped_arg = ped ? "--ped ${ped}" : ''
    """
    # code ${code}
    set -euo pipefail

    N=\$(bcftools query -l "${vcf}" | wc -l)
    # Header to a file ('grep -q' on a pipe would fail it under pipefail); its
    # ##PEDIGREE lines are the pedigree when --ped does not cover the samples.
    bcftools view -h "${vcf}" > header.txt
    : > all.stats
    : > chrx.stats
    if [ "\$N" -gt 0 ]; then        # sites-only VCFs have no samples to check
        bcftools stats --threads ${task.cpus} -s - -f PASS,. "${vcf}" > all.stats
        X=\$( (tabix -l "${vcf}" 2> /dev/null || bcftools view -h "${vcf}" | sed -n 's/^##contig=<ID=\\([^,>]*\\).*/\\1/p') \\
            | grep -m1 -xE 'chrX|X' || true)
        if [ -n "\$X" ]; then
            bcftools stats -s - -f PASS,. -r "\$X:${non_par}" "${vcf}" > chrx.stats
        fi
    fi

    # Relatedness / Mendelian errors: files with a few samples only.
    if [ "\$N" -ge 2 ] && [ "\$N" -le ${params.family_max_samples} ]; then
        FMT='%GT'
        for tag in GQ DP; do
            if grep -q "^##FORMAT=<ID=\$tag," header.txt; then FMT="\$FMT:%\$tag"; else FMT="\$FMT:."; fi
        done
        bcftools view --threads ${task.cpus} -m2 -M2 -v snps -f PASS,. -Ou "${vcf}" \\
            | bcftools query -f "%CHROM[\\t\$FMT]\\n" \\
            | sample_qc.py --id "${meta.id}" --stats all.stats --chrx chrx.stats --genotypes - ${ped_arg} \\
                --vcf-header header.txt --out ${meta.id}.sample_qc.json
    else
        sample_qc.py --id "${meta.id}" --stats all.stats --chrx chrx.stats ${ped_arg} --vcf-header header.txt \\
            --out ${meta.id}.sample_qc.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed -n '1s/^bcftools //p')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"id": "${meta.id}", "samples": [], "pairs": [], "trios": []}' > ${meta.id}.sample_qc.json
    printf '"${task.process}":\\n    bcftools: stub\\n' > versions.yml
    """
}
