/*
================================================================================
    GATHER_VCF  (only with --save_vcf)
    Merge the per-chunk ANNOVAR VCFs of one input file (chunks arrive in
    genomic order): bcftools concat --naive (no re-encoding) + tabix index.
================================================================================
*/

process GATHER_VCF {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcfs, stageAs: 'vcfs/*')

    output:
    tuple val(meta), path("${meta.id}.${meta.genome}_multianno.vcf.gz"), path("${meta.id}.${meta.genome}_multianno.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix   = "${meta.id}.${meta.genome}_multianno"
    // A single staged file arrives as a Path, not a List.
    def vcf_list = [vcfs].flatten().join(' ')
    """
    set -euo pipefail

    set -- ${vcf_list}
    if [ \$# -eq 1 ]; then
        ln -s "\$1" ${prefix}.vcf.gz
    else
        bcftools concat --no-version --naive --threads ${task.cpus} -Oz -o ${prefix}.vcf.gz "\$@"
    fi
    if ! bcftools index --threads ${task.cpus} --tbi ${prefix}.vcf.gz 2> index.err; then
        # Left-alignment can move a record across a chunk boundary: re-sort.
        echo "[GATHER_VCF] re-sorting concatenated VCF (\$(head -c 200 index.err))" >&2
        mkdir -p sort_tmp
        bcftools sort --temp-dir sort_tmp -Oz -o sorted.vcf.gz ${prefix}.vcf.gz
        rm -f ${prefix}.vcf.gz && mv sorted.vcf.gz ${prefix}.vcf.gz
        bcftools index --threads ${task.cpus} --tbi ${prefix}.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed -n '1s/^bcftools //p')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.genome}_multianno"
    """
    touch ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
    printf '"${task.process}":\\n    bcftools: stub\\n' > versions.yml
    """
}
