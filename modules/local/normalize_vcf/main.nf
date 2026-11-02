/*
================================================================================
    NORMALIZE_VCF
    Extract one genomic chunk and make it ANNOVAR-friendly:
      - split multi-allelic records into bi-allelic ones (bcftools norm -m -any)
      - optionally left-align / check REF against the reference FASTA
      - drop records ANNOVAR cannot annotate: spanning deletions (ALT='*'),
        symbolic alleles (<NON_REF>, <DEL>, <*> ...) and ALT='.' (gVCF blocks)
      - optionally keep FILTER=PASS only and/or drop genotype columns
================================================================================
*/

process NORMALIZE_VCF {
    tag "${meta.id}:${meta.chunk}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(index), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.${meta.chunk}.norm.vcf.gz"), emit: vcf
    path "versions.yml",                                           emit: versions

    script:
    def prefix  = "${meta.id}.${meta.chunk}"
    // --regions-overlap pos: a record belongs to the chunk that contains its
    // POS, so indels spanning a window boundary are never duplicated.
    // Regions go through a file: contig names may contain ':' and '*'
    // (HLA-DRB1*03:01:01:01), which bcftools cannot parse in --regions.
    // VCF contig names cannot contain quotes, so single-quoting is safe.
    def has_regions  = meta.regions && meta.regions != '.'
    def region_args  = has_regions ? meta.regions.tokenize(',').collect { item ->
        def i = item.lastIndexOf(':')
        def coords = item.substring(i + 1).tokenize('-')
        "'${item.substring(0, i)}' ${coords[0]} ${coords[1]}"
    }.join(' ') : ''
    def write_regions = has_regions ? "printf '%s\\t%s\\t%s\\n' ${region_args} > regions.tsv" : ''
    def regions = has_regions ? '--regions-file regions.tsv --regions-overlap pos' : ''
    def filters = params.pass_only.toString() == 'true' ? "--apply-filters 'PASS,.'" : ''
    def samples = params.drop_genotypes.toString() == 'true' ? '--drop-genotypes' : ''
    def ref     = fasta ? "--fasta-ref ${fasta} --check-ref ${params.check_ref}" : ''
    def min_ac  = params.drop_genotypes.toString() == 'true' ? '' : (params.min_ac.toString().toInteger() > 0 ? "--min-ac ${params.min_ac}" : '')
    """
    set -euo pipefail

    ${write_regions}

    bcftools view --no-version ${regions} ${filters} ${samples} -Ou "${vcf}" \\
        | bcftools norm --no-version --multiallelics -any ${ref} -Ou \\
        | bcftools view --no-version \\
            --exclude 'ALT="*" || ALT~"^<" || ALT="."' \\
            ${min_ac} \\
            -Oz -o ${prefix}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed -n '1s/^bcftools //p')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}.${meta.chunk}"
    """
    echo "" | gzip > ${prefix}.norm.vcf.gz
    printf '"${task.process}":\\n    bcftools: stub\\n' > versions.yml
    """
}
