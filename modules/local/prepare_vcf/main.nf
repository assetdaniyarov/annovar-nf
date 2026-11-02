/*
================================================================================
    PREPARE_VCF
    - detect / validate the genome build from the VCF header
    - make sure the VCF is BGZF-compressed, coordinate-sorted and indexed
      (plain .vcf, gzip .vcf.gz, BGZF without index, unsorted lift-over
      output and BCF are all accepted)
    - plan genomic chunks for parallel annotation from the index statistics
================================================================================
*/

process PREPARE_VCF {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(index)
    val code                // checksum of the bin/ scripts: -resume reruns after they change
    val parts               // CPUs per input file: one chunk each

    output:
    tuple val(meta), path("${meta.id}.prepared.vcf.gz"), path("${meta.id}.prepared.vcf.gz.{tbi,csi}"), env('GENOME_BUILD'), env('GENOME_SOURCE'), path("${meta.id}.chunks.tsv"), emit: vcf
    path "${meta.id}.build.json", emit: build_report
    path "versions.yml",          emit: versions

    script:
    def prefix      = meta.id
    def expected    = meta.genome ?: 'auto'
    def strict      = params.skip_build_check.toString() == 'true' ? '--no-strict' : ''
    def sort_memory = task.memory ? "${(task.memory.toMega() * 0.6).intValue()}M" : '768M'
    def index_name  = index ? index.name : ''
    def index_ext   = index_name.endsWith('.csi') ? 'csi' : 'tbi'
    """
    # code ${code}
    set -euo pipefail

    # ---- 1. genome build (header only: fast, fails before any heavy work) ---
    bcftools view --header-only "${vcf}" > header.txt
    GENOME_BUILD=\$(detect_genome_build.py \\
        --header   header.txt \\
        --filename "${vcf}" \\
        --expected "${expected}" \\
        --report   ${prefix}.build.json \\
        ${strict})
    # How the build was determined (contig_lengths, explicit, header_hints, filename).
    GENOME_SOURCE=\$(sed -n 's/.*"source": *"\\([^"]*\\)".*/\\1/p' ${prefix}.build.json)

    # ---- 2. BGZF + sorted + indexed ------------------------------------------
    is_bgzf_vcf() {
        # BGZF magic: gzip member with 'BC' extra subfield at bytes 12-13,
        # and the decompressed stream must start with '##' (VCF, not BCF).
        [ "\$(head -c 14 "\$1" | tail -c 2)" = "BC" ] || return 1
        [ "\$( (gzip -dc "\$1" 2>/dev/null || true) | head -c 2)" = "##" ]
    }

    out=${prefix}.prepared.vcf.gz
    if [ -n "${index_name}" ] && is_bgzf_vcf "${vcf}"; then
        # Already analysis-ready: link instead of rewriting (31 GB cohorts!).
        ln -s "${vcf}" "\$out"
        ln -s "${index_name}" "\$out.${index_ext}"
    else
        bcftools view --no-version --threads ${task.cpus} -Oz -o unsorted.vcf.gz "${vcf}"
        if bcftools index --threads ${task.cpus} --tbi unsorted.vcf.gz 2> index.err; then
            mv unsorted.vcf.gz "\$out"
            mv unsorted.vcf.gz.tbi "\$out.tbi"
        else
            echo "[PREPARE_VCF] input is not coordinate-sorted (\$(head -c 200 index.err)); sorting" >&2
            mkdir -p sort_tmp
            bcftools sort --max-mem ${sort_memory} --temp-dir sort_tmp -Oz -o "\$out" unsorted.vcf.gz
            bcftools index --threads ${task.cpus} --tbi "\$out"
            rm -rf unsorted.vcf.gz sort_tmp
        fi
    fi

    # ---- 3. chunk plan ---------------------------------------------------------
    if ! bcftools index --stats "\$out" > index_stats.tsv 2> /dev/null || [ ! -s index_stats.tsv ]; then
        # Some third-party indexes carry no record counts: count them instead.
        bcftools query -f '%CHROM\\n' "\$out" | uniq -c | awk -v OFS='\\t' '{print \$2, ".", \$1}' > index_stats.tsv
    fi
    plan_chunks.py \\
        --stats      index_stats.tsv \\
        --header     header.txt \\
        --chunk-size ${params.chunk_size} \\
        --parts      ${parts} \\
        --out        ${prefix}.chunks.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed -n '1s/^bcftools //p')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    GENOME_BUILD=${meta.genome ?: 'hg38'}
    GENOME_SOURCE=stub
    echo "" | gzip > ${prefix}.prepared.vcf.gz
    touch ${prefix}.prepared.vcf.gz.tbi
    printf '0\\tall\\t.\\t1\\n' > ${prefix}.chunks.tsv
    echo '{}' > ${prefix}.build.json
    printf '"${task.process}":\\n    bcftools: stub\\n' > versions.yml
    """
}
