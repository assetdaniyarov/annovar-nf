#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    annovar-nf: Nextflow pipeline for ANNOVAR variant annotation
    Supports hg19 and hg38 genome builds
    Author: Asset Daniyarov
    License: MIT
================================================================================
*/

// ─── Parameter defaults ───────────────────────────────────────────────────────
params.input         = null          // Path to samplesheet CSV or directory with VCFs
params.outdir        = './results'
params.genome        = 'hg38'        // 'hg19' or 'hg38'
params.annovar_db    = null          // Path to ANNOVAR humandb/
params.annovar_src   = null          // Path to annovar_src/ (contains table_annovar.pl)
params.protocol      = null          // Optional: override protocol string
params.operation     = null          // Optional: override operation string
params.threads       = 4             // Threads for table_annovar.pl (--thread)
params.help          = false

// ─── Help message ─────────────────────────────────────────────────────────────
def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║              annovar-nf  v1.1.0                         ║
    ║  Nextflow pipeline for multi-sample ANNOVAR annotation  ║
    ╚══════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf \\
            --input      samplesheet.csv \\
            --genome     hg38 \\
            --annovar_db /path/to/humandb \\
            --annovar_src /path/to/annovar_src \\
            --outdir     results/

    Input options:
        --input         Samplesheet CSV (columns: sample,vcf) OR directory with VCF files
        --genome        Genome build: hg19 or hg38 [default: hg38]
        --annovar_db    Path to ANNOVAR humandb directory [required]
        --annovar_src   Path to directory containing table_annovar.pl [required]
        --outdir        Output directory [default: ./results]
        --threads       Threads per sample for table_annovar.pl [default: 4]

    Advanced options:
        --protocol      Custom ANNOVAR protocol string (overrides defaults)
        --operation     Custom ANNOVAR operation string (overrides defaults)

    Profiles:
        -profile standard    Local execution (default)
        -profile slurm       SLURM cluster
        -profile docker      Use Docker containers
        -profile singularity Use Singularity containers

    Examples:
        # Annotate all VCFs in a directory (hg38, 50 threads per sample)
        nextflow run main.nf \\
            --input      /data/vcfs/ \\
            --genome     hg38 \\
            --annovar_db  /data/PublicData/annovar_src/annovar_20190101/humandb \\
            --annovar_src /data/PublicData/annovar_src/annovar_20190101 \\
            --threads    50 \\
            --outdir     results/

        # Use a samplesheet
        nextflow run main.nf --input samplesheet.csv --genome hg19 \\
            --annovar_db /db/humandb --annovar_src /tools/annovar

        # Resume a failed run
        nextflow run main.nf [options] -resume
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

// ─── Validate inputs ──────────────────────────────────────────────────────────
if (!params.input)       error "ERROR: --input is required"
if (!params.annovar_db)  error "ERROR: --annovar_db is required"
if (!params.annovar_src) error "ERROR: --annovar_src is required"
if (!['hg19','hg38'].contains(params.genome)) error "ERROR: --genome must be 'hg19' or 'hg38'"

// ─── Default protocol/operation per genome build ──────────────────────────────
//
//  hg38 — based on real production run (EP001_21.hard-filtered):
//      19 databases: refGene, knownGene, ensGene + 16 filter dbs
//      gnomAD 4.1, ClinVar 20250721, InterVar 20250721, dbNSFP 4.7a
//
//  hg19 — extended legacy set (40 databases)
//
def getProtocol(genome) {
    if (params.protocol) return params.protocol
    def gene_dbs = "refGene,knownGene,ensGene"
    if (genome == 'hg38') {
        def filter_dbs = [
            "avsnp151",
            "ALL.sites.2015_08",
            "AFR.sites.2015_08",
            "SAS.sites.2015_08",
            "EUR.sites.2015_08",
            "EAS.sites.2015_08",
            "cosmic70",
            "exac03",
            "gnomad41_exome",
            "gnomad41_genome",
            "revel",
            "nci60",
            "clinvar_20250721",
            "intervar_20250721",
            "dbnsfp47a_interpro",
            "regsnpintron"
        ].join(",")
        return "${gene_dbs},${filter_dbs}"
    } else {
        // hg19 — extended legacy set
        def filter_dbs = [
            "snp138",
            "avsnp138",
            "avsnp150",
            "ALL.sites.2015_08","AFR.sites.2015_08","AMR.sites.2015_08",
            "SAS.sites.2015_08","EUR.sites.2015_08","EAS.sites.2015_08",
            "esp6500siv2_ea","esp6500siv2_all","esp6500siv2_aa",
            "popfreq_all_20150413",
            "abraom","hrcr1","kaviar_20150923",
            "cg69",
            "dbnsfp35a","dbscsnv11",
            "kgXref",
            "exac03nonpsych","exac03nontcga",
            "gnomad_exome","gnomad_genome",
            "gme","mcap","revel","nci60",
            "icgc21",
            "cosmic68","cosmic70",
            "clinvar_20180603",
            "mitimpact24",
            "regsnpintron",
            "gerp++elem","gerp++gt2",
            "cadd13",
            "fathmm","gwava","eigen"
        ].join(",")
        return "${gene_dbs},${filter_dbs}"
    }
}

def getOperation(genome) {
    if (params.operation) return params.operation
    def protocol = getProtocol(genome)
    def n_gene   = 3
    def n_filter = protocol.split(",").size() - n_gene
    return (["g"] * n_gene + ["f"] * n_filter).join(",")
}

// ─── Helper: build channel from directory or samplesheet ─────────────────────
def buildSampleChannel() {
    def input_path = file(params.input)
    if (input_path.isDirectory()) {
        return Channel
            .fromPath("${params.input}/*.vcf")
            .map { vcf ->
                // Strip common suffixes: .hard-filtered, .FINAL, .filtered, .snp_indel
                def sample = vcf.name
                    .replaceAll(/\.hard-filtered$/, '')
                    .replaceAll(/\.FINAL$/, '')
                    .replaceAll(/\.filtered$/, '')
                    .replaceAll(/\.snp_indel$/, '')
                    .replaceAll(/\.vcf$/, '')
                tuple(sample, vcf)
            }
    } else {
        // samplesheet CSV: columns sample,vcf
        return Channel
            .fromPath(params.input)
            .splitCsv(header: true, strip: true)
            .map { row -> tuple(row.sample, file(row.vcf)) }
    }
}

// ─── Workflow ─────────────────────────────────────────────────────────────────
workflow {

    log.info """
    ╔══════════════════════════════════════════════════╗
    ║              annovar-nf  v1.1.0                 ║
    ╚══════════════════════════════════════════════════╝
    input      : ${params.input}
    genome     : ${params.genome}
    annovar_db : ${params.annovar_db}
    annovar_src: ${params.annovar_src}
    threads    : ${params.threads}
    outdir     : ${params.outdir}
    """.stripIndent()

    ch_vcf = buildSampleChannel()

    ANNOVAR_TABLE(
        ch_vcf,
        params.genome,
        params.annovar_db,
        params.annovar_src,
        getProtocol(params.genome),
        getOperation(params.genome),
        params.threads
    )

    FIX_HEADER(
        ANNOVAR_TABLE.out.vcf_txt_pair
    )

    ANNOTATION_SUMMARY(
        FIX_HEADER.out.fixed_txt.collect()
    )
}

// ─── Process 1: Run table_annovar.pl ─────────────────────────────────────────
process ANNOVAR_TABLE {
    tag          "${sample}"
    label        'process_high'
    publishDir   "${params.outdir}/annotated/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(vcf)
    val   genome
    val   annovar_db
    val   annovar_src
    val   protocol
    val   operation
    val   threads

    output:
    tuple val(sample), path("${sample}.${genome}_multianno.vcf"), path("${sample}.${genome}_multianno.txt"), emit: vcf_txt_pair
    path  "*.log",                                                                                            emit: logs

    script:
    """
    perl ${annovar_src}/table_annovar.pl \\
        ${vcf} \\
        ${annovar_db} \\
        --protocol    ${protocol} \\
        --operation   ${operation} \\
        --outfile     ${sample} \\
        --buildver    ${genome} \\
        --thread      ${threads} \\
        --remove \\
        --otherinfo \\
        --onetranscript \\
        --nastring    '.' \\
        --vcfinput \\
        2>&1 | tee ${sample}.annovar.log
    """
}

// ─── Process 2: Fix VCF header in ANNOVAR TXT output ─────────────────────────
process FIX_HEADER {
    tag          "${sample}"
    label        'process_low'
    publishDir   "${params.outdir}/final/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(anno_vcf), path(anno_txt)

    output:
    tuple val(sample), path("${sample}.annotated.header.txt"), emit: vcf_txt_pair
    path  "${sample}.annotated.header.txt",                    emit: fixed_txt
    path  "${sample}.annotated.header.txt"

    script:
    """
    fix_vcf_header.py \\
        --vcf  ${anno_vcf} \\
        --txt  ${anno_txt} \\
        --out  ${sample}.annotated.header.txt
    """
}

// ─── Process 3: Summary across all samples ───────────────────────────────────
process ANNOTATION_SUMMARY {
    label        'process_low'
    publishDir   "${params.outdir}/summary", mode: 'copy', overwrite: true

    input:
    path(fixed_txts)

    output:
    path "annotation_summary.tsv"
    path "annotation_summary.html"
    path "annotation_summary_totals.tsv"

    script:
    def inputs = fixed_txts instanceof List ? fixed_txts.join(' ') : fixed_txts
    """
    summarize_annotations.py \\
        --inputs ${inputs} \\
        --out    annotation_summary
    """
}
