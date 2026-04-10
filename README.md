# annovar-nf

**A Nextflow pipeline for scalable, reproducible multi-sample ANNOVAR variant annotation (hg19 / hg38)**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

`annovar-nf` wraps [ANNOVAR](https://annovar.openbioinformatics.org/) into a portable Nextflow DSL2 pipeline that:

- Annotates **multiple VCF samples in parallel**
- Supports both **hg19** and **hg38** with curated, up-to-date database sets
- Automatically **replaces the generic `Otherinfo` header** in ANNOVAR TXT output with real VCF column names
- Generates an **HTML + TSV summary** across all samples
- Produces Nextflow **timeline, report, trace, and DAG** for full provenance tracking
- Supports **multi-threading** per sample via `--threads`
- Runs locally, on **SLURM/LSF clusters**, or with **Docker/Singularity**

---

## Pipeline Overview

```
Input VCFs (directory or samplesheet CSV)
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  ANNOVAR_TABLE                                      │
│  table_annovar.pl                                   │
│  hg38: 19 databases (gnomAD 4.1, ClinVar 2025,    │
│         InterVar 2025, dbNSFP 4.7a, REVEL…)        │
│  hg19: 40 databases (legacy extended set)           │
│  Output: .multianno.vcf + .multianno.txt            │
└─────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  FIX_HEADER                                         │
│  Replace Otherinfo columns with real VCF            │
│  sample column names (#CHROM POS … FORMAT SAMPLE)  │
│  Output: .annotated.header.txt                      │
└─────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  MULTIQC_SUMMARY (all samples collected)            │
│  Per-sample counts: total, functional classes,      │
│  ClinVar pathogenic, REVEL ≥ 0.75, CADD ≥ 20       │
│  Output: annotation_summary.tsv + .html             │
└─────────────────────────────────────────────────────┘
```

---

## Requirements

| Tool | Version |
|------|---------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04.0 |
| [ANNOVAR](https://annovar.openbioinformatics.org/) | 2019-10-24 or later |
| Perl | ≥ 5.26 |
| Python | ≥ 3.10 |

> ANNOVAR requires a free academic license. Download from https://annovar.openbioinformatics.org/

---

## Installation

```bash
git clone https://github.com/adaniyarov/annovar-nf.git
cd annovar-nf
chmod +x bin/*.py
nextflow -version   # confirm ≥ 23.04
```

---

## Quick Start

### Option 1 — Directory of VCFs

```bash
nextflow run main.nf \
    --input      /data/vcfs/ \
    --genome     hg38 \
    --annovar_db  /data/PublicData/annovar_src/annovar_20190101/humandb \
    --annovar_src /data/PublicData/annovar_src/annovar_20190101 \
    --threads    10 \
    --outdir     results/
```

### Option 2 — Samplesheet CSV

```csv
sample,vcf
test,/data/vcfs/test.vcf
```

```bash
nextflow run main.nf \
    --input      samplesheet.csv \
    --genome     hg38 \
    --annovar_db  /data/PublicData/annovar_src/annovar_20190101/humandb \
    --annovar_src /data/PublicData/annovar_src/annovar_20190101 \
    --threads    10 \
    --outdir     results/
```

### Resume a failed run

```bash
nextflow run main.nf [same options] -resume
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | VCF directory or samplesheet CSV | *required* |
| `--genome` | Reference genome: `hg19` or `hg38` | `hg38` |
| `--annovar_db` | Path to ANNOVAR `humandb/` | *required* |
| `--annovar_src` | Path to directory with `table_annovar.pl` | *required* |
| `--threads` | Threads per sample (`--thread` in ANNOVAR) | `4` |
| `--outdir` | Output directory | `./results` |
| `--protocol` | Custom protocol string (overrides defaults) | auto |
| `--operation` | Custom operation string (overrides defaults) | auto |

---

## Default Databases

### hg38 — 19 databases (production-validated)

| Category | Databases |
|----------|-----------|
| Gene annotation | `refGene`, `knownGene`, `ensGene` |
| dbSNP | `avsnp151` |
| Population freq | `1000G 2015` (ALL, AFR, SAS, EUR, EAS) |
| Somatic | `cosmic70` |
| Exome freq | `exac03` |
| Population freq | `gnomad41_exome`, `gnomad41_genome` |
| Pathogenicity | `revel`, `nci60` |
| Clinical | `clinvar_20250721`, `intervar_20250721` |
| Functional pred | `dbnsfp47a_interpro` |
| Regulatory | `regsnpintron` |

### hg19 — 40 databases (legacy extended set)

All hg38 databases (older versions) plus: `snp138`, `avsnp138`, `ESP6500`, `popfreq_all`, `kaviar`, `hrcr1`, `abraom`, `GME`, `cg69`, `kgXref`, `ICGC 21`, `cosmic68`, `ClinVar 20180603`, `mitimpact24`, `GERP++`, `CADD 1.3`, `FATHMM`, `GWAVA`, `Eigen`

---

## Output Structure

```
results/
├── annotated/
│   └── test/
│       ├── test.hg38_multianno.vcf
│       ├── test.hg38_multianno.txt
│       └── test.annovar.log
├── final/
│   └── test/
│       └── test.annotated.header.txt   ← main output for downstream analysis
├── summary/
│   ├── annotation_summary.tsv
│   └── annotation_summary.html
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.html
```

---

## Sample Name Extraction

When using `--input` with a directory, sample names are derived automatically by stripping common suffixes:

| File name | Extracted sample |
|-----------|-----------------|
| `test.hard-filtered.vcf` | `test` |
| `SAMPLE_001.FINAL.vcf` | `SAMPLE_001` |
| `SAMPLE_002.filtered.vcf` | `SAMPLE_002` |

Use a samplesheet CSV for full control over sample naming.

---

## Running on SLURM

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --annovar_db /data/PublicData/annovar_src/annovar_20190101/humandb \
    --annovar_src /data/PublicData/annovar_src/annovar_20190101 \
    --threads 10 \
    --outdir results/ \
    -profile slurm
```

Adjust queue name and account in `nextflow.config`.

---

## Citation

If you use `annovar-nf` in your research, please cite:
> Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. *Nucleic Acids Research*. 2010;38(16):e164.

---

## License

MIT © Asset Daniyarov
