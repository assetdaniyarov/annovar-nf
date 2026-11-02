# annovar-nf

**Scalable, reproducible ANNOVAR annotation of VCF/BCF files (hg19 / hg38) with Nextflow.**

[![CI](https://github.com/assetdaniyarov/annovar-nf/actions/workflows/ci.yml/badge.svg)](https://github.com/assetdaniyarov/annovar-nf/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.04.6-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23880854.svg)](https://doi.org/10.5281/zenodo.23880854)

`annovar-nf` wraps [ANNOVAR](https://annovar.openbioinformatics.org/) `table_annovar.pl` in a Nextflow
DSL2 pipeline that accepts the VCFs you actually have — single-sample DRAGEN `*.hard-filtered.vcf.gz`,
joint-genotyped cohorts with hundreds of samples, lifted-over files, plain or compressed, sorted or not,
hg19 and hg38 mixed in one run — and turns them into consistently formatted annotation tables.

## Highlights

- **Any VCF-like input**: `.vcf`, `.vcf.gz` (BGZF *or* plain gzip), `.vcf.bgz`, `.gvcf.gz`, `.bcf`;
  indexed or not; sorted or not; with or without the `chr` prefix; single-sample or multi-sample.
- **Automatic genome-build detection per file** from `##contig` lengths, header hints and file name,
  so hg19 and hg38 files can be annotated in the same run. Contradictions (e.g. a file named `hg38`
  whose contigs have hg19 lengths) stop the run instead of silently producing wrong annotations.
- **ANNOVAR-friendly normalization** with bcftools: multi-allelic sites are split, spanning deletions
  (`*`), symbolic alleles (`<NON_REF>`, `<DEL>` ...) and empty ALTs are removed; optional
  left-alignment against a reference FASTA, PASS-only and sites-only modes.
- **Fast on many CPUs**: every input file is cut into one chunk per CPU, each annotated by its own
  single-threaded ANNOVAR process and merged back in genomic order — the GIAB HG002 genome (4
  million records, 20 hg38 databases) in **9.4 minutes on 48 CPUs**, 4.3 times faster than
  `table_annovar.pl --thread 48` with identical output, from cohorts of hundreds of samples down
  to single exomes ([benchmarks](#performance)).
- **Databases installed for you**: the first run on an empty humandb downloads the core clinical
  set (gene models, dbSNP, gnomAD, REVEL, ClinVar, InterVar, PubMind) before annotating;
  `--mode check_db` audits every protocol database (installed, intact, newer version, newer ClinVar
  at NCBI) and `--mode update_db` repairs damaged ones and builds the newest **ClinVar straight
  from NCBI** (`clinvar_latest` in a protocol always uses the newest installed release). Every
  install is verified and atomic.
- **Fail-fast pre-flight check** before annotation: missing tables, missing `*Mrna.fa`, stale `.idx`.
- **Clean outputs**: the generic `Otherinfo4..N` columns are replaced by the real VCF column names
  (`#CHROM POS ... FORMAT <samples>`), tables are BGZF-compressed, annotated VCFs are indexed.
- **Clinical report per input file** (HTML, plus TSV and JSON for laboratory systems): carried
  pathogenic / likely pathogenic variants and rare coding VUS with every genotype call's **quality
  flags** (low depth, low GQ, few alt reads, skewed VAF, DRAGEN homology / mosaic / FT filters) and
  a confidence, so questionable calls are marked (P/LP always listed; doubtful candidates and VUS
  behind a switch); HGVS notation (`c.590G>A p.(Arg197His)`); literature evidence from **PubMind**
  (hg38: what papers report about the variant, extracted by a language model); aligned tables whose
  genotype columns follow the input (single sample, one column per family member, cohort
  carriers); database versions; functional classes, SNV/indel, Ti/Tv — computed while streaming,
  with no second pass over the data.
- **Gene context for interpretation**: diseases and modes of inheritance of every gene
  ([GenCC](https://thegencc.org/)) and its loss-of-function constraint (gnomAD LOEUF) next to each
  variant; inheritance candidates checked against the gene's mode of inheritance; **gene panels**
  (gene lists or PanelApp downloads) that keep other genes out of view until asked; one click from a
  gene to ClinVar, gnomAD, Franklin, VarSome, OMIM, GenCC or the locus in IGV.
- **Samples and families**: per-sample QC (sex from chrX vs. the PED file, depth, het/hom, Ti/Tv)
  and what each sample carries (rare coding, P/LP, VUS, inheritance candidates),
  relatedness of the samples of a family (KING kinship, IBS0: parent–child, siblings, sample swaps),
  trio Mendelian errors and an **inheritance analysis** of the affected samples (de novo, homozygous,
  X-linked, compound heterozygous, dominant from an affected parent) with `--ped`.
- **Reproducibility & provenance**: pinned conda environment / Docker image, software versions,
  build-detection reports, Nextflow timeline/report/trace/DAG, `-resume` support.
- **Tested**: Python unit tests, `nextflow lint`, and end-to-end CI runs with a mock ANNOVAR
  (no licence needed) covering every input flavour, with conda and in the Docker image;
  [benchmarks](benchmarks/) against `table_annovar.pl` with a row-by-row comparison.

## Why not just run `table_annovar.pl`?

| | `table_annovar.pl` | annovar-nf |
|---|---|---|
| Input | uncompressed VCF; multi-allelic sites, `*` and symbolic alleles as they are | any VCF/BCF flavour, normalized for ANNOVAR |
| Genome build | given by hand | detected per file; contradictions stop the run |
| Many CPUs | `--thread`: threads wait for each other after every database; genes on one thread below 1 M records | one independent process per CPU, chunks of equal size |
| Output columns | `Otherinfo1..N` | real VCF column names, BGZF + index |
| Databases | downloaded by hand; missing or stale files found mid-run | missing ones installed before the run, pre-flight check, `check_db` / `update_db`, newest ClinVar from NCBI |
| Clinical review | — | per-file report: carried P/LP, VUS, inheritance, genotype quality flags |
| Provenance | — | software and database versions, Nextflow trace, `-resume` |

## Performance

Genome in a Bottle HG002 (public), hg38, 20 databases (gene models ×3, dbSNP, 1000 Genomes ×5,
gnomAD 4.1 exome + genome, ClinVar, InterVar, dbNSFP 5.4a, REVEL, PubMind, ...), 48 CPUs of a
2 × Xeon Gold 6226R server:

| Input | `table_annovar.pl --thread 48` | annovar-nf | Identical rows |
|---|---|---|---|
| whole genome, 4.0 M records | 40 min 45 s | **9 min 23 s** | 4,095,611 of 4,095,611 |
| coding exons (exome-sized), 44 k records | 6 min 0 s | **3 min 10 s** | 44,205 of 44,205 |

One single-threaded ANNOVAR process per CPU keeps every CPU busy (99 %), needs under 1 GB and
uses fewer CPU hours (5.8 against 9.1 for the genome); `--thread` keeps fewer than a third of its
CPUs busy, and one whole-genome run holds 9 GB. Method, per-database costs and the script to
reproduce: [benchmarks/](benchmarks/).

## Pipeline overview

```
input file(s) ─► PREPARE_VCF ──────► NORMALIZE_VCF ─────► ANNOVAR_TABLE ─────► GATHER_CHUNKS ─► annotation/<id>/
                 • detect build      • split multi-       • table_annovar.pl   • header once +
                 • BGZF/sort/index     allelics           • real VCF column      concatenated rows
                 • plan chunks       • drop * / <SYM>       names              • bcftools concat + index
                       │             (per chunk)          (per chunk)                │
                       ├──► DB_INSTALL (databases missing for the build, once) ──►   ▼
                       │    the chunks of that build wait for it        ANNOTATION_SUMMARY ─► summary/
                       └──► SAMPLE_QC (sex, relatedness, PED checks) ──► (QC, inheritance, P/LP, VUS)
```

## Requirements

| Software | Version | Notes |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | ≥ 26.04.6 | Java 17+ |
| [ANNOVAR](https://annovar.openbioinformatics.org/) | 2019Oct24 or later | licensed separately; not bundled |
| bcftools / htslib, Perl, Python ≥ 3.10 | see [`environment.yml`](environment.yml) | provided automatically with `-profile conda`, `docker` or `singularity` |

The Python helpers use only the standard library — no pandas required.

**Platforms.** Any x86-64 Linux with Bash, Java 17+ and one of conda/mamba/micromamba, Docker,
Singularity or Apptainer. The test profile gives the same tables on Ubuntu 20.04, 22.04 and 24.04,
Debian 12, Rocky Linux 8 and 9, Oracle Linux 8 and Ubuntu 22.04 under WSL2 (Windows). On CentOS 7
(glibc 2.17) the conda tools run, but conda's Java does not: use a Java build for glibc 2.17
(e.g. Eclipse Temurin) or a container profile. On WSL2, keep the data and the work directory
inside the Linux file system (`~/...`, not `/mnt/c/...`), and install Java 17+ if the
distribution's is older (e.g. `micromamba create -n nf -c conda-forge -c bioconda openjdk=21 nextflow`).

## Quick start

```bash
git clone https://github.com/assetdaniyarov/annovar-nf.git
cd annovar-nf

# 1. Smoke test with the bundled mock ANNOVAR (no licence, ~1 min)
nextflow run main.nf -profile test,docker        # or test,conda / test,singularity

# 2. Real data: a directory of VCFs, build detected per file
nextflow run main.nf -profile docker \
    --input       /data/vcfs/ \
    --annovar_src /opt/annovar \
    --annovar_db  /opt/annovar/humandb \
    --threads     32 \
    --outdir      results/
```

ANNOVAR itself is not redistributed: [register](https://www.openbioinformatics.org/annovar/annovar_download_form.php)
for a free academic licence, download it, and point `--annovar_src` at it. The container (or conda
environment) provides everything else; Nextflow mounts the ANNOVAR directories.

The databases come by themselves: on the first run the protocol databases missing from
`--annovar_db` (created if needed) are downloaded from ANNOVAR — ClinVar is built from the newest
NCBI release — and verified before annotation starts. Only the builds of your input files are
installed: about 30 GB to download for hg38 (185 GB on disk), 15 GB for hg19 (75 GB). Later runs
find them installed and start at once; `--auto_install_db false` stops instead.

`--input` accepts:

| Form | Example |
|---|---|
| directory | `--input /data/vcfs/` (pattern `--input_pattern`, add `--recursive` for sub-folders) |
| glob (quoted) | `--input '/data/run1/*.hard-filtered.vcf.gz'` |
| single file | `--input HG002.hard-filtered.vcf.gz` |
| samplesheet | `--input samples.csv` (or `.tsv`) |

Samplesheet columns: `sample` (optional; derived from the file name when empty), `vcf` (required;
relative paths are resolved against the samplesheet), `genome` (optional: `hg19`, `hg38`, `auto`).
See [`assets/samplesheet.csv`](assets/samplesheet.csv).

```csv
sample,vcf,genome
HG002,/data/vcfs/HG002.hard-filtered.vcf.gz,
cohort_wes,/data/cohorts/cohort_wes_hg19_dragen.vcf,hg19
panel_lift,/data/cohorts/panel_dragen_hg38_lift.vcf.gz,auto
```

Sample ids derived from file names drop the extension and common pipeline suffixes:
`HG002.hard-filtered.vcf.gz` → `HG002`, `cohort_hg19_v10_genome_dragen.vcf.gz` →
`cohort_hg19_v10_genome_dragen`.

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | samplesheet, directory, glob or file (**required**) |
| `--annovar_src` | — | directory containing `table_annovar.pl` (**required**) |
| `--annovar_db` | — | ANNOVAR `humandb/` directory (**required**; created and filled on the first run) |
| `--outdir` | `./results` | output directory |
| `--genome` | `auto` | `auto`, `hg19` or `hg38`; explicit values still get checked against `##contig` lengths |
| `--skip_build_check` | `false` | downgrade build contradictions to warnings (e.g. stale headers after lift-over) |
| `--protocol` / `--operation` | per build | override the default protocol for every build (operation inferred if omitted) |
| `--threads` | `8` | CPUs for ANNOVAR per input file: one chunk and ANNOVAR process each |
| `--chunk_size` | `1000000` | at most about this many records per chunk; `0` disables splitting |
| `--pass_only` | `false` | keep `FILTER=PASS` records only |
| `--drop_genotypes` | `false` | annotate sites only — strongly recommended for very large cohorts |
| `--min_ac` | `0` | drop alleles with AC below this after splitting (e.g. `1` removes alleles no sample carries) |
| `--fasta_hg19` / `--fasta_hg38` | — | reference FASTA for left-alignment + REF check (`--check_ref`, default `w`) |
| `--compress_output` | `true` | BGZF-compress the annotated tables |
| `--save_vcf` | `false` | also write the annotated VCF (`<id>.<build>_multianno.vcf.gz` + `.tbi`) |
| `--save_normalized` | `false` | publish normalized per-chunk VCFs |
| `--auto_install_db` | `true` | install missing protocol databases before annotating |
| `--skip_db_check` | `false` | skip the humandb pre-flight check (and the install) |
| `--ped` | — | PED file or its CSV/TSV spreadsheet export (family, sample, father, mother, sex, phenotype): inheritance analysis of the affected samples, sex and relationship checks. Optional: otherwise the VCF's `##PEDIGREE` (DRAGEN) or the relatedness of the samples gives the family |
| `--qc_min_dp` / `--qc_min_gq` / `--qc_min_alt_reads` | `10` / `10` / `4` | genotype quality flags `LOW_DP`, `LOW_GQ`, `LOW_ALT_READS` |
| `--qc_het_vaf` / `--qc_min_hom_vaf` | `0.2,0.8` / `0.8` | `VAF_SKEW`: heterozygous VAF outside MIN,MAX; homozygous below |
| `--family_max_samples` | `10` | files with at most this many samples get every genotype listed, relatedness and the inheritance analysis |
| `--gene_panel` | — | gene panel(s), comma-separated: gene lists or PanelApp TSV downloads; the report shows other genes on request |
| `--gene_info` | `true` | diseases, modes of inheritance (GenCC) and LOEUF (gnomAD) per gene; the gene table is built in humandb on the first run |
| `--skip_sample_qc` | `false` | no per-sample QC |

Full list: `nextflow run main.nf --help`, or [`nextflow_schema.json`](nextflow_schema.json) (types,
allowed values; the launch form of Seqera Platform). Details: [docs/usage.md](docs/usage.md).

### ANNOVAR protocols

Default protocols per build live in [`conf/annovar_protocols.config`](conf/annovar_protocols.config)
as `'<database>:<operation>'` entries: the core set for clinical germline review (refGene,
knownGene, dbSNP 151, gnomAD exomes and genomes, REVEL, ClinVar, InterVar, PubMind for hg38).
Override them for your site with a small config file; databases you add are installed on the next
run:

```groovy
// my_site.config — use with: -c my_site.config
params {
    annovar_src = '/opt/annovar'
    annovar_db  = '/opt/annovar/humandb'
    annovar_protocols.hg38 = ['refGene:g', 'gnomad41_genome:f', 'clinvar_latest:f', 'dbnsfp54a:f']
}
```

[`conf/dragen.config`](conf/dragen.config) (`-profile dragen`) is a complete example for a single
shared server.

### Keeping databases up to date

```bash
nextflow run main.nf -profile docker -c my_site.config --mode check_db     # read-only audit
nextflow run main.nf -profile docker -c my_site.config --mode update_db    # repair / newest ClinVar
```

See [docs/databases.md](docs/databases.md) for every check and action.

## Outputs

```
results/
├── annotation/<id>/
│   ├── <id>.<build>_multianno.txt.gz        # main table: annotations + real VCF columns
│   └── <id>.<build>_multianno.vcf.gz(.tbi)  # only with --save_vcf
├── summary/
│   ├── <id>.annotation_summary.html         # report per file: samples, P/LP, inheritance, VUS, versions
│   ├── <id>.annotation_summary.json         # the same for machines (laboratory systems)
│   ├── <id>.annotation_summary_samples.tsv  # per sample: QC, sex check, PED data, P/LP / VUS carried
│   ├── <id>.annotation_summary_plp.tsv       # carried pathogenic / likely pathogenic (per carrier, with call quality)
│   ├── <id>.annotation_summary_vus.tsv       # carried rare coding VUS / conflicting
│   ├── <id>.annotation_summary_candidates.tsv  # inheritance candidates (de novo, recessive, ...)
│   ├── <id>.annotation_summary_totals.tsv
│   ├── <id>.annotation_summary.tsv          # Func × ExonicFunc counts
│   └── annotation_summary.*                 # all files together (only with several inputs)
├── normalized/<id>/                          # only with --save_normalized
└── pipeline_info/
    ├── software_versions.yml                 # tools + every ANNOVAR database used
    └── execution_{report,timeline,trace}_*.html|txt, pipeline_dag_*.html
```

Column layout of `*_multianno.txt.gz`: `Chr Start End Ref Alt <annotation columns> Otherinfo1
Otherinfo2 Otherinfo3 #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample columns>`. Read it with
`zcat`, `pandas.read_csv(..., sep="\t")` or `data.table::fread()`. More: [docs/output.md](docs/output.md).

## Execution profiles

| Profile | Purpose |
|---|---|
| `standard` | local executor, tools from `PATH` (default) |
| `conda` / `mamba` | create the pinned environment from `environment.yml` |
| `docker` / `singularity` / `apptainer` | use the container image (below) |
| `slurm` | submit to SLURM; set `--slurm_queue`, `--cluster_options` |
| `test` | mock ANNOVAR + tiny VCFs covering every input flavour |
| `dragen` | example site profile (paths and database versions of one server) |

Combine them: `-profile slurm,singularity`, `-profile test,docker`.

**Container image.** `ghcr.io/assetdaniyarov/annovar-nf:<version>` (GitHub Container Registry) holds
bcftools/htslib, Perl and Python from [`environment.yml`](environment.yml) — not ANNOVAR, which
Nextflow mounts from `--annovar_src`. The `docker`, `singularity` and `apptainer` profiles pull
the image of the pipeline's version (`--container` overrides it); every release tag builds and
publishes it ([`.github/workflows/docker.yml`](.github/workflows/docker.yml)). To build it
yourself, e.g. without internet access on the compute nodes:

```bash
docker build -t ghcr.io/assetdaniyarov/annovar-nf:0.1.1 .
```

## Large cohorts

For joint-genotyped WGS cohorts (tens of millions of records, hundreds of samples):

- use `--drop_genotypes` unless you need genotypes in the ANNOVAR table — per-sample columns
  dominate file size (genotypes stay available in the original VCF);
- CPUs and RAM are detected automatically; on a shared server cap total usage with
  `--max_cpus` / `--max_memory`. A file gets one chunk per `--threads` CPU (and at least one per
  `--chunk_size` records), so give a large file all the CPUs you can spare.

## Development

```bash
pip install pytest ruff && pytest && ruff check      # Python helpers
nextflow lint .                                       # Nextflow code
nextflow run main.nf -profile test,conda              # end-to-end with the mock ANNOVAR
nextflow run main.nf -profile test -stub              # wiring only
benchmarks/benchmark.sh -h                            # time and compare with table_annovar.pl
```

See [CONTRIBUTING.md](CONTRIBUTING.md) and [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use annovar-nf, please cite it ([CITATION.cff](CITATION.cff); GitHub shows the citation
under "Cite this repository") together with ANNOVAR, Nextflow, BCFtools and the databases of your
protocol and of the gene table (e.g. ClinVar, gnomAD, InterVar, PubMind, GenCC; see
[docs/databases.md](docs/databases.md)):

- Daniyarov A. annovar-nf: scalable ANNOVAR annotation of VCF/BCF files with Nextflow. Zenodo.
  [doi:10.5281/zenodo.23880854](https://doi.org/10.5281/zenodo.23880854) (all versions; the DOI of
  each release is listed on the Zenodo page, e.g. 0.1.0:
  [10.5281/zenodo.23880855](https://doi.org/10.5281/zenodo.23880855)).
- Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from
  high-throughput sequencing data. *Nucleic Acids Res.* 2010;38(16):e164.
- Di Tommaso P, et al. Nextflow enables reproducible computational workflows.
  *Nat Biotechnol.* 2017;35:316–319.
- Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008.

## License

annovar-nf is MIT-licensed © Asset Daniyarov. It is an independent project, not affiliated with or
endorsed by the authors of ANNOVAR, and it contains no ANNOVAR code or data: it runs the copy of
ANNOVAR you point it to (the tests use a mock).

- **ANNOVAR** is free for personal, academic and non-profit use and may not be redistributed;
  commercial users and service providers need a licence from QIAGEN
  ([download page](https://www.openbioinformatics.org/annovar/annovar_download_form.php)).
  Every user of annovar-nf needs their own ANNOVAR download.
- **Annotation databases** have their own terms (e.g. COSMIC, dbNSFP components; PubMind is for
  academic, non-commercial use); check them for your use, especially clinical or commercial use.
