# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- `--mode update_db` installs the 1000 Genomes packages (`.zip`, e.g. `1000g2014oct`) in the
  container image: the archives are downloaded and unpacked by `annovar_db.py`, which needs neither
  `wget` nor `unzip`.

## [0.1.1] - 2026-11-02

### Added
- dbSNP rsID of every listed variant: in the report under the change, linked to dbSNP, in the
  gene link menu, and as `rsid` in the TSV and JSON reports (from the protocol's `avsnp` database,
  else the VCF `ID` column).

### Fixed
- Mitochondrial variants are annotated by every filter database whatever name it uses for the
  mitochondrion (`M` or `MT`): `chrM` variants get dbSNP, PubMind, CADD and the ClinVar files
  distributed by ANNOVAR, which list it as `MT` (hg38, and hg19 files with the rCRS sequence).

## [0.1.0] - 2026-11-01

### Added
- One code base for hg19 and hg38; requires Nextflow ≥ 26.04.6 (strict syntax).
- Input as samplesheet (CSV/TSV, optional `genome` column), directory (`--input_pattern`,
  `--recursive`), glob or single file; `.vcf`, `.vcf.gz` (BGZF or gzip), `.vcf.bgz`, `.gvcf.gz`, `.bcf`.
- `PREPARE_VCF`: per-file genome-build detection (`##contig` lengths, header hints, file name) with
  conflict checks; BGZF compression, sorting and indexing only when needed.
- `NORMALIZE_VCF`: multi-allelic splitting; removal of `*`, symbolic and empty ALT alleles; optional
  left-alignment (`--fasta_hg19/--fasta_hg38`), `--pass_only`, `--drop_genotypes`, `--min_ac`.
  Chunk regions are passed to bcftools as a regions file, so contig names with `:` and `*`
  (HLA alleles) work.
- Scatter/gather with one chunk per CPU (`--threads`), each annotated by its own single-threaded
  ANNOVAR process at the same time: the primary chromosomes are cut into pieces of equal records,
  other contigs packed separately; at most about `--chunk_size` and at least 1,000 records per
  chunk. ANNOVAR task memory follows the chunk size. `benchmarks/benchmark.sh` times a run against a
  single `table_annovar.pl` and compares the tables row by row.
- Pre-flight check of every ANNOVAR database in the protocol (`--skip_db_check`) and of the ANNOVAR
  scripts (executable bit).
- Database maintenance: `--mode check_db` / `--mode update_db` and `bin/annovar_db.py` (audit,
  verify, download, clinvar, install); protocol databases that are not installed are installed
  before annotation (`--auto_install_db`, on by default), so the first run on an empty humandb
  downloads them; the database modes work in containers (humandb mounted). Verified atomic
  installs; resumable parallel downloads (32 MB ranges over 8 connections) with size checks;
  ClinVar built from the NCBI VCF; `<family>_latest` protocol aliases (e.g. `clinvar_latest`);
  newer versions suggested only if published on the ANNOVAR download page.
- Protocols per build in `conf/annovar_protocols.config`: a core set for clinical germline review
  (refGene, knownGene, dbSNP 151, gnomAD exomes + genomes, REVEL, ClinVar, InterVar, PubMind for
  hg38); an extended set in `conf/dragen.config`; `--operation` inferred from `--protocol`.
- Outputs: `annotation/<id>/<id>.<build>_multianno.txt.gz` (real VCF column names, BGZF);
  annotated VCF + index with `--save_vcf`. Genes in `Gene.*` / `GeneDetail.*` are separated by `;`
  with every ANNOVAR release.
- Report per input file `summary/<id>.annotation_summary.html` (+ TSV files; combined report for
  several input files): carried pathogenic / likely pathogenic variants (ClinVar, plus InterVar-only
  calls) and rare coding ClinVar VUS / conflicting, with zygosity, VAF, depth, FILTER, review stars,
  condition, ACMG criteria, PubMind literature evidence (hg38), gnomAD frequency, REVEL, CADD,
  AlphaMissense and ClinVar / gnomAD links;
  HGVS notation of the first RefSeq transcript (`c.590G>A p.(Arg197His)`, splice sites from
  `GeneDetail`) and Sequence Ontology consequences; variant tables with one column set of fixed
  widths, genotypes laid out for the input (single sample, family member columns with one row per
  variant, cohort carriers) and doubtful candidates / VUS behind a switch;
  reference calls (`0/0`) at pathogenic sites are counted, not listed; genome-build evidence;
  genomic regions and exonic consequences; annotation sources with database versions (key versions
  also in the report header).
- Genotype call quality on every listed variant: zygosity, VAF, depth, GQ, alt reads, flags
  (`LOW_DP`, `LOW_GQ`, `LOW_ALT_READS`, `VAF_SKEW`, DRAGEN `MOSAIC` / `HOMOLOGY` / `FT`,
  `ALT_READS_IN_REF`) and a confidence (high / caution / low); thresholds `--qc_*`. Mitochondrial
  calls are reported as heteroplasmy levels (no `VAF_SKEW`).
- Multi-sample files: every sample's call for families (`--family_max_samples`), carriers for
  cohorts; one TSV row per variant and carrier sample.
- `SAMPLE_QC`: per-sample counts, sex from chrX, pairwise KING kinship / IBS0 (parent–child,
  siblings, unrelated), trio Mendelian errors; checked against `--ped`.
- Gene information in the report (`--gene_info`): diseases and modes of inheritance from GenCC and
  gnomAD v4.1 LOEUF / pLI / missense Z next to every gene, built into humandb as
  `annovar_nf/genes.tsv` by `annovar_db.py genes` (first run, `--mode update_db` after 30 days);
  inheritance candidates marked and ordered by whether the gene's mode of inheritance fits the
  model; HGNC previous symbols resolved.
- Gene panels (`--gene_panel`): gene lists or PanelApp TSV downloads (green genes); variants in
  other genes hidden in the report until asked, panel counts in the section headings, `panels`
  column in the TSV files.
- Link menu on every gene name: the variant in ClinVar, gnomAD, Franklin and VarSome, the gene in
  OMIM and GenCC, the locus in a running IGV.
- What every sample of a multi-sample file carries — rare coding variants, ClinVar P/LP,
  InterVar-only P/LP, rare coding VUS, inheritance candidates per model — in the report
  (*Findings per sample*), the samples TSV and the JSON report.
- Inheritance analysis of the affected samples (`--ped`; a single sample without it): de novo,
  homozygous, X-linked, compound heterozygous (phased by the parents), dominant from an affected
  parent, with segregation against unaffected relatives. Without `--ped` the pedigree comes from
  the VCF's `##PEDIGREE` (DRAGEN) or is inferred from the relatedness, and the children are the
  probands; DRAGEN de novo calls (`DN` / `DQ`) are shown. Without sequenced parents, homozygous
  candidates are damaging variants only; low-confidence calls and Mendelian inconsistencies are
  marked doubtful.
- `<id>.annotation_summary.json`: machine-readable report (versions, QC, samples, variants with
  genotype calls and inheritance) for a laboratory information system.
- `pipeline_info/software_versions.yml` with tool versions and every ANNOVAR database used per build
  (resolved name, file, file date, size).
- Resources: `--max_cpus` / `--max_memory` default to the machine's CPUs / RAM, larger requests are
  capped; `--threads` per input file. Deterministic errors are not retried.
- Boolean parameters take `true` / `false` on the command line (`--onetranscript false`); other
  values stop the run with a message.
- ANNOVAR is stopped together with its task (Ctrl+C, kill).
- `-resume` reruns a step after its `bin/` scripts change (their checksum is part of the task), so
  an updated pipeline never reuses stale tables or reports.
- Final console summary sized for Nextflow's live display: the database check lists every database
  in columns when it fits the window (else one line per build), at most 10 input files are listed;
  full tables with `-ansi-log false` and in the report files.
- Profiles: `conda`, `mamba`, `docker`, `singularity`, `apptainer`, `slurm`, `test`, `dragen`.
- `nextflow_schema.json` describing every parameter (checked against `nextflow.config` by the tests).
- `environment.yml` (GSL pinned to the version bcftools 1.21 links), `Dockerfile`, mock ANNOVAR and
  test VCFs, Python unit tests, GitHub Actions CI (pytest, ruff, `nextflow lint`, end-to-end runs
  with conda and in the Docker image, stub run), container image `ghcr.io/assetdaniyarov/annovar-nf`
  built and published for every release tag. The test profile gives the same tables on Ubuntu
  20.04–24.04, Debian 12, Rocky Linux 8/9, Oracle Linux 8 and WSL2.
- Benchmarks on public data (GIAB HG002 genome and its coding exons) against `table_annovar.pl`.
- Report tables show long text (conditions, HGVS) shortened and in full on hover; tooltips appear
  only when they add something. Chunk progress lines in plain logs only (the live display counts
  the ANNOVAR tasks on one line).
