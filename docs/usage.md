# annovar-nf: usage

## Inputs

### What files are accepted

| Input | Handling |
|---|---|
| `.vcf` (plain text) | compressed with BGZF, indexed |
| `.vcf.gz` BGZF + `.tbi`/`.csi` next to it | used in place (symlinked, never rewritten) |
| `.vcf.gz` BGZF without index | indexed |
| `.vcf.gz` compressed with plain `gzip` | re-compressed with BGZF, indexed |
| `.bcf` | converted to BGZF VCF, indexed |
| unsorted records (typical after lift-over) | sorted with `bcftools sort` |
| `chr1` or `1` contig names | both supported (ANNOVAR strips the `chr` prefix) |
| gVCF (`.g.vcf.gz`, `.gvcf.gz`) | reference blocks and `<NON_REF>` alleles are removed |
| multi-sample / joint-genotyped VCF | supported; genotypes kept unless `--drop_genotypes` |

Index files are found automatically when they sit next to the VCF as `<file>.tbi` or `<file>.csi`.

### Directory and glob inputs

```bash
--input /data/vcfs/                                   # *.{vcf,vcf.gz,vcf.bgz,gvcf.gz,bcf}
--input /data/vcfs/ --recursive                       # include sub-directories
--input /data/vcfs/ --input_pattern '*.hard-filtered.vcf.gz'
--input '/data/run*/vcf/*.vcf.gz'                     # quote globs!
```

Duplicate sample ids stop the run; use a samplesheet to name files explicitly.

### Samplesheet

CSV (or TSV with a `.tsv` extension) with a header line:

| Column | Required | Description |
|---|---|---|
| `vcf` | yes | path to the variant file; relative paths are resolved against the samplesheet |
| `sample` | no | output id; derived from the file name when empty |
| `genome` | no | `hg19`, `hg38` (aliases: `GRCh37`, `b37`, `GRCh38`) or `auto` |

## Genome build

With `--genome auto` (default) every file is classified independently:

1. **`##contig` lengths** of chromosomes 1–22, X, Y — authoritative.
2. **Header hints** — `##reference=...hg38.fa`, DRAGEN/GATK command lines, `assembly=` tags.
3. **File-name tokens** — `hg19`, `GRCh37`, `b37`, `hs37d5`, `hg38`, `GRCh38`, `hs38DH`
   (as separate words: `cohort_hg19_dragen.vcf` matches, `sample_hg190.vcf` does not).

If sources disagree, or nothing identifies the build, the run stops with an explanation. A build
set explicitly (`--genome hg38` or the samplesheet `genome` column) is still compared with the
contig lengths. `--skip_build_check` turns contradictions into warnings — use it only when you know
the header is stale (e.g. a lift-over tool that did not rewrite `##contig` lines).

The build of every file, how it was determined and the evidence (`##contig` length votes, header
hints, file-name tokens, warnings) are shown in the "Genome build" section of
`summary/<id>.annotation_summary.html`.

## Normalization

Per chunk, before ANNOVAR:

```text
bcftools view   [--regions <chunk> --regions-overlap pos] [--apply-filters PASS,.] [--drop-genotypes]
bcftools norm   --multiallelics -any [--fasta-ref <fasta> --check-ref w]
bcftools view   --exclude 'ALT="*" || ALT~"^<" || ALT="."' [--min-ac N]
```

- Without `--fasta_hg19/--fasta_hg38` indels are not left-aligned (DRAGEN and GATK output is
  already left-aligned). With a FASTA, REF mismatches are reported (`--check_ref w`), excluded
  (`x`) or fixed (`s`).
- `--min_ac 1` removes alleles that no sample carries after splitting multi-allelic sites —
  useful for cohort subsets.

## Parallelism

`PREPARE_VCF` reads per-contig record counts from the index and plans one chunk per CPU
(`--threads`), and at least `total / --chunk_size` chunks, all of about the same number of
records:

- the primary chromosomes are cut, in header order, into pieces of equal records (records are
  assumed to spread evenly along a chromosome); a cut close to a chromosome end moves to that end.
  Names show the first and last piece: `chr1_3-chr2_1` (end of chr1, start of chr2),
  `chr13-chr16` (whole chromosomes);
- the remaining contigs (alt, decoy, HLA, EBV ...) are never cut and get chunks of their own
  (`other_1` ...);
- no chunk is made smaller than 1,000 records to use more CPUs, and a file below that is one
  chunk (`all`); `--chunk_size 0` never splits.

Each chunk is one `NORMALIZE_VCF` + one `ANNOVAR_TABLE` task, and all chunks of a file run at
the same time, each as one single-threaded ANNOVAR process. This is the fastest way to run
ANNOVAR on many CPUs:

- `table_annovar.pl --thread N` runs Perl threads that wait for each other after every
  database, and annotates genes with a single thread below a million records, so on the GIAB
  HG002 genome its 48 threads kept fewer than a third of their CPUs busy on average; independent
  processes stay at 100 %;
- each process looks up its own genomic range through the database indexes, while a process
  that needs more than 90 % of a database's index bins (a whole genome at once) reads the
  whole database file instead;
- every process loads the database indexes once (about a minute with a large hg38 protocol),
  which the 1,000-record minimum keeps reasonable next to the annotation itself.

A chunk gets more than one CPU (and ANNOVAR `--thread`) only when a file has fewer chunks than
CPUs. Results are concatenated in genomic order. [benchmarks/](../benchmarks/) measures the
layouts against each other and against a single `table_annovar.pl` run.

## ANNOVAR protocol

Defaults: [`conf/annovar_protocols.config`](../conf/annovar_protocols.config). Three ways to change it:

```bash
# 1. one run, all builds
--protocol refGene,gnomad41_genome,clinvar_20240917 --operation g,f,f

# 2. operation inferred: gene models → g, everything else → f
--protocol refGene,ensGene,gnomad41_genome,clinvar_20240917

# 3. per build, in a config file (-c site.config)
params.annovar_protocols.hg19 = ['refGene:g', 'gnomad211_genome:f', 'clinvar_20240917:f']
```

Further `table_annovar.pl` options go into `--annovar_args`, e.g.
`--annovar_args '--intronhgvs 50'`. `--onetranscript` is on by default (`--onetranscript false`
reports all transcripts).

Every entry needs `<annovar_db>/<build>_<database>.txt`. Entries that are not installed are
installed before the files of that build are annotated (ANNOVAR packages; ClinVar from NCBI), so
a database added to the protocol arrives with the next run; `--auto_install_db false` stops with
the list instead. A database ANNOVAR does not distribute (your own tables) has to be placed in
humandb yourself. See [databases.md](databases.md) for the sizes of the default set.

## Families and genotype quality

A multi-sample VCF (e.g. DRAGEN joint calling of a family) is one input file; its report shows
every sample's genotype call for each listed variant, one column per family member. Give the
family structure in a PED file:

```
# family  sample  father  mother  sex  phenotype   (sex: 1 male, 2 female, 0 unknown;
FAM1      F1_01   0       0       1    1            #  phenotype: 2 affected, 1 unaffected,
FAM1      F1_02   F1_01   0       2    2            #  0 unknown; parents not sequenced: 0)
FAM1      F1_03   F1_01   0       2    1
```

```bash
nextflow run main.nf -profile docker -c my_site.config --input FAM1.vcf.gz --ped families.ped
```

One PED file can describe the families of every input file. A CSV or TSV export of a spreadsheet
works too, as long as its first six columns are the PED columns (a header row is skipped). Samples
are matched to the VCF by name, ignoring case and the separators `.`, `_`, `-` (`F1.01` in the
PED matches `F1_01` in the VCF). The PED data drive the inheritance analysis (affected samples
are the probands) and the checks of sex and relationships against the data.

`--ped` is optional. For the samples it does not cover, the pedigree is taken from the VCF header
(`##PEDIGREE=<Child=…,Mother=…,Father=…>`, written by DRAGEN joint calling with a pedigree), else
inferred from the relatedness of the samples (a sample in a parent–child relationship with two
mutually unrelated samples is their child; a sample in a parent–child relationship with two full
siblings is their parent). Neither source knows who is affected, so the children are analysed as
probands, and nobody is excluded as unaffected. A single-sample file is its own proband
(homozygous, hemizygous, possible compound heterozygous; without parents homozygous and unphased
compound heterozygous candidates are damaging variants only). DRAGEN's own de novo calls
(`FORMAT/DN`, `DQ`) are shown with the genotype and in the de novo candidates. Candidates resting
on a low-confidence call or on genotypes that break Mendelian inheritance (a homozygous child of a
homozygous-reference parent) are marked doubtful and hidden in the HTML report until asked.

Genotype quality flags (`--qc_min_dp`, `--qc_min_gq`, `--qc_min_alt_reads`, `--qc_het_vaf`,
`--qc_min_hom_vaf`) mark calls to confirm before reporting; see [output.md](output.md#genotype-quality).
Files with more than `--family_max_samples` samples (cohorts) list carriers only and skip the
relatedness and inheritance analysis.

## Gene panels and gene information

```bash
nextflow run main.nf -profile docker -c my_site.config --input FAM1.vcf.gz --ped families.ped \
    --gene_panel epilepsy_panelapp.tsv
```

`--gene_panel` takes one or more files, comma-separated: a list of gene symbols (one per line;
other columns, `#` comments and a header line are ignored) or a PanelApp TSV download (*Download
panel (tsv)*; the green genes are used, and the panel is named after it). The report then hides
variants in other genes until asked, and counts the variants in the panel genes in every section
heading; the TSV and JSON files keep every variant (see
[output.md](output.md#gene-information-links-and-gene-panels)).

The diseases, modes of inheritance and LOEUF shown next to every gene come from the gene table in
humandb (`--gene_info`, default `true`). The first run builds it (about 140 MB downloaded from
GenCC, gnomAD and HGNC, 2 MB on disk); `--mode update_db` rebuilds it when it is older than 30
days. `--gene_info false` leaves the gene information out.

## Resources

Labels and defaults are in [`conf/base.config`](../conf/base.config). Jobs killed for memory or
time (exit codes 130–145, 104, 175) are retried with more resources (×2, then ×3); other errors
stop the run. Global ceilings: `--max_cpus`, `--max_memory`, `--max_time`.

An `ANNOVAR_TABLE` task asks for 2 GB + 1 GB per CPU + 20 × the size of its compressed chunk: a
single-threaded ANNOVAR process needs under 1 GB, plus its query lines, which carry every sample
column of the input. Small requests let the local executor run one chunk per CPU at the same
time; a whole-genome chunk of a single sample asks for about 3 GB.

`--max_cpus` and `--max_memory` default to the CPUs and RAM of the machine Nextflow runs on
(`-profile dragen`: 3/4 of them, leaving headroom for other users; `-profile slurm`: 128 CPUs /
512 GB). Any request above the ceiling is capped, so the same command works on a 16-core
workstation and a 128-core server: `--threads 120` simply becomes 64 on a 64-core machine, and
the run banner shows the effective value. With the local executor the ceilings are also the total
capacity tasks are scheduled into.

Override per process in a config file:

```groovy
process {
    withName: 'ANNOVAR_TABLE' { memory = 64.GB; time = 96.h }
}
```

## Running on a cluster

```bash
nextflow run main.nf -profile slurm,singularity \
    --slurm_queue normal --cluster_options '--account=genomics' \
    --input samples.csv --annovar_src /shared/annovar --annovar_db /shared/annovar/humandb
```

Run Nextflow itself inside `screen`/`tmux` or with `-bg`, and add `-resume` after any
interruption — finished chunks are not recomputed.
