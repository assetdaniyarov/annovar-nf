# annovar-nf: output

## `annotation/<id>/`

### `<id>.<build>_multianno.txt.gz`

Tab-separated ANNOVAR table, BGZF-compressed (`--compress_output false` writes plain text).
One row per variant allele after splitting multi-allelic sites.

| Columns | Content |
|---|---|
| `Chr Start End Ref Alt` | ANNOVAR coordinates (indels: `-` for the empty allele, `Start` shifted for deletions) |
| annotation columns | one or more per database in the protocol (e.g. `Func.refGene`, `gnomad41_genome_AF`, `CLNSIG`) |
| `Otherinfo1` | allele frequency/zygosity computed by ANNOVAR from the genotypes |
| `Otherinfo2` | QUAL |
| `Otherinfo3` | read depth |
| `#CHROM ... FORMAT <samples>` | the original VCF record, including every sample column |

Missing values are `.`.
Several genes in `Gene.*` / `GeneDetail.*` (overlapping genes, both neighbours of an intergenic
variant) are separated by `;` with every ANNOVAR release, e.g. `SHOX;CRLF2` and `dist=1234;dist=5678`;
in the annotated VCF (`--save_vcf`) the `;` is written as `\x3b`. Transcripts in `AAChange.*` stay
separated by `,`.

```bash
zcat HG002.hg38_multianno.txt.gz | head -3 | cut -f1-12
```

```python
import pandas as pd
df = pd.read_csv("HG002.hg38_multianno.txt.gz", sep="\t", low_memory=False)
```

```r
df <- data.table::fread("HG002.hg38_multianno.txt.gz")
```

**How ANNOVAR matches filter databases.** `annotate_variation.pl -filter` looks a variant up by
chromosome, position and alternate allele, not the reference allele, and reports the first
database record it reads. A database that lists one position and alternate allele with several
reference alleles can therefore report another allele's value, and which one depends on the
order the records are read — how ANNOVAR batches the input. This is ANNOVAR's behaviour, in
`table_annovar.pl` as in annovar-nf. Example: hg38 `regsnpintron` lists `T>C` and `A>C` at
chr1:148,541,475 (NBPF14), so with two layouts of the same whole genome a few `regsnp_fpr` values
in the 1q21 segmental duplications can differ.

**Mitochondrial variants.** ANNOVAR matches the chromosome name exactly (without `chr`), and its
databases name the mitochondrion either `M` (gene models, dbNSFP, gnomAD) or `MT` (avsnp, PubMind,
CADD, the ClinVar files ANNOVAR distributes), so `table_annovar.pl` leaves a `chrM` variant without
the `MT`-named annotations, and an `MT` variant without the `M`-named ones. annovar-nf annotates
every mitochondrial record under both names and fills the empty filter-database columns of the row
from the second lookup (gene-based columns come from the record's own name; a value found under
both names is kept as is). This applies to hg38 and to hg19 files whose mitochondrion is the rCRS
sequence (16,569 bp, e.g. b37 `MT`); UCSC hg19 `chrM` is another sequence and is annotated under its
own name only. The ClinVar database that `--mode update_db` builds lists mitochondrial records
under both names.

### `<id>.<build>_multianno.vcf.gz` (+ `.tbi`, only with `--save_vcf`)

The normalized input VCF with ANNOVAR annotations added to `INFO` (between `ANNOVAR_DATE` and
`ALLELE_END`), coordinate-sorted and tabix-indexed — convenient for `bcftools`, IGV or VEP-style tools.

`table_annovar.pl` logs are not published; when ANNOVAR fails, the end of its log is printed in the
error message (the full log stays in the task work directory).

## `summary/`

One set of files per input file, named after it:

| File | Content |
|---|---|
| `<id>.annotation_summary.html` | self-contained report; the header names the key database versions (ClinVar release, gnomAD, dbNSFP, InterVar, dbSNP, ANNOVAR). Sections in this order: totals; input files; samples (QC, sex check, what each sample carries, relatedness, trios); pathogenic / likely pathogenic variants; inheritance candidates; variants of uncertain significance; genome build (how it was determined and the evidence); genomic regions; exonic consequences; annotation sources (every ANNOVAR database used: source, version, type, file date, size); all metrics |
| `<id>.annotation_summary.json` | the same content for machines, e.g. a laboratory information system (see below) |
| `<id>.annotation_summary_samples.tsv` | one row per sample: PED data, sex declared / inferred from chrX and the check, mean depth, het / hom SNVs, Ti/Tv, indels, singletons, missing genotypes; what the sample carries: rare coding variants, ClinVar P/LP, InterVar-only P/LP, rare coding VUS, inheritance candidates (`NA` when the database is not in the protocol) |
| `<id>.annotation_summary_plp.tsv` | carried pathogenic / likely pathogenic variants (see below) |
| `<id>.annotation_summary_vus.tsv` | carried rare coding/splicing ClinVar VUS and conflicting classifications |
| `<id>.annotation_summary_candidates.tsv` | inheritance candidates (see below) |
| `<id>.annotation_summary_totals.tsv` | records, SNV/indel/MNV, Ti/Tv, PASS, region classes, rare coding, P/LP and VUS counts, REVEL ≥ 0.75, CADD phred ≥ 20 (`NA` when the database is not in the protocol) |
| `<id>.annotation_summary.tsv` | `Func.refGene × ExonicFunc.refGene` counts |

With more than one input file the same files without the `<id>.` prefix (and without the JSON)
cover all files together (one column / row / card per file).

Counts are variant records per input file (not per genotype); for a multi-sample VCF one record
is counted once regardless of how many samples carry it. The same counts per sample — the records
whose genotype has the alternative allele — are in the samples TSV and, for more than one sample,
in the *Findings per sample* table of the report. In the variant TSVs `file` is the input file and
`sample` the VCF sample: there is one row per variant and carrier sample.

### Report layout

The sections keep the same order for every input; the variant tables (P/LP, inheritance
candidates, VUS) follow its type:

| Input | Genotype column(s) | Rows |
|---|---|---|
| single sample | one: zygosity, VAF, depth, flags | one per variant |
| family (2 – `--family_max_samples` samples) | one per member; the header gives sex (♂ / ♀) and role | one per variant, the models of every proband in it |
| cohort, or the combined report of several files | the carriers (at most 12 shown) | one per variant |
| sites-only VCF | none | one per variant |

All variant tables have the same columns with the same fixed widths, so a column lines up from
one table to the next; long values end in `…` and the tooltip holds the full text (for genotypes
also GQ, alt reads and what each flag means). A single-sample or trio report fits a 1366-pixel
screen; wider families scroll inside the table. Doubtful rows — inheritance candidates and VUS
whose calls have low confidence, Mendelian inconsistencies — are hidden behind a *show doubtful*
switch above the table; P/LP variants are always shown. For a single sample the samples table
leaves out role, singletons and missing genotypes, and there is no *Findings per sample* table
(the file card gives the same counts).

Variant changes are written in HGVS for the first RefSeq transcript of ANNOVAR's refGene:
`c.590G>A p.(Arg197His)` (from `c.G590A:p.R197H`), splice sites from `GeneDetail.refGene`
(`c.194+2T>C`); consequences use Sequence Ontology words (missense, stop gained, frameshift,
inframe deletion, splice site ...). The variant TSVs and the JSON report carry `transcript`,
`hgvs_c` and `hgvs_p` next to ANNOVAR's `aachange`. Transcript versions are not in refGene.

### Genotype quality

Every listed genotype call carries its zygosity (`ref` / `het` / `hom` / `hemi` / `missing`),
variant allele fraction (VAF, from `AD`, else `FORMAT/AF`), depth, GQ, alt-allele reads, quality
flags and a confidence. The thresholds are parameters.

| Flag | Meaning | Default |
|---|---|---|
| `LOW_DP` | read depth below `--qc_min_dp` | 10 |
| `LOW_GQ` | genotype quality below `--qc_min_gq` | 10 |
| `LOW_ALT_READS` | carrier with fewer than `--qc_min_alt_reads` reads for the variant allele | 4 |
| `VAF_SKEW` | heterozygous VAF outside `--qc_het_vaf`, homozygous / hemizygous below `--qc_min_hom_vaf` | 0.2–0.8, 0.8 |
| `MOSAIC` | DRAGEN `FORMAT/MSC=1`: possible mosaic variant | |
| `HOMOLOGY` | DRAGEN `FORMAT/HML=1`: region of homology (e.g. a paralogous gene) | |
| `ALT_READS_IN_REF` | reference call with ≥ 2 reads for the variant allele (relevant for parents of a de novo) | |
| `FT:<filter>` | the caller filtered this genotype (`FORMAT/FT`, DRAGEN joint calling) | |

Confidence: **low** when the caller filtered the site (`FILTER`) or the genotype (`FT`), or with
two or more flags; **caution** with one flag; **high** otherwise. Low-confidence calls are listed
in the TSV and JSON files and need confirmation, e.g. Sanger or IGV, before they are reported;
the HTML report hides doubtful candidates and VUS until asked (see Report layout).

On the mitochondrial genome (`chrM` / `MT`) the VAF is the heteroplasmy level: such calls get no
`VAF_SKEW` and the report shows them as *heteroplasmy* / *homoplasmy* with the level in percent.

### Samples, relatedness and pedigree

`SAMPLE_QC` runs `bcftools stats` per input file (PASS records): per-sample counts and the
chrX ratio of heterozygous to homozygous + hemizygous SNVs outside the pseudo-autosomal regions
(male ≤ 0.25, female ≥ 0.6). For files with at most `--family_max_samples` samples (default 10)
every pair gets the KING-robust kinship and the opposite-homozygote (IBS0) rate on autosomal SNVs
with GQ ≥ 20 and depth ≥ 10: about 0.5 duplicate, 0.25 first-degree (IBS0 near 0: parent–child,
otherwise full siblings), 0.125 second-degree, 0 unrelated. The pedigree comes from `--ped`, else
from the VCF header (`##PEDIGREE`, DRAGEN), else it is inferred from these relationships (see
[usage.md](usage.md#families-and-genotype-quality)); the report says which. Declared sex and
relationships are checked against the data (✘ mismatch: sample swap or wrong pedigree) and every
complete trio gets its Mendelian error rate.

### Inheritance candidates

For single samples and families (at most `--family_max_samples` samples), the rare protein-altering
or splicing variants of every affected sample (`--ped`, phenotype 2; with a pedigree from the VCF
header or inferred: the children; without any: the sample of a single-sample file) are matched
against the inheritance models below. A candidate must segregate: an unaffected relative with the
same genotype excludes it.

| Model | Rule | gnomAD max AF |
|---|---|---|
| de novo | heterozygous (hemizygous on chrX in males), father and mother called reference without reads for the variant; DRAGEN's de novo call (`DN=DeNovo`) when the VCF has one, else no low-confidence call in the trio | < 0.1 % |
| homozygous (AR) | homozygous; no unaffected relative and no parent (unless affected) homozygous; with neither parent sequenced, damaging variants only (loss of function, splicing, REVEL ≥ 0.5 or CADD ≥ 20) | < 1 % |
| X-linked (hemizygous) | male, chrX outside the PARs; no unaffected male carrier | < 1 % |
| compound heterozygous | ≥ 2 heterozygous variants in one gene, one from each parent; variants homozygous in a parent (unless affected) left out | < 1 % |
| possible compound heterozygous | the same with one parent (one variant from it, one not) or none sequenced (phase unknown; damaging variants only: loss of function, splicing, REVEL ≥ 0.5 or CADD ≥ 20), or one variant de novo | < 1 % |
| dominant (affected parent) | heterozygous, carried by an affected parent | < 0.1 % |

Splicing means a splice site of a protein-coding transcript (`ncRNA_splicing` is left out).
Sorted by model, then whether the gene's mode of inheritance fits the model (see Gene information),
then the confidence of the proband's call (high first), then loss of function, then
REVEL ≥ 0.5 / CADD ≥ 20. The P/LP and VUS tables show the model of a variant in their first
column.

A candidate is **doubtful** — kept in the TSV (`doubtful` column) and JSON, hidden in the HTML
report until asked — when the proband's call has low confidence or the genotypes break Mendelian
inheritance: a homozygous child of a parent called homozygous reference, usually a genotyping error
and rarely a deletion on the other allele or uniparental disomy.

### JSON report

`<id>.annotation_summary.json` (schema `annovar-nf.report`, version 1) holds: pipeline and ANNOVAR
versions with every database (source, version, file date); the input file (build, how it was
determined, samples); the genotype-quality thresholds; totals; per-sample QC with the PED data;
what every sample carries (`sample_findings`, with the inheritance candidates per model);
relatedness and trios; the sources of the gene table (`gene_information`) and the gene panels; and the `pathogenic`, `uncertain_significance` and
`inheritance_candidates` variants with their annotation (including `transcript`, `hgvs_c`,
`hgvs_p`), every genotype call (quality flags, confidence) and inheritance models per sample;
a candidate's `doubtful` gives the reason it is doubtful, else `null`.

### Pathogenic / likely pathogenic and VUS

Only variants the sample carries are listed: the genotype has an alternative allele. Reference
calls (`0/0`, `./.`) at known pathogenic sites — DRAGEN's targeted callers, for example, write
such records with FILTER `TargetedWT` / `TargetedLowVQL` — are counted
(`clinvar_plp_not_carried`) but not listed. For files without genotypes every record counts.

| List | Selection |
|---|---|
| P/LP | ClinVar `Pathogenic`, `Likely_pathogenic`, `Pathogenic/Likely_pathogenic` (not conflicting), source `ClinVar`; plus InterVar `Pathogenic` / `Likely pathogenic` calls that ClinVar does not classify as P/LP, source `InterVar` |
| VUS | ClinVar `Uncertain_significance` and `Conflicting_classifications_of_pathogenicity`, exonic or splicing, gnomAD frequency below 1 % in every population (or absent) |

Columns of the HTML tables: source (P/LP) and inheritance model; gene (see Gene information
below); HGVS change, consequence, transcript, VCF coordinates and the dbSNP rsID linked to dbSNP (a
non-PASS FILTER highlighted);
genotype calls (see Report layout and Genotype quality); ClinVar significance with review stars
(0–4); condition (ClinVar), with the GenCC diseases of the gene below it; InterVar call with the triggered ACMG criteria; PubMind (hg38, when
`pubmind_db_*` is in the protocol): the pathogenicity reported in the literature, the paper-level
score (0 benign, 1 pathogenic) and the evidence confidence (0–3), with the PubMind variant ids in
the tooltip; highest gnomAD frequency (overall and population maximum, exome and genome); REVEL;
CADD phred; AlphaMissense (with dbNSFP). The TSV files also hold `rsid`, FILTER, dbscSNV ADA (hg19)
and `pubmind`, `pubmind_score`, `pubmind_confidence`, `pubmind_pvid`. The scores are taken from the
protocol databases and are `.` when a database is not used. The rsID comes from the dbSNP database
of the protocol (`avsnp151`, matched per alternate allele), else from the VCF `ID` column when it
holds rsIDs.

PubMind is extracted from papers by a language model: it points to the literature on a variant and
does not classify it. The P/LP list stays ClinVar and InterVar only. In the GIAB HG002 genome, 16
rare coding variants are pathogenic or likely pathogenic in PubMind but not P/LP in ClinVar or
InterVar (11 of them absent from ClinVar, 4 VUS or conflicting, 1 benign); check the papers
(search the PVID on [pubmind.wglab.org](https://pubmind.wglab.org/)) before using such a call.

### Gene information, links and gene panels

With `--gene_info` (default) the report uses the gene table that `annovar_db.py genes` builds in
humandb (see [databases.md](databases.md#gene-table)): for every gene the diseases from
[GenCC](https://thegencc.org/) with their mode of inheritance and the strength of the evidence, and
the loss-of-function constraint of gnomAD v4.1.

- **Gene column**: the modes of inheritance of the gene's diseases with at least supportive
  evidence (`AD`, `AR`, `XL`, `XLR`, `SD` semidominant, `MT` ...) and its LOEUF, in red below 0.6
  (intolerant of loss of function, gnomAD's threshold). The tooltip lists every disease with its
  mode of inheritance, GenCC classification and number of submitters, plus pLI and the missense Z
  score.
- **Condition column**: the ClinVar condition of the variant, and under it the gene's
  best-supported GenCC disease (`+2`: two more in the tooltip).
- **Inheritance candidates**: next to the model, the gene's modes of inheritance in green when
  they fit the model (de novo / dominant: AD, SD, X-linked; homozygous / compound heterozygous:
  AR, SD; X-linked hemizygous: X-linked) and struck through in red when they do not, e.g. a de
  novo variant in a gene known only for recessive disease. Within a model, fitting candidates come
  first, then those of genes without GenCC data, then the others.
- **Links**: clicking a gene name opens a menu with the variant in ClinVar, dbSNP, gnomAD,
  [Franklin](https://franklin.genoox.com/) and [VarSome](https://varsome.com/) (VCF coordinates of
  the file's build), the gene in OMIM and GenCC, and **IGV**: the locus in an IGV desktop running on
  the same computer (View > Preferences > Advanced > Enable port 60151, the default).

The TSV files add `gene_moi`, `gene_diseases` (`disease (MOI, classification); ...`), `loeuf` and
`panels`; the candidates TSV adds `moi_fit` (`yes`, `no`, `.` without GenCC data); the JSON report
carries the same fields, `moi_fit` as `true` / `false` / `null`, and the sources of the gene table.
Gene names are matched through HGNC, so a gene named with a previous symbol in refGene (`AARS`)
gets the information of its approved symbol (`AARS1`).

**Gene panels** (`--gene_panel epilepsy.txt[,cardio.tsv]`): a text file with one gene symbol per
line, or a [PanelApp](https://panelapp.genomicsengland.co.uk/) / PanelApp Australia TSV download, of
which the green genes (`GEL_Status` 3–4) are used. The report names the panels at the top and
counts the variants in their genes in every section heading (`P/LP (5 · 2 in the gene panel)`);
variants in other genes are hidden behind a *show … outside the gene panel* switch above each
table, which keeps incidental findings out of view until asked for. The TSV and JSON files list
every variant, with the panels of its gene in `panels`.

## `normalized/<id>/` (only with `--save_normalized`)

Per-chunk VCFs exactly as given to ANNOVAR.

## `pipeline_info/`

| File | Content |
|---|---|
| `software_versions.yml` | versions of ANNOVAR, bcftools, htslib, Perl, Python per process, and every ANNOVAR database used per build (see below) |
| `execution_report_*.html` | Nextflow resource usage report |
| `execution_timeline_*.html` | task timeline |
| `execution_trace_*.txt` | per-task trace (CPU, memory, I/O) |
| `pipeline_dag_*.html` | workflow graph |

Databases in `software_versions.yml`, one block per genome build annotated in the run, in protocol
order. The key is the database actually used (`clinvar_latest` resolved to its release); the file
date is the ANNOVAR release date for files fetched by `annotate_variation.pl`, otherwise the date
the file was built or installed:

```yaml
"ANNOVAR databases (hg38)":
    humandb: "/data/annovar/humandb"
    "refGene": {operation: "g", file: "hg38_refGene.txt", file_date: "2019-12-05", size: "22 MB"}
    "avsnp151": {operation: "f", file: "hg38_avsnp151.txt", file_date: "2024-08-22", size: "27.6 GB"}
    "clinvar_20261025": {protocol: "clinvar_latest", operation: "f", file: "hg38_clinvar_20261025.txt", file_date: "2026-10-28", size: "878 MB"}
```
