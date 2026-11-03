# annovar-nf: ANNOVAR databases

The protocols in [`conf/annovar_protocols.config`](../conf/annovar_protocols.config) (or your site
config) are the single source of truth: they define what is annotated **and** what has to be
installed. Three pipeline modes work from them:

| Mode | What it does | Writes to humandb |
|---|---|---|
| `--mode annotate` (default) | installs the protocol databases that are missing, quick pre-flight check, then annotation | only missing databases (`--auto_install_db false`: no) |
| `--mode check_db` | full audit of every protocol database | no |
| `--mode update_db` | install missing / damaged databases, build the newest ClinVar, re-audit, rebuild the gene table when older than 30 days | yes |

## First run: databases installed automatically

Annotation needs nothing but ANNOVAR itself: when protocol databases are missing from
`--annovar_db` (an empty or new directory is fine; it is created), the run installs them first —
ANNOVAR packages verified and moved in atomically, `clinvar_latest` built from the newest NCBI
release (`--clinvar_release`) — and then annotates. This happens once per genome build, only for
builds that input files actually have; the files of that build wait for it, everything else runs
in the meantime. Installed databases are never touched (updates are `--mode update_db`), and a
missing database that ANNOVAR does not distribute stops the run with a message. Two runs started on
the same new humandb do not collide: the second waits and then finds the databases installed. The
install log is `<outdir>/database_check/logs/<build>.db_install.log`.

The default protocols ([`conf/annovar_protocols.config`](../conf/annovar_protocols.config)) are a
core set for clinical germline review:

| Database | Purpose | hg38: download / on disk | hg19: download / on disk |
|---|---|---|---|
| `refGene`, `knownGene` | gene models; HGVS from refGene | 0.2 GB / 0.1 GB | 0.1 GB / 0.1 GB |
| `avsnp151` | dbSNP rs numbers | 8.3 GB / 30 GB | 8.0 GB / 29 GB |
| `gnomad41_*` / `gnomad211_*` (exome + genome) | population frequencies: the rare-variant filter of the report | 19 GB / 122 GB | 6.3 GB / 34 GB |
| `revel` | missense pathogenicity (ClinGen-calibrated for PP3/BP4) | 0.4 GB / 2.5 GB | 0.4 GB / 2.5 GB |
| `clinvar_latest` | ClinVar, newest NCBI release | 0.2 GB / 0.9 GB | 0.2 GB / 0.9 GB |
| `intervar_20250721` / `intervar_20180118` | automated ACMG classification | 1.1 GB / 30 GB | 0.5 GB / 11 GB |
| `pubmind_db_v111` (hg38) | pathogenicity reported in the literature | 0.01 GB / 0.08 GB | — |
| **total** | | **about 30 GB / 185 GB** | **about 15 GB / 75 GB** |

At a few MB/s from the ANNOVAR server (8 connections per file) the hg38 set takes one to several
hours. Without gnomAD the report cannot tell rare from common variants, so keep it in any protocol
used for clinical review. [`conf/dragen.config`](../conf/dragen.config) shows an extended set
(dbNSFP 5.4a with CADD, AlphaMissense and dozens of other scores; regSNP-intron; 1000 Genomes;
...).

```bash
# read-only report (my_site.config sets annovar_src / annovar_db, see the README)
nextflow run main.nf -profile docker -c my_site.config --mode check_db

# fix everything the report flags
nextflow run main.nf -profile docker -c my_site.config --mode update_db

# only hg38, a specific ClinVar release
nextflow run main.nf -profile docker -c my_site.config --mode update_db --db_builds hg38 --clinvar_release 20261025
```

Reports are written to `<outdir>/database_check/<build>.db_audit.{before,after}.{tsv,json,txt}`.
In an interactive terminal Nextflow's live display cannot redraw a block taller than the window,
so the run ends with every database in as many columns as the window is wide when that fits its
height, and otherwise with one line per build plus the databases that need action. The full tables
are in the `.txt` reports and are printed at the end with `-ansi-log false` (and in logs).

## What is checked

For every database of every protocol:

| Check | Problem detected |
|---|---|
| `<build>_<db>.txt` exists and is not empty | missing database |
| file ends with a newline | truncated download |
| `.idx` header size equals the table size | incomplete download or edited table (ANNOVAR would silently ignore the index) |
| large table (> 1 GB) without `.idx` | warning: annotation will be slow |
| plausible layout (≥ 5 columns; genePred for gene-based) | wrong or corrupt file |
| `<build>_<db>Mrna.fa` for gene-based databases | exonic consequences cannot be computed |
| leftover `*.txt.gz` / `*.idx.gz` | interrupted download |
| newer version in the same family at ANNOVAR | informational (e.g. `dbnsfp35a` → `dbnsfp47a`) |
| newer ClinVar release at NCBI | `clinvar_latest` is updated automatically |

The quick pre-flight check before every annotation run covers existence, `*Mrna.fa` and `.idx`
consistency, so a broken database stops the run before any work is done.

Versions within a family are compared by ANNOVAR release date, because version numbers are not
comparable across schemes (`gnomad211` = 2.1.1 is older than `gnomad41` = 4.1).

## Actions taken by `update_db`

| Audit result | Action |
|---|---|
| missing, distributed by ANNOVAR | download |
| damaged, distributed by ANNOVAR | download again and replace |
| `clinvar_latest` older than NCBI, or a ClinVar entry missing | build from the NCBI VCF |
| newer version available in the family | **none** — switching versions changes results, edit the protocol |
| missing, not distributed by ANNOVAR (custom databases you built yourself) | **none** — reported as `manual` |

Downloads go into a private staging directory inside humandb
(`.annovar-nf-staging.<build>_<package>`); the files are verified and only then moved into place, so
an interrupted or failed download never leaves a half-written database. Already installed files are
never overwritten unless the audit marked them as damaged. Parallel downloads of different databases:
`--db_download_forks` (default 2).

The ANNOVAR server throttles every connection to well under 1 MB/s, and `annotate_variation.pl`
gives up after a single 30-second stall, which makes large packages (hg19 `dbnsfp47a` is 22 GB)
practically impossible to fetch with it. Packages listed in the ANNOVAR `avdblist` (`.gz` files, and
the `.zip` archives of 1000 Genomes) are therefore downloaded and unpacked by `annovar_db.py`
itself, without `wget` or `unzip`: 32 MB byte ranges over 8 parallel connections, each retried and
resumed, with the size checked against the listing. Downloaded parts stay in the staging directory
when something fails, so rerunning `--mode update_db` (with or without `-resume`) continues where it
stopped. Packages missing from the listing are tried with `annotate_variation.pl`.

Protocol names that come from a differently named package are mapped automatically
(`ALL.sites.2015_08`, `EUR.sites.2015_08`, ... → `1000g2015aug`).

## ClinVar

NCBI publishes a new ClinVar release every week; ANNOVAR repackages one only now and then.
annovar-nf builds ANNOVAR databases directly from the NCBI VCF:

1. download `clinvar_<date>.vcf.gz` for GRCh37 (hg19) or GRCh38 (hg38) and verify its MD5;
2. convert alleles to ANNOVAR coordinates (shared leading bases stripped, `-` for empty alleles);
3. keep the columns of the 2024+ layout: `CLNALLELEID CLNDN CLNDISDB CLNREVSTAT CLNSIG`
   plus oncogenicity (`ONC*`) and somatic clinical impact (`SCI*`); commas are escaped as `\x2c`
   like in ANNOVAR's own files;
4. write `MT` records also as `M` so both `chrM` (UCSC) and `MT` (b37) VCFs are matched;
5. sort, write `<build>_clinvar_<date>.txt` and its `.idx`, verify, install atomically.

Use `clinvar_latest` in the protocol to always annotate with the newest installed release (the
ClinVar row of the run banner shows the resolved name, e.g. `clinvar_latest → clinvar_20261025`). A pinned entry such as
`clinvar_20240917` is never changed automatically.

## Gene table

The report's gene information (see [output.md](output.md#gene-information-links-and-gene-panels))
comes from `annovar_nf/genes.tsv` in humandb, which `annovar_db.py genes` builds from three open
sources (in a folder of its own, so that rebuilding it does not make `-resume` re-run the ANNOVAR
tasks):

| Source | What | Download |
|---|---|---|
| [GenCC](https://thegencc.org/) submissions (CC0) | gene–disease validity from ClinGen, Genomics England PanelApp, Orphanet, Ambry, Invitae, Illumina and others: disease, mode of inheritance, classification (Definitive, Strong, Moderate, Supportive, Limited; disputed and refuted ones are left out) | 27 MB |
| [gnomAD v4.1](https://gnomad.broadinstitute.org/) constraint metrics (CC0) | LOEUF, pLI and missense Z of the MANE Select transcript | 96 MB |
| [HGNC](https://www.genenames.org/) complete set | approved symbols, gene names, previous symbols | 17 MB |

One row per gene with GenCC or constraint data, under its approved symbol. A gene's modes of
inheritance are those of its diseases with at least supportive evidence; all diseases, with their
classification and number of GenCC submitters, are kept. The first annotation run builds the table
when it is missing (`--auto_install_db`; a failed download is only a warning, the report then has
no gene information), and `--mode update_db` rebuilds it when it is older than 30 days. The header
lines record the build date and the date of each source, and the report repeats them under
*Annotation sources*.

```bash
annovar_db.py genes --humandb /data/annovar/humandb                 # build if missing
annovar_db.py genes --humandb /data/annovar/humandb --max-age 30    # rebuild if older
annovar_db.py genes --humandb humandb --gencc gencc.tsv --constraint constraint.tsv --hgnc hgnc.txt  # offline
```

Cite GenCC (DiStefano MT, et al. *Genet Med.* 2022;24(8):1732–1742) and gnomAD constraint
(Karczewski KJ, et al. *Nature.* 2020;581:434–443) when you use them.

## Stand-alone use

The same logic is available without Nextflow:

```bash
bin/annovar_db.py audit    --humandb /opt/annovar/humandb --build hg38 \
                           --entries refGene:g,gnomad41_genome:f,clinvar_latest:f
bin/annovar_db.py verify   --humandb /opt/annovar/humandb --build hg19 --entries refGene:g,avsnp151:f
bin/annovar_db.py download --humandb /opt/annovar/humandb --build hg38 \
                           --annovar-src /opt/annovar gnomad41_genome dbnsfp47a
bin/annovar_db.py clinvar  --humandb /opt/annovar/humandb --build hg38 --date latest
```

## Notes on specific databases

- `kgXref` is a cross-reference table used by `table_annovar.pl --xreffile` for gene-based
  annotation; it is distributed by ANNOVAR for both builds (no need for UCSC downloads).
- gnomAD 4.x (`gnomad40_*`, `gnomad41_*`) exists for hg38 only; use `gnomad211_*` for hg19.
- `pubmind_db_v111` (hg38 only, 79 MB): [PubMind](https://pubmind.wglab.org/) 1.1.1 (WGLab;
  papers up to September 2026), the pathogenicity of small variants reported in the literature,
  extracted from full-text papers by a large language model. One row per variant: `PVID` (PubMind
  variant ids, `|`-separated), `pathogenicity_sum` (Pathogenic, Likely pathogenic, Uncertain
  significance, Conflicting, Benign ...), `paper_level_pathogenicity_score` (0 benign to 1
  pathogenic, averaged over papers; empty for Uncertain significance) and `confidence` (0–3, how
  much evidence was collected). It is a pointer to the literature, not a curated classification;
  see [output.md](output.md). Free for academic, non-commercial use (licence of the Children's
  Hospital of Philadelphia); cite Wang & Wang, *Nat Commun* 17, 9962 (2026),
  [doi:10.1038/s41467-026-76834-4](https://doi.org/10.1038/s41467-026-76834-4).
  ANNOVAR marks `pubmind_db_v110` as obsolete.
- "Newer available" only suggests releases listed on the
  [ANNOVAR download page](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)
  (read from its source, `WGLab/doc-ANNOVAR` `docs/user-guide/download.md` on GitHub).
  The server's file list (`avdblist`) also carries unannounced uploads — e.g. hg19 `dbnsfp54a`,
  or `intervar_20250801` without an `.idx` (every annotation would scan the whole table) — which
  are ignored. If the page is unreachable, every `avdblist` package counts.
- dbNSFP 5.4a (`dbnsfp54a`, hg38, ~25 GB download, ~95 GB unpacked, 133 columns) includes the
  InterPro domains, so `dbnsfp47a_interpro` is not needed next to it. Compared with 4.7a it adds
  e.g. MutPred2, MisFit, PHACTboost, MutFormer, MutScore, popEVE and GPN-MSA, and drops LRT,
  FATHMM, fathmm-MKL, MutPred, EVE, GenoCanyon, fitCons, LINSIGHT, SiPhy and the GTEx eQTL/sQTL
  columns.
- InterVar: `intervar_20250721` (hg38) has the same columns as `intervar_20180118`;
  `intervar_20180118` is the newest hg19 release.
- Very large hg19 databases (`cadd13`, `fathmm`, `eigen`: 250–350 GB each) dominate disk usage and
  annotation time; drop them from the protocol if you do not use them.
