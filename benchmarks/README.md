# Benchmarks

[`benchmark.sh`](benchmark.sh) runs annovar-nf on a VCF, then a single `table_annovar.pl --thread N`
on the same normalized records with the same protocol, and compares the two tables row by row:

```bash
benchmarks/benchmark.sh -i HG002.vcf.gz -s /opt/annovar -d /opt/annovar/humandb -t 48 -o bench/ \
    -- -profile conda --genome hg38
column -t bench/benchmark.tsv && cat bench/concordance.txt
```

`benchmark.tsv` lists, per run, the wall time, the CPU hours actually used (Nextflow trace: real
time × %CPU; `/usr/bin/time` for the baseline), the largest peak memory of one process and the
number of records. `-b` skips the baseline run.

## Setup

- **Data**: Genome in a Bottle HG002 small-variant calls, GRCh38, chromosomes 1–22
  ([`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/)):
  4,048,331 records, 4,095,611 rows after splitting multi-allelic sites. The exome-sized input is
  the same file restricted to the coding exons of refGene ±50 bp (52 Mb): 44,034 records,
  44,205 rows.
- **Protocol**: the hg38 protocol of [`conf/dragen.config`](../conf/dragen.config), 20 databases —
  refGene, knownGene, ensGene, avsnp151, 1000 Genomes ×5, COSMIC 70, ExAC, gnomAD 4.1 exome and
  genome, REVEL, NCI-60, ClinVar (NCBI 2026-10-25), InterVar 2025-07-21, dbNSFP 5.4a, regSNP-intron,
  PubMind 1.1.1.
- **Machine**: 2 × Xeon Gold 6226R (32 cores, 64 threads), 503 GB RAM, databases on NVMe; 48 CPUs
  given to each run (`--threads 48`, `-profile dragen,conda`; the baseline with `--thread 48`).

## Results

| Input | Run | Wall time | CPU hours | Peak memory of one process |
|---|---|---|---|---|
| whole genome | `table_annovar.pl --thread 48`, one run | 40 min 45 s | 9.09 | 9.1 GB |
| | **annovar-nf, 48 chunks** | **9 min 23 s** | **5.79** | **0.9 GB** |
| coding exons (exome-sized) | `table_annovar.pl --thread 48`, one run | 6 min 0 s | 2.70 | 0.9 GB |
| | **annovar-nf, 44 chunks** | **3 min 10 s** | **1.68** | **0.8 GB** |

**Concordance:** 4,095,611 of 4,095,611 and 44,205 of 44,205 rows identical.

The ANNOVAR tasks of annovar-nf ran at 99 % CPU (median); per chunk, the whole genome took a
median of 7.3 min (slowest 8.9 min), the exome-sized input 2.3 min (slowest 2.7 min). The
annovar-nf wall time includes Nextflow start-up, preparation, normalization, gathering and the
report.

`table_annovar.pl --thread N` keeps its CPUs far less busy: its Perl threads wait for each other
after every database, and it annotates genes with a single thread below one million records. On
the whole genome its 48 threads used 9.1 CPU hours in 41 minutes — fewer than a third of the
CPUs on average — and one process held 9 GB.

## Cost of one ANNOVAR process

Loading the database indexes costs every ANNOVAR process about a minute with this protocol,
whatever its number of records (`avsnp151` alone 32 s); gene annotation adds 16 s. The annotation
itself grows with the genomic range a process covers. Measured on HG002 with one thread, on 20
records (fixed cost) and on chr2:1–25,000,000 (39,638 records):

| Database | Fixed (s) | Per 25 Mb of a genome (s) |
|---|---|---|
| refGene / knownGene / ensGene | 3.8 / 8.2 / 4.4 | 2.1 / 2.5 / 2.3 |
| avsnp151 | 32.2 | 5.4 |
| 1000 Genomes (each of 5) | 3.2 | 2.5 |
| gnomAD 4.1 genome | 3.9 | 25.7 |
| gnomAD 4.1 exome | 0.7 | 4.6 |
| InterVar | 0.4 | 6.3 |
| dbNSFP 5.4a | 0.8 | 2.3 |
| PubMind | 0.3 | 0.9 |
| others (COSMIC, ExAC, REVEL, NCI-60, ClinVar, regSNP-intron) | 3.4 | 6.7 |

This is why annovar-nf plans one chunk per CPU but none below 1,000 records.
