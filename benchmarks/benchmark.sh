#!/usr/bin/env bash
# Benchmark annovar-nf against a single table_annovar.pl run on the same input.
#
#   benchmarks/benchmark.sh -i sample.vcf.gz -s /opt/annovar -d /opt/annovar/humandb -t 48 -o bench/ \
#       [-b] [-- <nextflow args, e.g. -profile conda>]
#
# 1. annovar-nf, timed, with --save_normalized so that both runs annotate the
#    same normalized records;
# 2. the baseline (skipped with -b): table_annovar.pl --thread <t> once on those
#    records - the usual way to run ANNOVAR on a large file - with the same protocol;
# 3. concordance: the baseline table, post-processed like annovar-nf's, compared
#    with annovar-nf's table row by row (differing cells counted per column).
#
# Writes <outdir>/benchmark.tsv (wall time, CPU hours, largest peak memory of a
# process, records) and <outdir>/concordance.txt. Needs nextflow, perl, python3.
set -euo pipefail

usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
INPUT='' ANNOVAR_SRC='' ANNOVAR_DB='' THREADS=8 OUTDIR=bench SKIP_BASELINE=false
while getopts 'i:s:d:t:o:bh' opt; do
    case $opt in
        i) INPUT=$(readlink -f "$OPTARG") ;;
        s) ANNOVAR_SRC=$(readlink -f "$OPTARG") ;;
        d) ANNOVAR_DB=$(readlink -f "$OPTARG") ;;
        t) THREADS=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        b) SKIP_BASELINE=true ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
shift $((OPTIND - 1))
[ "${1:-}" = '--' ] && shift
[ -n "$INPUT" ] && [ -n "$ANNOVAR_SRC" ] && [ -n "$ANNOVAR_DB" ] || usage 1

PIPELINE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
cd "$OUTDIR"

# ---- 1. annovar-nf ------------------------------------------------------------
start=$(date +%s)
nextflow run "$PIPELINE/main.nf" --input "$INPUT" --annovar_src "$ANNOVAR_SRC" --annovar_db "$ANNOVAR_DB" \
    --threads "$THREADS" --outdir annovar_nf --save_normalized --skip_sample_qc \
    --trace_report_suffix bench "$@"
nf_wall=$(( $(date +%s) - start ))
table=$(ls annovar_nf/annotation/*/*_multianno.txt.gz | head -n 1)
trace=annovar_nf/pipeline_info/execution_trace_bench.txt

# CPU hours and the largest peak memory of one task, from the Nextflow trace.
read -r nf_cpu nf_rss < <(awk -F'\t' '
    function sec(s,  t, n, i, u) { t = 0; n = split(s, a, " ")
        for (i = 1; i <= n; i++) { u = a[i]; if (u ~ /ms$/) t += u / 1000; else if (u ~ /s$/) t += u + 0
            else if (u ~ /m$/) t += u * 60; else if (u ~ /h$/) t += u * 3600; else if (u ~ /d$/) t += u * 86400 }
        return t }
    function mb(s) { return s ~ /GB/ ? s * 1024 : s ~ /KB/ ? s / 1024 : s + 0 }
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    { c += sec($h["realtime"]) * ($h["%cpu"] + 0) / 100; r = mb($h["peak_rss"]); if (r > m) m = r }
    END { printf "%.2f %.0f\n", c / 3600, m }' "$trace")
records=$(gzip -dc "$table" | tail -n +2 | wc -l)
printf 'run\twall_s\tcpu_h\tpeak_rss_mb\trecords\n' > benchmark.tsv
printf 'annovar-nf (--threads %s)\t%s\t%s\t%s\t%s\n' "$THREADS" "$nf_wall" "$nf_cpu" "$nf_rss" "$records" >> benchmark.tsv
"$SKIP_BASELINE" && { column -t -s $'\t' benchmark.tsv; exit 0; }

# ---- 2. baseline: one table_annovar.pl run -------------------------------------
# Same records (the normalized chunks, back in header contig order) and the
# same table_annovar.pl arguments as the annovar-nf tasks.
id=$(basename "$(dirname "$table")")
build=$(basename "$table" | sed -E 's/.*\.(hg[0-9]+)_multianno.*/\1/')
chunks=(annovar_nf/normalized/"$id"/*.norm.vcf.gz)
mkdir -p baseline
gzip -dc "${chunks[0]}" | grep '^#' > baseline/input.vcf
for c in "${chunks[@]}"; do gzip -dc "$c" | grep -v '^#' || true; done \
    | awk -F'\t' -v OFS='\t' '
        NR == FNR { if (sub(/^##contig=<ID=/, "")) { sub(/[,>].*/, ""); order[$0] = ++n } next }
        { print (($1 in order) ? order[$1] : n + 1), $0 }' baseline/input.vcf - \
    | LC_ALL=C sort -s -t $'\t' -k1,1n -k3,3n | cut -f 2- >> baseline/input.vcf
args=$(grep -h -oE -- '--(protocol|operation) +[^ ]+' work/*/*/.command.sh | sort -u)
protocol=$(echo "$args" | awk '$1 == "--protocol" {print $2}' | head -n 1)
operation=$(echo "$args" | awk '$1 == "--operation" {print $2}' | head -n 1)
[ -n "$protocol" ] && [ -n "$operation" ] || { echo "protocol of the annovar-nf run not found" >&2; exit 1; }

cd baseline
/usr/bin/time -f '%e %M %U %S' -o time.txt perl "$ANNOVAR_SRC/table_annovar.pl" input.vcf "$ANNOVAR_DB" \
    --buildver "$build" --outfile base --protocol "$protocol" --operation "$operation" \
    --nastring . --vcfinput --otherinfo --remove --onetranscript --thread "$THREADS" > annovar.log 2>&1
read -r base_wall base_rss base_user base_sys < time.txt
base_cpu=$(awk -v u="$base_user" -v s="$base_sys" 'BEGIN {print (u + s) / 3600}')

# ---- 3. concordance ---------------------------------------------------------------
python3 "$PIPELINE/bin/format_annovar_table.py" --txt "base.${build}_multianno.txt" --vcf input.vcf \
    --header-out header.tsv --stats stats.json --sample "$id" --chunk all --genome "$build" > body.tsv
cd ..
printf 'table_annovar.pl --thread %s\t%.0f\t%.2f\t%.0f\t%s\n' "$THREADS" "$base_wall" "$base_cpu" \
    "$(awk -v r="$base_rss" 'BEGIN {print r / 1024}')" "$(wc -l < baseline/body.tsv)" >> benchmark.tsv
column -t -s $'\t' benchmark.tsv

# Row-by-row comparison; differing cells are counted per column. ANNOVAR
# matches filter databases on position and alternate allele only and keeps the
# first record it reads, so a database listing one position and alt allele
# with several reference alleles can give different values depending on how
# the records are batched (see docs/output.md).
python3 - "$table" baseline/body.tsv <<'EOF' | tee concordance.txt
import collections, gzip, sys
with gzip.open(sys.argv[1], "rt") as a, open(sys.argv[2]) as b:
    header = a.readline().rstrip("\n").split("\t")
    rows = diff = 0
    cols = collections.Counter()
    for x, y in zip(a, b, strict=True):
        rows += 1
        if x != y:
            diff += 1
            x, y = x.rstrip("\n").split("\t"), y.rstrip("\n").split("\t")
            cols.update(header[i] for i in range(len(header)) if x[i] != y[i])
print(f"concordance: {rows - diff:,} of {rows:,} rows identical")
for col, n in cols.most_common():
    print(f"  {col}: {n} row(s) differ")
EOF
