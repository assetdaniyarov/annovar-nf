#!/usr/bin/env python3
"""
plan_chunks.py
--------------
Split an indexed VCF into genomic chunks that can be annotated in parallel.

Input is the per-contig record count produced by `bcftools index --stats`
(columns: contig, length, n_records). Contig lengths missing from the index
('.') are taken from the VCF header when possible.

Strategy
  * the file gets one chunk per CPU (--parts) and at least
    ceil(total / chunk_size) chunks, of about equal size, but no chunk below
    --min-records unless the file itself is smaller; --chunk-size 0 or a
    single chunk: one chunk 'all' that covers the whole file.
  * the primary chromosomes (1-22, X, Y, M/MT, with or without the 'chr'
    prefix) are cut, in header order, into pieces of equal records, assuming
    records spread evenly along each chromosome; a cut close to a chromosome
    end moves to that end. Chunk names show the first and last piece:
    chr1_3-chr2_1 (end of chr1 and start of chr2), chr13-chr16 (whole
    chromosomes). The remaining contigs (alt, decoy, random, HLA, EBV ...)
    are never cut and are packed into chunks of their own (other_1 ...).

One single-threaded ANNOVAR process per CPU is the fastest layout:
table_annovar.pl --thread runs Perl threads that wait for each other after
every database (and annotate genes with one thread below a million records),
so they keep only half of their CPUs busy, while independent processes run
flat out. Each process loads the database indexes once (about a minute with
the full hg38 protocol), which --min-records keeps small next to the work.
Never more chunks than CPUs when the CPUs set the count: the extra chunks
would run as a second, mostly idle wave.

Chunks are emitted in VCF-header contig order, so concatenating the per-chunk
results in chunk order yields a coordinate-sorted file. Regions must be used
with `bcftools view --regions-overlap pos` so that an indel spanning a cut is
assigned to exactly one chunk.

Output TSV (no header): index, chunk_name, regions, n_records
    regions is '.' for the whole file, otherwise a comma-separated list of
    'contig:start-end' (chr1:1-62239105,chr2:1-242193529,...). Coordinates are
    always explicit; NORMALIZE_VCF splits each item at its last ':' and passes
    them to bcftools as a regions file, because contig names may contain ':'
    and '*' (HLA-DRB1*03:01:01:01) that bcftools cannot parse in --regions.
    VCF contig names cannot contain ','.
"""

from __future__ import annotations

import argparse
import math
import re
import sys

PRIMARY_RE = re.compile(r"^(chr)?([0-9]+|X|Y|M|MT)$", re.I)
# End coordinate for contigs whose length is unknown.
MAX_POS = 2**31 - 1
# Smallest chunk worth its own ANNOVAR process (index loading takes about a minute).
MIN_RECORDS = 1_000


def region(contig: str, length: int, start: int = 1, end: int | None = None) -> str:
    """'contig:start-end' with explicit coordinates, so that contig names
    containing ':' (HLA-A*01:01:01:01) can be split off at the last ':'."""
    return f"{contig}:{start}-{end or length or MAX_POS}"


def read_header_lengths(path: str | None) -> dict[str, int]:
    """Map contig -> length from ##contig header lines (order preserved)."""
    lengths: dict[str, int] = {}
    if not path:
        return lengths
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if not line.startswith("##contig="):
                continue
            cid = re.search(r"[<,]ID=([^,>]+)", line)
            clen = re.search(r"[<,]length=(\d+)", line)
            if cid:
                lengths[cid.group(1)] = int(clen.group(1)) if clen else 0
    return lengths


def read_stats(path: str, header_lengths: dict[str, int]) -> list[tuple[str, int, int]]:
    """Return [(contig, length, n_records)] for contigs that have records."""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[0]:
                continue
            contig, length, count = parts[0], parts[1], parts[2]
            try:
                n = int(count)
            except ValueError:
                continue
            if n <= 0:
                continue
            length_i = int(length) if length.isdigit() else header_lengths.get(contig, 0)
            rows.append((contig, length_i, n))
    return rows


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", text)


def plan(rows: list[tuple[str, int, int]], chunk_size: int, parts: int = 1,
         min_records: int = MIN_RECORDS) -> list[tuple[str, str, int]]:
    """Return [(chunk_name, regions, n_records)]."""
    total = sum(n for _, _, n in rows)
    if total == 0:
        return []
    by_size = math.ceil(total / chunk_size) if chunk_size > 0 else 1
    by_cpus = min(parts, total // max(min_records, 1)) if chunk_size > 0 else 1
    n = max(by_size, by_cpus)
    if n <= 1:
        return [("all", ".", total)]
    target = total / n

    # Runs of consecutive primary / other contigs, in header order.
    runs: list[list[tuple[str, int, int]]] = []
    for row in rows:
        if runs and is_primary(runs[-1][0][0]) == is_primary(row[0]):
            runs[-1].append(row)
        else:
            runs.append([row])
    other = {i: pack_whole(run, target) for i, run in enumerate(runs) if not is_primary(run[0][0])}
    primary = {i: sum(r[2] for r in run) for i, run in enumerate(runs) if is_primary(run[0][0])}
    # Pieces for the primary runs: what the CPUs leave after the other contigs,
    # else the records over the target; shared in proportion to records.
    n_other = sum(len(c) for c in other.values())
    n_primary = n - n_other if by_cpus > by_size else round(sum(primary.values()) / target)
    pieces = share(primary, max(n_primary, len(primary)))

    chunks: list[tuple[str, str, int]] = []
    other_idx = 0
    for i, run in enumerate(runs):
        if i in other:
            for regions, count in other[i]:
                other_idx += 1
                chunks.append((f"other_{other_idx}", regions, count))
        else:
            chunks.extend(cut(run, pieces[i], snap=0.1 * target))
    return chunks


def is_primary(contig: str) -> bool:
    return bool(PRIMARY_RE.match(contig))


def share(weights: dict[int, int], n: int) -> dict[int, int]:
    """Split n into integers proportional to weights (largest remainder), each at least 1."""
    total = sum(weights.values()) or 1
    exact = {k: n * w / total for k, w in weights.items()}
    out = {k: max(1, math.floor(x)) for k, x in exact.items()}
    for k in sorted(exact, key=lambda k: exact[k] - math.floor(exact[k]), reverse=True):
        if sum(out.values()) >= n:
            break
        out[k] += 1
    return out


def pack_whole(run: list[tuple[str, int, int]], target: float) -> list[tuple[str, int]]:
    """Whole contigs packed in order into groups of about `target` records: [(regions, n_records)]."""
    groups: list[list[tuple[str, int, int]]] = []
    for row in run:
        size = sum(r[2] for r in groups[-1]) if groups else 0
        if not groups or abs(size + row[2] - target) > abs(size - target):
            groups.append([])
        groups[-1].append(row)
    return [(",".join(region(c, length) for c, length, _ in g), sum(r[2] for r in g)) for g in groups]


def cut(run: list[tuple[str, int, int]], pieces: int, snap: float) -> list[tuple[str, str, int]]:
    """Cut consecutive primary contigs into `pieces` chunks of about equal
    records, assuming records spread evenly along each contig. A cut that
    falls within `snap` records of a contig end moves to that end, and a
    contig of unknown length is never cut."""
    total = sum(n for _, _, n in run)
    # Cut points as (contig index, position): the cut follows that position.
    cuts: list[tuple[int, int]] = []
    before = 0
    j = 1
    for i, (_, length, n) in enumerate(run):
        while j < pieces and total * j / pieces <= before + n:
            into = total * j / pieces - before
            if length <= 0:
                pos = 0 if into < n / 2 else MAX_POS
            elif into <= snap:
                pos = 0
            elif n - into <= snap:
                pos = length
            else:
                pos = min(length - 1, max(1, round(length * into / n)))
            cuts.append((i, pos))
            j += 1
        before += n

    # Segments of each chunk: (contig index, start, end), 1-based inclusive.
    chunks: list[list[tuple[int, int, int]]] = [[]]
    cuts_by_contig: dict[int, list[int]] = {}
    for i, pos in cuts:
        cuts_by_contig.setdefault(i, []).append(pos)
    for i, (_, length, _) in enumerate(run):
        end_pos = length if length > 0 else MAX_POS
        start = 1
        for pos in cuts_by_contig.get(i, []):
            if pos >= start:
                chunks[-1].append((i, start, pos))
                start = pos + 1
            chunks.append([])
        if start <= end_pos:
            chunks[-1].append((i, start, end_pos))
    chunks = [c for c in chunks if c]

    # Window labels: chr1_1, chr1_2 ... for a contig split between chunks.
    split = {}
    for chunk in chunks:
        for i, _, _ in chunk:
            split[i] = split.get(i, 0) + 1
    seen: dict[int, int] = {}
    out = []
    for chunk in chunks:
        labels, regions, count = [], [], 0
        for i, start, end in chunk:
            contig, length, n = run[i]
            seen[i] = seen.get(i, 0) + 1
            labels.append(safe_name(contig) + (f"_{seen[i]}" if split[i] > 1 else ""))
            regions.append(region(contig, length, start, end if length > 0 else None))
            count += n if split[i] == 1 or length <= 0 else round(n * (end - start + 1) / length)
        name = labels[0] if len(labels) == 1 else f"{labels[0]}-{labels[-1]}"
        out.append((name, ",".join(regions), count))
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--stats", required=True, help="output of 'bcftools index --stats'")
    p.add_argument("--header", help="VCF header, used to fill in missing contig lengths")
    p.add_argument("--chunk-size", type=int, default=1_000_000,
                   help="at most about this many records per chunk; 0 disables splitting [1000000]")
    p.add_argument("--parts", type=int, default=1, help="CPUs for the file: one chunk each [1]")
    p.add_argument("--min-records", type=int, default=MIN_RECORDS,
                   help=f"smallest chunk made to use more CPUs [{MIN_RECORDS}]")
    p.add_argument("--out", required=True, help="output TSV")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    rows = read_stats(args.stats, read_header_lengths(args.header))
    chunks = plan(rows, args.chunk_size, args.parts, args.min_records)
    with open(args.out, "w", encoding="utf-8") as out:
        for i, (name, regions, n) in enumerate(chunks):
            out.write(f"{i}\t{name}\t{regions}\t{n}\n")
    total = sum(n for _, _, n in rows)
    print(f"[plan_chunks] {total} records on {len(rows)} contigs -> {len(chunks)} chunk(s)", file=sys.stderr)
    if not chunks:
        print("[plan_chunks] WARNING: the VCF contains no records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
