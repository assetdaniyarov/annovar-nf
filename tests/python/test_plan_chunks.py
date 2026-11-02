"""Tests for bin/plan_chunks.py."""


def test_small_file_is_one_chunk(plan):
    rows = [("chr1", 248956422, 10), ("chr2", 242193529, 5)]
    assert plan.plan(rows, 100) == [("all", ".", 15)]


def test_chunking_disabled(plan):
    assert plan.plan([("chr1", 1000, 5_000_000)], 0) == [("all", ".", 5_000_000)]


def test_empty(plan):
    assert plan.plan([], 100) == []


def test_primary_contigs_cut_into_equal_pieces(plan):
    rows = [("chr1", 1000, 25), ("chr2", 500, 5)]
    chunks = plan.plan(rows, 10)
    assert chunks == [("chr1_1", "chr1:1-400", 10), ("chr1_2", "chr1:401-800", 10),
                      ("chr1_3-chr2", "chr1:801-1000,chr2:1-500", 10)]


def test_windows_cover_contig_without_gaps(plan):
    chunks = plan.plan([("1", 249250621, 2_664_646)], 1_000_000)
    bounds = [tuple(map(int, c[1].split(":")[1].split("-"))) for c in chunks]
    assert bounds[0][0] == 1 and bounds[-1][1] == 249250621
    assert all(b[0] == a[1] + 1 for a, b in zip(bounds, bounds[1:]))


def test_non_primary_contigs_grouped(plan):
    rows = [("chr1", 1000, 8), ("chrUn_KI270742v1", 186739, 3), ("chrEBV", 171823, 2),
            ("HLA-A*01:01:01:01", 3503, 1)]
    chunks = plan.plan(rows, 10)
    names = [c[0] for c in chunks]
    assert names == ["chr1", "other_1"]
    assert chunks[1][1] == "chrUn_KI270742v1:1-186739,chrEBV:1-171823,HLA-A*01:01:01:01:1-3503"
    # each item splits back into contig and coordinates at its last ':'
    contigs = [item.rsplit(":", 1)[0] for item in chunks[1][1].split(",")]
    assert contigs == ["chrUn_KI270742v1", "chrEBV", "HLA-A*01:01:01:01"]


def test_small_chromosomes_packed(plan):
    # A genome per chromosome would be 6 chunks; 30 records / 10 per chunk -> 3.
    rows = [("chr1", 1000, 9), ("chr2", 1000, 6), ("chr3", 1000, 5), ("chr4", 1000, 4),
            ("chrX", 1000, 4), ("chrM", 16569, 2)]
    chunks = plan.plan(rows, 10)
    assert [(c[0], c[2]) for c in chunks] == [("chr1", 9), ("chr2-chr3", 11), ("chr4-chrM", 10)]
    assert chunks[1][1] == "chr2:1-1000,chr3:1-1000"
    assert sum(c[2] for c in chunks) == 30


WGS_COUNTS = [420, 430, 350, 370, 330, 320, 300, 280, 230, 260, 260, 250, 190, 170, 160, 170, 150,
              160, 120, 110, 80, 80, 210, 3, 1]
WGS_NAMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def wgs_rows(scale=1000):
    """A DRAGEN whole genome: about 5.4 million records (typical per-contig counts)."""
    return ([(c, 100_000_000, n * scale) for c, n in zip(WGS_NAMES, WGS_COUNTS)]
            + [("chrUn_GL000220v1", 161802, 9000)])


def covered(chunks):
    """Regions of all chunks as {contig: [(start, end), ...]}."""
    out = {}
    for _, regions, _ in chunks:
        for item in regions.split(","):
            contig, span = item.rsplit(":", 1)
            out.setdefault(contig, []).append(tuple(map(int, span.split("-"))))
    return out


def test_whole_genome_by_chunk_size(plan):
    chunks = plan.plan(wgs_rows(), 1_000_000)
    assert len(chunks) == 7 and chunks[-1][0] == "other_1"
    assert all(800_000 < c[2] < 1_000_000 for c in chunks[:-1])
    assert set(covered(chunks)) == {r[0] for r in wgs_rows()}


def test_one_chunk_per_cpu(plan):
    chunks = plan.plan(wgs_rows(), 1_000_000, parts=48)
    assert len(chunks) == 48 and chunks[-1][0] == "other_1"
    sizes = [c[2] for c in chunks[:-1]]
    assert max(sizes) < 1.15 * sum(sizes) / len(sizes)
    assert chunks[0] == ("chr1_1", "chr1:1-27375887", 114979)
    # every primary chromosome covered once, without gaps
    for contig, spans in covered(chunks).items():
        if contig.startswith("chrUn"):
            continue
        assert spans[0][0] == 1 and spans[-1][1] == 100_000_000
        assert all(b[0] == a[1] + 1 for a, b in zip(spans, spans[1:]))


def test_cpus_never_exceeded(plan):
    for parts in (2, 7, 8, 16, 31, 48, 64):
        assert len(plan.plan(wgs_rows(), 1_000_000, parts=parts)) <= max(parts, 7)


def test_small_file_not_split_below_min_records(plan):
    exome = [(c, 100_000_000, n * 25) for c, n in zip(WGS_NAMES, WGS_COUNTS)]
    chunks = plan.plan(exome, 1_000_000, parts=48, min_records=10_000)
    assert len(chunks) == 13 and min(c[2] for c in chunks) > 0.9 * 10_000
    assert plan.plan(exome[:1], 1_000_000, parts=48, min_records=20_000) == [("all", ".", 10_500)]
    # default minimum: 1,000 records
    assert len(plan.plan(exome, 1_000_000, parts=48)) == 48
    assert plan.plan([("chr1", 1000, 1_999)], 1_000_000, parts=48) == [("all", ".", 1_999)]


def test_cut_snaps_to_chromosome_ends(plan):
    # 101 records: the cut at 50.5 falls 0.5 records into chr2 and moves to its start.
    rows = [("chr1", 1000, 50), ("chr2", 1000, 51)]
    assert [c[:2] for c in plan.plan(rows, 1_000_000, parts=2, min_records=10)] == [
        ("chr1", "chr1:1-1000"), ("chr2", "chr2:1-1000")]


def test_unknown_length_kept_whole(plan):
    assert plan.plan([("chr1", 0, 50), ("chr2", 0, 1)], 10)[0] == ("chr1", f"chr1:1-{plan.MAX_POS}", 50)


def test_read_stats_uses_header_lengths(plan, tmp_path):
    stats = tmp_path / "stats.tsv"
    stats.write_text("chr1\t.\t10\nchr2\t242193529\t0\n")
    header = tmp_path / "h.txt"
    header.write_text("##contig=<ID=chr1,length=248956422>\n#CHROM\n")
    rows = plan.read_stats(str(stats), plan.read_header_lengths(str(header)))
    assert rows == [("chr1", 248956422, 10)]


def test_cli(plan, tmp_path):
    stats = tmp_path / "stats.tsv"
    stats.write_text("chr1\t1000\t30\n")
    out = tmp_path / "chunks.tsv"
    assert plan.main(["--stats", str(stats), "--chunk-size", "10", "--out", str(out)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:2] == ["0", "chr1_1"] and len(lines) == 3
    assert plan.main(["--stats", str(stats), "--chunk-size", "100", "--parts", "3", "--min-records", "10",
                      "--out", str(out)]) == 0
    assert len(out.read_text().splitlines()) == 3
