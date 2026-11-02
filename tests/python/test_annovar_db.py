"""Tests for bin/annovar_db.py (no network access)."""

from pathlib import Path

import pytest
from conftest import load_script

db = load_script("annovar_db")

GENEPRED = "\t".join(["1", "NM_1", "chr1", "+", "10", "100", "20", "90", "1", "10,", "100,", "0", "G", "cmpl"]) + "\n"


def make_db(humandb, build, name, lines, idx_size=None, mrna=False):
    txt = humandb / f"{build}_{name}.txt"
    txt.write_text("".join(lines))
    if idx_size is not None:
        (humandb / f"{build}_{name}.txt.idx").write_text(f"#BIN\t1000\t{idx_size}\n1\t0\t0\t10\n")
    if mrna:
        (humandb / f"{build}_{name}Mrna.fa").write_text(">NM_1\nACGT\n")
    return txt


@pytest.mark.parametrize("pos,ref,alt,expected", [
    (100, "A", "G", (100, 100, "A", "G")),            # SNV
    (100, "AG", "A", (101, 101, "G", "-")),           # 1-bp deletion
    (100, "ACTG", "A", (101, 103, "CTG", "-")),       # 3-bp deletion
    (100, "A", "ATT", (100, 100, "-", "TT")),         # insertion
    (100, "ACG", "ATT", (101, 102, "CG", "TT")),      # block substitution
    (100, "AC", "GT", (100, 101, "AC", "GT")),        # MNV
])
def test_vcf_to_annovar(pos, ref, alt, expected):
    assert db.vcf_to_annovar(pos, ref, alt) == expected


def test_parse_entries():
    assert db.parse_entries("refGene:g, avsnp151, ensGene,clinvar_latest:f") == [
        ("refGene", "g"), ("avsnp151", "f"), ("ensGene", "g"), ("clinvar_latest", "f")]


def test_family_and_version():
    assert db.family_and_version("gnomad41_exome") == ("gnomad#_exome", (41,))
    assert db.newest_in_family("clinvar_20200316", ["clinvar_20240611", "clinvar_20221231", "cosmic70"]) == \
        "clinvar_20240611"


def test_download_name():
    assert db.download_name("AFR.sites.2015_08") == "1000g2015aug"
    assert db.download_name("gnomad41_genome") == "gnomad41_genome"


def test_newest_uses_release_dates():
    dates = {"gnomad211_exome": "20191205", "gnomad41_exome": "20240531", "gnomad40_exome": "20231127"}
    assert db.newest_in_family("gnomad41_exome", dates, dates) == "gnomad41_exome"
    rows = db.audit(Path("/nonexistent"), "hg38", [("gnomad211_exome", "f")], dates, None)
    assert rows[0].newest_official == "gnomad41_exome"
    # 1000G tables are dated by their package, so 2012_04 is not "newer" than 2015_08
    kg = {"ALL.sites.2012_04": "20191205", "1000g2015aug": "20191205"}
    assert db.newest_in_family("ALL.sites.2015_08", list(kg) + ["ALL.sites.2015_08"], kg) == "ALL.sites.2015_08"


def test_parse_avdblist():
    text = ("hg38_avsnp151.txt.gz\t20240524\t1\nhg38_avsnp151.txt.idx.gz\t20240524\t1\n"
            "hg38_1000g2015aug.zip\t20191205\t1\nhg19_other.txt.gz\t1\t1\n")
    assert db.parse_avdblist(text, "hg38") == {"avsnp151": "20240524", "1000g2015aug": "20191205"}


def test_resolve_latest(tmp_path):
    for d in ("clinvar_20240917", "clinvar_20261025", "clinvar_2026"):
        (tmp_path / f"hg38_{d}.txt").write_text("x\n")
    (tmp_path / "hg38_clinvar_20990101.txt.idx").write_text("x\n")
    assert db.resolve(tmp_path, "hg38", "clinvar_latest") == "clinvar_20261025"
    assert db.resolve(tmp_path, "hg19", "clinvar_latest") is None
    assert db.resolve(tmp_path, "hg38", "revel") == "revel"


def test_verify_ok_filter(tmp_path):
    line = "1\t100\t100\tA\tG\t0.5\n"
    make_db(tmp_path, "hg38", "revel", ["#Chr\tStart\tEnd\tRef\tAlt\tREVEL\n", line],
            idx_size=len("#Chr\tStart\tEnd\tRef\tAlt\tREVEL\n") + len(line))
    assert db.summarize_checks(db.verify_database(tmp_path, "hg38", "revel", "f")) == ("ok", "")


def test_verify_detects_truncation_and_stale_index(tmp_path):
    make_db(tmp_path, "hg38", "revel", ["1\t100\t100\tA\tG\t0.5"], idx_size=999)
    level, msg = db.summarize_checks(db.verify_database(tmp_path, "hg38", "revel", "f"))
    assert level == "error" and "newline" in msg and ".idx expects 999" in msg


def test_verify_gene_requires_mrna(tmp_path):
    make_db(tmp_path, "hg19", "refGene", [GENEPRED])
    level, msg = db.summarize_checks(db.verify_database(tmp_path, "hg19", "refGene", "g"))
    assert level == "error" and "Mrna.fa" in msg
    (tmp_path / "hg19_refGeneMrna.fa").write_text(">x\nA\n")
    assert db.summarize_checks(db.verify_database(tmp_path, "hg19", "refGene", "g"))[0] == "ok"


def test_verify_missing(tmp_path):
    assert db.verify_database(tmp_path, "hg38", "nope", "f")[0].message == "not installed"


def test_audit_actions(tmp_path):
    make_db(tmp_path, "hg38", "avsnp150", ["1\t1\t1\tA\tG\trs1\n"])
    make_db(tmp_path, "hg38", "clinvar_20240917", ["1\t1\t1\tA\tG\t1\n"])
    make_db(tmp_path, "hg38", "revel", ["1\t1\t1\tA\tG"])                        # truncated
    avdb = {"avsnp150": "1", "avsnp151": "1", "revel": "1", "gnomad41_genome": "1", "clinvar_20240611": "1"}
    entries = [("avsnp150", "f"), ("revel", "f"), ("gnomad41_genome", "f"), ("clinvar_latest", "f"),
               ("intervar_20250721", "f"), ("clinvar_20240917", "f")]
    rows = {r.entry: r for r in db.audit(tmp_path, "hg38", entries, avdb, "20261025")}
    assert rows["avsnp150"].status == "newer available: avsnp151" and rows["avsnp150"].action == "none"
    assert rows["revel"].status == "damaged" and rows["revel"].action == "redownload"
    assert rows["gnomad41_genome"].action == "download"
    assert rows["clinvar_latest"].database == "clinvar_20240917" and rows["clinvar_latest"].action == "clinvar"
    assert rows["clinvar_20240917"].action == "none"          # pinned: reported, not changed
    assert rows["intervar_20250721"].action == "manual"


def test_parse_published():
    page = """| Build | Table Name | Explanation | Date |
| ----- | ---------- | ----------- | ---- |
| hg38 | dbnsfp54a | dbNSFP version 5.4a (updated to include interpro) | 20260925 |
| hg19 | dbnsfp47a | dbNSFP version 4.7a | 20240525 |
| hg38 | `intervar_20250721` | InterVar: clinical interpretation (indels not supported)| 20250721 |
| hg19/hg38 | cosmic71, 72, ..., 80 | read here |
hg38 dbnsfp99a mentioned in prose, not in a table
| hg38 | 1000g2015aug (6 data sets) | 1000 Genomes |"""
    assert db.parse_published(page, "hg38") == {"dbnsfp54a", "intervar_20250721", "1000g2015aug"}
    assert db.parse_published(page, "hg19") == {"dbnsfp47a"}


def test_audit_suggests_only_published_versions(tmp_path):
    # avdblist also carries unannounced uploads (hg19 dbnsfp54a, intervar_20250801 without .idx):
    # only versions listed on the ANNOVAR download page count as newer releases.
    make_db(tmp_path, "hg19", "dbnsfp47a", ["1\t1\t1\tA\tG\t1\n"])
    make_db(tmp_path, "hg38", "intervar_20250721", ["1\t1\t1\tA\tG\t1\n"])
    avdb19 = {"dbnsfp47a": "20240525", "dbnsfp54a": "20261001"}
    avdb38 = {"intervar_20180118": "20191205", "intervar_20250721": "20250721", "intervar_20250801": "20250904"}
    row = db.audit(tmp_path, "hg19", [("dbnsfp47a", "f")], avdb19, None, {"dbnsfp47a"})[0]
    assert (row.status, row.newest_official) == ("ok", "dbnsfp47a")
    row = db.audit(tmp_path, "hg38", [("intervar_20250721", "f")], avdb38, None,
                   {"intervar_20180118", "intervar_20250721"})[0]
    assert (row.status, row.action, row.newest_official) == ("ok", "none", "intervar_20250721")
    # without the page (unreachable) every avdblist package counts, as before
    row = db.audit(tmp_path, "hg19", [("dbnsfp47a", "f")], avdb19, None)[0]
    assert row.status == "newer available: dbnsfp54a"


def test_audit_clinvar_built_from_ncbi_is_current(tmp_path):
    # Our NCBI build is newer than ANNOVAR's repackaged clinvar_20240611 and must not be
    # reported as outdated by ANNOVAR's release dates.
    make_db(tmp_path, "hg38", "clinvar_20240917", ["1\t1\t1\tA\tG\t1\n"])
    make_db(tmp_path, "hg38", "clinvar_20261025", ["1\t1\t1\tA\tG\t1\n"])
    avdb = {"clinvar_20240611": "20240701", "clinvar_20240917": "20241001"}
    row = db.audit(tmp_path, "hg38", [("clinvar_latest", "f")], avdb, "20261025")[0]
    assert row.database == "clinvar_20261025"
    assert (row.status, row.action, row.newest_official) == ("ok", "none", "clinvar_20261025")
    assert "built from NCBI ClinVar" in row.notes


def test_clinvar_build_and_index(tmp_path):
    vcf = [
        "##fileformat=VCFv4.1\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "1\t1500\t3\tG\tA\t.\t.\tALLELEID=30;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter\n",
        "1\t999\t1\tCAG\tC\t.\t.\tALLELEID=10;CLNSIG=Pathogenic;CLNDN=Disease_A\n",
        "1\t999\t2\tC\tT\t.\t.\tALLELEID=20;CLNSIG=Likely_benign\n",
        "MT\t3243\t4\tA\tG\t.\t.\tALLELEID=40;CLNSIG=Pathogenic\n",
        "1\t2000\t5\tA\t.\t.\t.\tALLELEID=50\n",
        "NW_1\t5\t6\tA\tG\t.\t.\tALLELEID=60\n",
    ]
    out = tmp_path / "hg38_clinvar_20990101.txt"
    n = db.write_database(db.clinvar_records(vcf), out)
    assert n == 5                                   # MT duplicated as M, '.' ALT and NW_ contig skipped
    lines = out.read_text().splitlines()
    assert lines[0].startswith("#Chr\tStart\tEnd\tRef\tAlt\tCLNALLELEID")
    body = [ln.split("\t") for ln in lines[1:]]
    # sorted by start after conversion: SNV 999 (bin 0), deletion 1000-1001 (bin 1000), SNV 1500
    assert [(r[0], r[1], r[3], r[4]) for r in body[:3]] == [("1", "999", "C", "T"), ("1", "1000", "AG", "-"),
                                                            ("1", "1500", "G", "A")]
    assert body[2][8] == "criteria_provided\\x2c_single_submitter"
    assert [r[0] for r in body[3:]] == ["M", "MT"]

    idx = (tmp_path / (out.name + ".idx")).read_text().splitlines()
    assert idx[0] == f"#BIN\t1000\t{out.stat().st_size}"
    raw = out.read_bytes()
    for entry in idx[1:]:
        chrom, b, o0, o1 = entry.split("\t")
        chunk = raw[int(o0):int(o1)].decode().splitlines()
        assert chunk and all(c.split("\t")[0] == chrom and int(c.split("\t")[1]) // 1000 * 1000 == int(b)
                             for c in chunk)
    assert db.summarize_checks(db.verify_database(tmp_path, "hg38", "clinvar_20990101", "f"))[0] == "ok"


def test_package_files():
    listing = ("hg19_dbnsfp47a.txt.gz\t20240525\t300\nhg19_dbnsfp47a.txt.idx.gz\t20240525\t20\n"
               "hg19_dbnsfp47a_interpro.txt.gz\t20240525\t9\nhg19_refGene.txt.gz\t1\t5\n"
               "hg19_refGeneMrna.fa.gz\t1\t6\nhg19_refGeneVersion.txt.gz\t1\t7\nhg19_refGeneWithVer.txt.gz\t1\t8\n"
               "hg19_1000g2015aug.zip\t1\t9\n")
    assert db.package_files(listing, "hg19", "dbnsfp47a") == {"hg19_dbnsfp47a.txt.gz": 300,
                                                              "hg19_dbnsfp47a.txt.idx.gz": 20}
    assert set(db.package_files(listing, "hg19", "refGene")) == {
        "hg19_refGene.txt.gz", "hg19_refGeneMrna.fa.gz", "hg19_refGeneVersion.txt.gz"}
    assert db.package_files(listing, "hg19", "1000g2015aug") == {}     # .zip: annotate_variation.pl
    assert db.package_files(listing, "hg19", "nope") == {}


class FakeResponse:
    """Minimal urlopen() response serving data[start:], optionally dropping after *cut* bytes."""

    def __init__(self, data, start, cut=None):
        self.status = 206 if start else 200
        self.headers = {"Content-Length": str(len(data) - start)}
        self.body, self.cut, self.sent = data[start:], cut, 0

    def read(self, n):
        if self.cut is not None and self.sent >= self.cut:
            raise ConnectionResetError("connection reset by peer")
        end = len(self.body) if self.cut is None else min(len(self.body), self.cut)
        chunk = self.body[self.sent:min(self.sent + min(n, 4), end)]
        self.sent += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_download_resumes_after_interruption(tmp_path, monkeypatch):
    data = bytes(range(256)) * 4
    ranges = []

    def fake_urlopen(request, timeout):
        rng = request.headers.get("Range", "bytes=0-")
        start = int(rng.split("=")[1].rstrip("-"))
        ranges.append(start)
        return FakeResponse(data, start, cut=300 if len(ranges) < 3 else None)   # drop twice

    monkeypatch.setattr(db.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(db.time, "sleep", lambda s: None)
    dest = tmp_path / "x.txt.gz"
    db.download_file("http://example/x.txt.gz", dest, size=len(data))
    assert dest.read_bytes() == data
    assert ranges == [0, 300, 600]          # each retry continued where the last one stopped


def test_download_gives_up_on_404(tmp_path, monkeypatch):
    def fake_urlopen(request, timeout):
        raise db.urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(db.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(OSError):
        db.download_file("http://example/missing.txt.gz", tmp_path / "m.gz", size=10)


def test_install_package_keeps_partial_archive(tmp_path, monkeypatch):
    calls = []

    def fake_fetch(staging, files):
        calls.append(sorted(p.name for p in staging.iterdir()))
        return False

    monkeypatch.setattr(db, "fetch_package", fake_fetch)
    staging = tmp_path / ".annovar-nf-staging.hg19_dbnsfp47a"
    staging.mkdir()
    (staging / "hg19_dbnsfp47a.txt.gz").write_bytes(b"partial")
    (staging / "hg19_dbnsfp47a.txt").write_text("stale unpacked table")
    files = {"hg19_dbnsfp47a.txt.gz": 100}
    assert not db.install_package(tmp_path, "hg19", "dbnsfp47a", tmp_path, False, files)
    assert calls == [["hg19_dbnsfp47a.txt.gz"]]                 # stale table dropped, archive kept
    assert (staging / "hg19_dbnsfp47a.txt.gz").read_bytes() == b"partial"


def test_download_parallel_joins_ranges(tmp_path, monkeypatch):
    data = bytes(range(256)) * 40

    def fake_urlopen(request, timeout):
        first, last = request.headers["Range"].split("=")[1].split("-")
        body = data[: int(last) + 1]
        resp = FakeResponse(body, int(first))
        resp.status = 206
        return resp

    monkeypatch.setattr(db.urllib.request, "urlopen", fake_urlopen)
    dest = tmp_path / "big.txt.gz"
    (tmp_path / "big.txt.gz.part1").write_bytes(data[2560:2600])     # partial part from an earlier run
    db.download_parallel("http://example/big.txt.gz", dest, size=len(data), connections=4, part_size=2560,
                         report_every=0.01)
    assert dest.read_bytes() == data
    assert not list(tmp_path.glob("*.part*"))


def test_missing_entries(tmp_path):
    make_db(tmp_path, "hg38", "refGene", [GENEPRED], mrna=True)
    make_db(tmp_path, "hg38", "clinvar_20261025", ["1\t1\t1\tA\tG\t1\n"])
    entries = [("refGene", "g"), ("clinvar_latest", "f"), ("pubmind_db_v111", "f"), ("ALL.sites.2015_08", "f")]
    assert db.missing_entries(tmp_path, "hg38", entries) == ["pubmind_db_v111", "ALL.sites.2015_08"]
    assert db.missing_entries(tmp_path, "hg19", entries[:2]) == ["refGene", "clinvar_latest"]


def test_install_missing(tmp_path, monkeypatch):
    """Missing ANNOVAR packages are downloaded once each, ClinVar is built, custom databases are reported."""
    make_db(tmp_path, "hg38", "refGene", [GENEPRED], mrna=True)
    listing = ("hg38_pubmind_db_v111.txt.gz\t20261028\t1000\nhg38_pubmind_db_v111.txt.idx.gz\t20261028\t10\n"
               "hg38_1000g2015aug.zip\t20150824\t2000\n")
    installed, clinvar = [], []
    monkeypatch.setattr(db, "avdblist_text", lambda build: listing)
    monkeypatch.setattr(db, "install_package",
                        lambda humandb, build, package, src, force, files: installed.append((package, files)) or True)
    monkeypatch.setattr(db, "build_clinvar", lambda humandb, build, release: clinvar.append(release) or "clinvar_x")
    entries = db.parse_entries("refGene:g,pubmind_db_v111:f,ALL.sites.2015_08:f,EUR.sites.2015_08:f,"
                               "clinvar_latest:f,my_lab_db:f")
    failed = db.install_missing(tmp_path, "hg38", entries, tmp_path, clinvar_release="20261025")
    assert failed == ["my_lab_db"]
    assert sorted(installed, key=lambda p: p[0]) == [
        ("1000g2015aug", None),                       # .zip: annotate_variation.pl, one package for both tables
        ("pubmind_db_v111", {"hg38_pubmind_db_v111.txt.gz": 1000, "hg38_pubmind_db_v111.txt.idx.gz": 10})]
    assert clinvar == ["20261025"]


def test_install_missing_reports_failed_downloads(tmp_path, monkeypatch):
    monkeypatch.setattr(db, "avdblist_text", lambda build: "hg38_revel.txt.gz\t20160603\t100\n")
    monkeypatch.setattr(db, "install_package", lambda *args: False)
    assert db.install_missing(tmp_path, "hg38", [("revel", "f")], tmp_path) == ["revel"]
    assert db.install_missing(tmp_path, "hg19", [], tmp_path) == []
