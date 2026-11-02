"""Tests for bin/detect_genome_build.py."""

import json

import pytest

HG19_CONTIGS = ["##contig=<ID=chr1,length=249250621>", "##contig=<ID=chrX,length=155270560>"]
HG38_CONTIGS = ["##contig=<ID=chr1,length=248956422>", "##contig=<ID=chrX,length=156040895>"]
B37_CONTIGS = ["##contig=<ID=1,length=249250621,assembly=b37>", "##contig=<ID=2,length=243199373>"]


def run(detect, header, filename="sample.vcf.gz", expected="auto", strict=True):
    return detect.resolve(header, filename, expected, strict)


@pytest.mark.parametrize("header,build", [(HG19_CONTIGS, "hg19"), (HG38_CONTIGS, "hg38"), (B37_CONTIGS, "hg19")])
def test_contig_lengths(detect, header, build):
    assert run(detect, header)[0] == build


def test_contig_lengths_are_the_source(detect):
    build, report = run(detect, HG38_CONTIGS, "HG002.hard-filtered.vcf.gz")
    assert build == "hg38" and report["source"] == "contig_lengths"


@pytest.mark.parametrize("filename,build", [
    ("cohort_wes_hg19_dragen_78g.vcf", "hg19"),
    ("panel_dragen_436_hg38_lift.vcf.gz", "hg38"),
    ("cohort.GRCh37.vcf.gz", "hg19"),
    ("sites_only_grch38.bcf", "hg38"),
])
def test_filename_fallback(detect, filename, build):
    found, report = run(detect, ["##contig=<ID=chr1>"], filename)
    assert found == build and report["source"] == "filename" and report["warnings"]


def test_tokens_need_word_boundaries(detect):
    with pytest.raises(detect.BuildError):
        run(detect, [], "sample_hg190x.vcf")


def test_header_reference_hint(detect):
    header = ["##reference=file:///ref/hg38/hg38.fa"]
    assert run(detect, header, "x.vcf")[0] == "hg38"


def test_conflict_contig_vs_filename_is_error(detect):
    with pytest.raises(detect.BuildError, match="file name"):
        run(detect, HG19_CONTIGS, "cohort_hg38_lift.vcf.gz")


def test_conflict_downgraded_with_no_strict(detect):
    build, report = run(detect, HG19_CONTIGS, "cohort_hg38_lift.vcf.gz", strict=False)
    assert build == "hg19" and report["warnings"]


def test_explicit_build_contradicted_by_contigs(detect):
    with pytest.raises(detect.BuildError, match="contig lengths indicate 'hg19'"):
        run(detect, HG19_CONTIGS, expected="hg38")


def test_explicit_build_wins_over_filename(detect):
    build, report = run(detect, [], "x_hg19.vcf", expected="hg38")
    assert build == "hg38" and report["source"] == "explicit" and report["warnings"]


def test_no_evidence(detect):
    with pytest.raises(detect.BuildError, match="cannot determine"):
        run(detect, ["##fileformat=VCFv4.2"], "sample.vcf")


def test_cli_prints_build_and_report(detect, tmp_path, capsys):
    header = tmp_path / "header.txt"
    header.write_text("\n".join(["##fileformat=VCFv4.2", *HG38_CONTIGS, "#CHROM\tPOS"]) + "\n")
    report = tmp_path / "r.json"
    rc = detect.main(["--header", str(header), "--filename", "a.vcf.gz", "--report", str(report)])
    assert rc == 0
    assert capsys.readouterr().out.strip() == "hg38"
    assert json.loads(report.read_text())["build"] == "hg38"


def test_cli_aliases(detect, tmp_path, capsys):
    header = tmp_path / "header.txt"
    header.write_text("##fileformat=VCFv4.2\n")
    assert detect.main(["--header", str(header), "--filename", "a.vcf", "--expected", "GRCh37"]) == 0
    assert capsys.readouterr().out.strip() == "hg19"
