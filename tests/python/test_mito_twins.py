"""Mitochondrial twin records (chrM <-> MT) in bin/format_annovar_table.py."""

import json

import pytest
from test_format_and_summary import ANNO, VCF_COLS, row, txt_header, vcf_header

TWIN = "annovar_nf_mito_twin"
ID_COL = len(ANNO) + 3 + 2          # VCF ID inside the Otherinfo columns


def vcf_text(contigs, records):
    head = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length={n}>" for c, n in contigs]
    return "\n".join(head + ["\t".join(VCF_COLS)] + ["\t".join(r) for r in records]) + "\n"


def rec(chrom, pos, ref="A", alt="G"):
    return [chrom, str(pos), ".", ref, alt, "50", "PASS", "DP=10", "GT", "0/1", "0/0"]


def twin_lines(path):
    return [line.split("\t") for line in path.read_text().splitlines() if TWIN in line]


def test_twins_hg38(fmt, tmp_path):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(vcf_text([("chr1", 248956422), ("chrM", 16569)], [rec("chr1", 100), rec("chrM", 3578, "T", "G")]))
    assert fmt.add_mito_twins(str(vcf), "hg38") == 1
    twins = twin_lines(vcf)
    assert [t[:5] for t in twins] == [["chrMT", "3578", TWIN, "T", "G"]]
    assert twins[0][5:] == rec("chrM", 3578, "T", "G")[5:]      # the rest of the record is kept


@pytest.mark.parametrize("chrom, length, expected", [
    ("chrM", 16571, 0),     # UCSC hg19 chrM: another sequence, positions differ
    ("chrM", 16569, 1),     # rCRS named chrM (e.g. DRAGEN hg19)
    ("MT", 16569, 1),       # b37
    ("MT", None, 1),        # MT without a declared length: rCRS by convention
    ("chrM", None, 0),
])
def test_twins_hg19_only_for_rcrs(fmt, tmp_path, chrom, length, expected):
    contigs = [(chrom, length)] if length else []
    vcf = tmp_path / "in.vcf"
    vcf.write_text(vcf_text(contigs, [rec(chrom, 73)]))
    assert fmt.add_mito_twins(str(vcf), "hg19") == expected


def test_twins_none_without_mito(fmt, tmp_path):
    vcf = tmp_path / "in.vcf"
    text = vcf_text([], [rec("chr1", 100)])
    vcf.write_text(text.rstrip("\n"))                 # no final newline: file left untouched
    assert fmt.add_mito_twins(str(vcf), "hg38") == 0
    assert vcf.read_text() == text.rstrip("\n")


def test_twins_cli_appends_after_last_line(fmt, tmp_path, capsys):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(vcf_text([], [rec("MT", 15, "C", "T")]).rstrip("\n"))
    assert fmt.main(["mito-twins", str(vcf), "--genome", "hg38"]) == 0
    assert capsys.readouterr().out.strip() == "1"
    lines = vcf.read_text().splitlines()
    assert lines[-2].startswith("MT\t15\t.") and lines[-1].startswith(f"M\t15\t{TWIN}\t")


def with_id(line, value):
    f = line.split("\t")
    f[ID_COL] = value
    return "\t".join(f)


def test_format_merges_and_drops_twins(fmt, tmp_path):
    rows = [
        row("chr1", 100, "A", "G", "exonic", "nonsynonymous SNV", "Benign", "0.1"),
        # chrM: gene annotation present, ClinVar and REVEL empty (the database calls it MT)
        row("chrM", 3578, "T", "G", "exonic", "nonsynonymous SNV", ".", ".", gene="MT-RNR2", af="0.2"),
        # its twin: no genes on 'MT', but the MT-named filter databases matched
        with_id(row("chrMT", 3578, "T", "G", "intergenic", ".", "Likely_pathogenic", "0.8", gene=".", af="0.9"),
                TWIN),
    ]
    txt = tmp_path / "x.hg38_multianno.txt"
    txt.write_text(txt_header() + "\n" + "\n".join(rows) + "\n")
    vcf = tmp_path / "in.vcf"
    vcf_header(vcf)
    body, stats = tmp_path / "b.tsv", tmp_path / "s.json"
    assert fmt.main(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(tmp_path / "h.tsv"),
                     "--out", str(body), "--stats", str(stats), "--mito-twins"]) == 0
    lines = [line.split("\t") for line in body.read_text().splitlines()]
    assert [f[0] for f in lines] == ["chr1", "chrM"]                 # twin dropped
    mito = dict(zip(ANNO, lines[1]))
    assert (mito["Func.refGene"], mito["Gene.refGene"]) == ("exonic", "MT-RNR2")   # gene columns kept
    assert (mito["CLNSIG"], mito["REVEL"]) == ("Likely_pathogenic", "0.8")         # empty ones filled
    assert mito["gnomad41_genome_AF"] == "0.2"                                     # non-empty ones kept
    s = json.loads(stats.read_text())
    assert s["counts"]["total"] == 2 and [v["vcf"] for v in s["plp"]] == ["chrM-3578-T-G"]


def test_vcf_merges_and_drops_twins(fmt, tmp_path):
    info = "ANNOVAR_DATE=2025-03-02;Func.refGene={f};Gene.refGene={g};CLNSIG={c};ALLELE_END"
    vcf_rows = [
        ["chrM", "3578", ".", "T", "G", "50", "PASS", info.format(f="exonic", g="MT-RNR2", c="."), "GT", "0/1"],
        ["chrMT", "3578", TWIN, "T", "G", "50", "PASS", info.format(f="intergenic", g=".", c="Likely_pathogenic"),
         "GT", "0/1"],
    ]
    annovar_vcf = tmp_path / "x.hg38_multianno.vcf"
    annovar_vcf.write_text("##fileformat=VCFv4.2\n" + "\n".join("\t".join(r) for r in vcf_rows) + "\n")
    txt = tmp_path / "x.txt"
    txt.write_text(txt_header() + "\n")
    vcf = tmp_path / "in.vcf"
    vcf_header(vcf)
    out = tmp_path / "out.vcf"
    assert fmt.main(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(tmp_path / "h.tsv"),
                     "--out", str(tmp_path / "b.tsv"), "--annovar-vcf", str(annovar_vcf),
                     "--vcf-out", str(out), "--mito-twins"]) == 0
    records = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    assert len(records) == 1 and records[0][0] == "chrM"
    assert records[0][7] == info.format(f="exonic", g="MT-RNR2", c="Likely_pathogenic")


def test_without_flag_twins_are_ordinary_rows(fmt, tmp_path):
    rows = [row("chrM", 15, "C", "T", "exonic", ".", ".", "."),
            with_id(row("chrMT", 15, "C", "T", "intergenic", ".", "Pathogenic", "."), TWIN)]
    txt = tmp_path / "x.txt"
    txt.write_text(txt_header() + "\n" + "\n".join(rows) + "\n")
    vcf = tmp_path / "in.vcf"
    vcf_header(vcf)
    body = tmp_path / "b.tsv"
    assert fmt.main(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(tmp_path / "h.tsv"),
                     "--out", str(body)]) == 0
    assert len(body.read_text().splitlines()) == 2
