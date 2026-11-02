"""Tests for bin/format_annovar_table.py and bin/summarize_annotations.py."""

import gzip
import json
import re

import pytest

VCF_COLS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"]
ANNO = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
        "CLNSIG", "REVEL", "gnomad41_genome_AF"]


def vcf_header(path, gz=False):
    text = "##fileformat=VCFv4.2\n" + "\t".join(VCF_COLS) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


def row(chr_, pos, ref, alt, func, exonic, clnsig, revel, filt="PASS", af=".", gts=("0/1", "0/0"), gene="G"):
    other = ["0.5", "50", "10", chr_, str(pos), ".", ref, alt, "50", filt, "DP=10", "GT", *gts]
    return "\t".join([chr_, str(pos), str(pos), ref, alt, func, gene, exonic, clnsig, revel, af, *other])


def txt_header():
    return "\t".join(ANNO + [f"Otherinfo{i}" for i in range(1, 3 + len(VCF_COLS) + 1)])


@pytest.fixture
def annovar_txt(tmp_path):
    header = txt_header()
    rows = [
        row("chr1", 100, "A", "G", "exonic", "nonsynonymous SNV", "Pathogenic", "0.9"),
        row("chr1", 200, "C", "T", "intronic", ".", "Conflicting_classifications_of_pathogenicity", "."),
        row("chr1", 300, "-", "AT", "exonic", "frameshift insertion", "Likely_pathogenic", ".", filt="LowQual"),
        row("chr2", 400, "G", "T", "splicing", ".", "Benign", "0.2", af="0.3"),
        # reference call at a known pathogenic site (e.g. DRAGEN targeted caller): not carried
        row("chr2", 500, "C", "A", "exonic", "nonsynonymous SNV", "Pathogenic", ".", filt="TargetedWT",
            gts=("0/0", "./.")),
        row("chr3", 600, "G", "A", "exonic", "nonsynonymous SNV", "Uncertain_significance", "0.6", af="0.0001",
            gts=("1/1", "0/1")),
    ]
    path = tmp_path / "x.hg38_multianno.txt"
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return path


def test_build_header_new_annovar(fmt):
    txt = ANNO + [f"Otherinfo{i}" for i in range(1, 15)]
    new = fmt.build_header(txt, VCF_COLS, len(txt))
    assert new == ANNO + ["Otherinfo1", "Otherinfo2", "Otherinfo3"] + VCF_COLS


def test_build_header_old_annovar_single_otherinfo(fmt):
    # Pre-2019 ANNOVAR: one 'Otherinfo' header column, wider data rows.
    new = fmt.build_header(ANNO + ["Otherinfo"], VCF_COLS, len(ANNO) + 3 + len(VCF_COLS))
    assert new[-len(VCF_COLS):] == VCF_COLS and new[len(ANNO)] == "Otherinfo1"


def test_build_header_mismatch(fmt):
    with pytest.raises(ValueError, match="does not match"):
        fmt.build_header(ANNO + ["Otherinfo1"], VCF_COLS, len(ANNO) + 2)


def test_build_header_requires_otherinfo(fmt):
    with pytest.raises(ValueError, match="Otherinfo"):
        fmt.build_header(ANNO, VCF_COLS, len(ANNO))


@pytest.mark.parametrize("value,expected", [
    ("Pathogenic", True), ("Likely_pathogenic", True), ("Pathogenic/Likely_pathogenic", True),
    ("Conflicting_classifications_of_pathogenicity", False), ("Benign", False), (".", False),
])
def test_is_plp(fmt, value, expected):
    assert fmt.is_plp(value) is expected


@pytest.mark.parametrize("value,expected", [
    ("Uncertain_significance", "VUS"), ("Uncertain_significance|association", "VUS"),
    ("Conflicting_classifications_of_pathogenicity", "conflicting"),
    ("Conflicting_classifications_of_pathogenicity|risk_factor", "conflicting"),
    ("Uncertain_risk_allele", None), ("Pathogenic", None), (".", None),
])
def test_vus_kind(fmt, value, expected):
    assert fmt.vus_kind(value) == expected


def test_genotype_info(fmt):
    het = fmt.genotype_info("GT:AD:DP:GQ", ["0/1:30,10:40:50"], ["S1"])
    call = het["calls"][0]
    assert het["carried"] and (call["sample"], call["zyg"], call["vaf"], call["dp"], call["gq"]) == \
        ("S1", "het", 0.25, 40, 50)
    assert call["alt_reads"] == 10 and call["flags"] == [] and call["conf"] == "high"
    assert fmt.genotype_info("GT:AD", ["1|1:0,22"])["calls"][0]["zyg"] == "hom"
    hemi = fmt.genotype_info("GT:AF", ["1:0.97"])["calls"][0]
    assert (hemi["zyg"], hemi["vaf"]) == ("hemi", 0.97)
    assert fmt.genotype_info("GT:AD:DP", ["0/0:35,0:35"])["carried"] is False
    assert fmt.genotype_info("GT", ["./."])["carried"] is False
    # families: every call kept; cohorts: carriers only
    cohort = fmt.genotype_info("GT", ["0/1", "1/1", "0/0", "./."], ["A", "B", "C", "D"])
    assert (cohort["het"], cohort["hom"], [c["zyg"] for c in cohort["calls"]]) == \
        (1, 1, ["het", "hom", "ref", "missing"])
    carriers = fmt.genotype_info("GT", ["0/1", "1/1", "0/0"], ["A", "B", "C"], keep_all=False)
    assert [c["sample"] for c in carriers["calls"]] == ["A", "B"]
    assert fmt.genotype_info("", [])["carried"] is None          # sites-only VCF


@pytest.mark.parametrize("sample,fmt_keys,site_pass,flags,conf", [
    ("0/1:20,18:38:60", "GT:AD:DP:GQ", True, [], "high"),
    ("0/1:4,3:7:30", "GT:AD:DP:GQ", True, ["LOW_DP", "LOW_ALT_READS"], "low"),
    ("0/1:40,5:45:5", "GT:AD:DP:GQ", True, ["LOW_GQ", "VAF_SKEW"], "low"),
    ("0/1:30,6:36:40", "GT:AD:DP:GQ", True, ["VAF_SKEW"], "caution"),
    ("1/1:10,20:30:40", "GT:AD:DP:GQ", True, ["VAF_SKEW"], "caution"),
    ("0/1:20,18:38:60:LowGQ:1:0", "GT:AD:DP:GQ:FT:HML:MSC", True, ["FT:LowGQ", "HOMOLOGY"], "low"),
    ("0/1:20,18:38:60:PASS:0:1", "GT:AD:DP:GQ:FT:HML:MSC", True, ["MOSAIC"], "caution"),
    ("0/0:30,3:33:60", "GT:AD:DP:GQ", True, ["ALT_READS_IN_REF"], "caution"),
    ("0/1:20,18:38:60", "GT:AD:DP:GQ", False, [], "low"),          # site FILTER not PASS
    ("./.:0,0:0:0", "GT:AD:DP:GQ", True, [], "missing"),
])
def test_call_quality_flags(fmt, sample, fmt_keys, site_pass, flags, conf):
    keys = {k: i for i, k in enumerate(fmt_keys.split(":"))}
    call = fmt.parse_call(keys, sample, fmt.Qc(), site_pass)
    assert (call["flags"], call["conf"]) == (flags, conf)


def test_call_dragen_de_novo(fmt):
    keys = {k: i for i, k in enumerate("GT:AD:DP:GQ:DN:DQ".split(":"))}
    call = fmt.parse_call(keys, "0/1:20,18:38:60:DeNovo:25.3", fmt.Qc())
    assert (call["dn"], call["dq"]) == ("DeNovo", 25.3)
    assert "dn" not in fmt.parse_call(keys, "0/1:20,18:38:60:.:.", fmt.Qc())


def test_call_mitochondrial_heteroplasmy(fmt):
    keys = {"GT": 0, "AD": 1, "DP": 2, "GQ": 3}
    # 4% heteroplasmy at 1026x: a skewed fraction on a nuclear chromosome, expected on chrM
    assert fmt.parse_call(keys, "0/1:981,45:1026:60", fmt.Qc())["flags"] == ["VAF_SKEW"]
    call = fmt.parse_call(keys, "0/1:981,45:1026:60", fmt.Qc(), mito=True)
    assert (call["flags"], call["conf"], call["vaf"]) == ([], "high", 0.044)


def test_splicing_is_protein_coding_only(fmt):
    assert fmt.splicing("splicing") and fmt.splicing("exonic;splicing")
    assert not fmt.splicing("ncRNA_splicing") and not fmt.splicing("intronic")


def test_call_quality_thresholds(fmt):
    keys = {"GT": 0, "AD": 1, "DP": 2}
    strict = fmt.Qc(min_dp=50)
    assert fmt.parse_call(keys, "0/1:20,18:38", strict)["flags"] == ["LOW_DP"]
    args = fmt.parse_args(["--txt", "t", "--vcf", "v", "--header-out", "h", "--qc-het-vaf", "0.3,0.7"])
    assert (args.qc.min_het_vaf, args.qc.max_het_vaf) == (0.3, 0.7)


@pytest.mark.parametrize("gz", [False, True])
def test_format_end_to_end(fmt, tmp_path, annovar_txt, gz):
    vcf = tmp_path / ("in.vcf.gz" if gz else "in.vcf")
    vcf_header(vcf, gz)
    header_out, body_out, stats_out = tmp_path / "h.tsv", tmp_path / "b.tsv", tmp_path / "s.json"
    rc = fmt.main(["--txt", str(annovar_txt), "--vcf", str(vcf), "--header-out", str(header_out),
                   "--out", str(body_out), "--stats", str(stats_out), "--sample", "S", "--chunk", "chr1",
                   "--genome", "hg38"])
    assert rc == 0

    header = header_out.read_text().rstrip("\n").split("\t")
    assert header[-len(VCF_COLS):] == VCF_COLS
    body = body_out.read_text().splitlines()
    assert len(body) == 6 and all(len(line.split("\t")) == len(header) for line in body)

    stats = json.loads(stats_out.read_text())
    c = stats["counts"]
    assert (c["total"], c["snv"], c["indel"], c["ti"], c["tv"]) == (6, 5, 1, 3, 2)
    assert c["pass"] == 4 and c["revel_ge_0.75"] == 1
    assert stats["func"]["exonic\tnonsynonymous SNV"] == 3
    assert stats["n_genotype_columns"] == 2 and stats["genome"] == "hg38"
    # rare = AF < 1% or absent from gnomAD; coding = protein-altering or splicing
    assert c["rare"] == 5 and c["coding_rare"] == 4
    # P/LP: only carried records are listed; the 0/0 reference call is counted apart
    assert (c["clinvar_plp"], c["clinvar_plp_not_carried"], c["vus_candidates"]) == (2, 1, 1)
    assert [v["clnsig"] for v in stats["plp"]] == ["Pathogenic", "Likely_pathogenic"]
    first = stats["plp"][0]
    assert (first["source"], first["exonic"], first["vcf"], first["filter"]) == \
        ("ClinVar", "nonsynonymous SNV", "chr1-100-A-G", "PASS")
    assert first["carriers"] == "1 het / 0 hom"
    # two samples (a family): every genotype listed
    assert [(g["sample"], g["zyg"]) for g in first["genotypes"]] == [("S1", "het"), ("S2", "ref")]
    assert stats["plp"][1]["filter"] == "LowQual"
    assert stats["plp"][1]["genotypes"][0]["conf"] == "low"          # site filtered by the caller
    # VUS: rare exonic, carried; the intronic conflicting record is not listed
    assert [(v["kind"], v["gene"], v["revel"], v["gnomad_af"]) for v in stats["vus"]] == \
        [("VUS", "G", "0.6", "0.0001")]
    assert stats["vus"][0]["carriers"] == "1 het / 1 hom"
    # rare coding/splicing variants carried by a sample, for the inheritance analysis
    assert [v["vcf"] for v in stats["rare"]] == ["chr1-100-A-G", "chr1-300---AT", "chr3-600-G-A"]
    assert stats["samples"] == ["S1", "S2"] and stats["qc"]["min_dp"] == 10
    # per sample, counted like the totals: S2 carries only the rare VUS (1/1 is S1, 0/1 is S2)
    assert stats["sample_counts"] == {
        "S1": {"coding_rare": 3, "clinvar_plp": 2, "vus_candidates": 1},
        "S2": {"coding_rare": 1, "vus_candidates": 1}}


def test_carriers_read_gt_only(fmt):
    header = ANNO + ["Otherinfo1", "Otherinfo2", "Otherinfo3"] + VCF_COLS[:-2] + ["A", "B", "C", "D"]
    stats = fmt.Stats(header)
    line = row("chr1", 5, "A", "G", "exonic", "stopgain", ".", ".", gts=("0|1", "./.", "1", "0/0"))
    f = line.split("\t", stats.max_idx + 1)
    assert stats.carriers(f) == ["A", "C"]
    f[stats.idx["format"]] = "DP:GT"                     # GT anywhere in FORMAT
    f[-1] = f[-1].replace("0|1", "9:0|1").replace("./.", "9:./.")
    assert stats.carriers(f) == ["A"]
    f[stats.idx["format"]] = "DP"                        # no genotypes
    assert stats.carriers(f) == []


def test_gene_separator(fmt):
    header = ["Chr", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "AAChange.refGene", "Gene_ensGene"]
    cols = fmt.gene_columns(header)
    assert cols == [2, 3, 5]
    line = "chr1\tupstream;downstream\tA,B;C\tdist=1,dist=2;dist=3\tA:NM_1:c.1A>G,A:NM_2:c.1A>G\tX,Y\tS1,S2\n"
    assert fmt.gene_separator(line, cols) == \
        "chr1\tupstream;downstream\tA;B;C\tdist=1;dist=2;dist=3\tA:NM_1:c.1A>G,A:NM_2:c.1A>G\tX;Y\tS1,S2\n"
    single = "chr1\texonic\tA\t.\tA:NM_1:c.1A>G,A:NM_2:c.1A>G\tX\n"
    assert fmt.gene_separator(single, cols) is single
    assert fmt.gene_separator(line, []) is line


def test_vcf_gene_separator(fmt):
    info = (r"DP=10;Func.refGene=intergenic;Gene.refGene=SHOX,CRLF2;GeneDetail.refGene=dist\x3d1,dist\x3d2;"
            "AAChange.refGene=A:B,C:D;ALLELE_END")
    line = "\t".join(["chrX", "1", ".", "A", "G", "50", "PASS", info, "GT", "0/1"]) + "\n"
    out = fmt.vcf_gene_separator(line).split("\t")
    assert r"Gene.refGene=SHOX\x3bCRLF2;" in out[7] and r"GeneDetail.refGene=dist\x3d1\x3bdist\x3d2;" in out[7]
    assert "AAChange.refGene=A:B,C:D;" in out[7] and out[8:] == ["GT", "0/1\n"]
    assert fmt.vcf_gene_separator("#CHROM\tPOS\n") == "#CHROM\tPOS\n"
    sites_only = "\t".join(["chr1", "1", ".", "A", "G", ".", "PASS", "Gene.refGene=A,B"]) + "\n"
    assert fmt.vcf_gene_separator(sites_only).endswith("Gene.refGene=A\\x3bB\n")


def test_format_gene_separator_table_and_vcf(fmt, tmp_path):
    txt = tmp_path / "x.hg38_multianno.txt"
    line = row("chrX", 100, "A", "G", "intergenic", ".", ".", ".", gene="SHOX,CRLF2")
    txt.write_text(txt_header() + "\n" + line + "\n")
    vcf = tmp_path / "in.vcf"
    vcf_header(vcf)
    annovar_vcf = tmp_path / "x.hg38_multianno.vcf"
    annovar_vcf.write_text("##fileformat=VCFv4.2\n" + "\t".join(VCF_COLS) + "\n" + "\t".join(
        ["chrX", "100", ".", "A", "G", "50", "PASS", "Gene.refGene=SHOX,CRLF2;ALLELE_END", "GT", "0/1", "0/0"]) + "\n")
    h, b, v = tmp_path / "h.tsv", tmp_path / "b.tsv", tmp_path / "out.vcf"
    assert fmt.main(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(h), "--out", str(b),
                     "--annovar-vcf", str(annovar_vcf), "--vcf-out", str(v)]) == 0
    assert b.read_text().split("\t")[6] == "SHOX;CRLF2"
    assert v.read_text().splitlines()[1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"
    assert r"Gene.refGene=SHOX\x3bCRLF2;ALLELE_END" in v.read_text()
    with pytest.raises(SystemExit):
        fmt.parse_args(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(h), "--annovar-vcf", "x.vcf"])


def test_format_empty_table(fmt, tmp_path):
    txt = tmp_path / "empty.txt"
    txt.write_text("\t".join(ANNO + ["Otherinfo1"]) + "\n")
    vcf = tmp_path / "in.vcf"
    vcf_header(vcf)
    h, b = tmp_path / "h.tsv", tmp_path / "b.tsv"
    assert fmt.main(["--txt", str(txt), "--vcf", str(vcf), "--header-out", str(h), "--out", str(b)]) == 0
    assert b.read_text() == ""


def test_summary(summarize, tmp_path):
    chunks = [
        {"sample": "A", "chunk": "chr1", "n_genotype_columns": 2,
         "genome": "hg38", "counts": {"total": 3, "snv": 2, "indel": 1, "ti": 1, "tv": 1, "pass": 3, "clinvar_plp": 1,
                                      "clinvar_plp_not_carried": 173, "vus_candidates": 1},
         "func": {"exonic\tnonsynonymous SNV": 2, "intronic\t.": 1}, "columns": {"clnsig": "CLNSIG"},
         "samples": ["A1"], "qc": {"min_dp": 10, "min_gq": 10},
         "plp": [{"source": "ClinVar", "chr": "chr1", "start": "100", "ref": "A", "alt": "G", "vcf": "chr1-100-A-G",
                  "gene": "TTN", "func": "exonic", "exonic": "nonsynonymous SNV", "aachange": "TTN:NM_1:p.X1Y",
                  "clnsig": "Pathogenic", "clnrevstat": "reviewed_by_expert_panel", "clndn": "DCM",
                  "clnalleleid": "12345", "filter": "PASS", "gnomad_af": ".", "revel": "0.9",
                  "carriers": "1 het / 0 hom",
                  "genotypes": [{"sample": "A1", "gt": "0/1", "zyg": "het", "vaf": 0.48, "dp": 35, "gq": 60,
                                 "alt_reads": 17, "flags": [], "conf": "high"}]}],
         "vus": [{"kind": "VUS", "chr": "chr1", "start": "900", "ref": "C", "alt": "T", "vcf": "chr1-900-C-T",
                  "gene": "MYH7", "func": "exonic", "exonic": "nonsynonymous SNV", "aachange": "MYH7:p.R1W",
                  "clnsig": "Uncertain_significance", "clnrevstat": "criteria_provided,_single_submitter",
                  "clndn": "HCM", "filter": "PASS", "gnomad_af": "1e-05", "revel": "0.7",
                  "genotypes": [{"sample": "A1", "gt": "0/1", "zyg": "het", "vaf": 0.3, "dp": 12, "gq": 8,
                                 "alt_reads": 4, "flags": ["LOW_GQ"], "conf": "caution"}]}],
         "rare": [{"chr": "chr1", "start": "900", "ref": "C", "alt": "T", "vcf": "chr1-900-C-T", "gene": "MYH7",
                   "func": "exonic", "exonic": "stopgain", "gnomad_af": "1e-05",
                   "genotypes": [{"sample": "A1", "zyg": "hom", "conf": "high", "flags": []}]}]},
        {"sample": "A", "chunk": "chr2", "counts": {"total": 1, "snv": 1, "ti": 1}, "func": {"intergenic\t.": 1},
         "columns": {}},
        {"sample": "B", "chunk": "all", "counts": {"total": 0}, "func": {}, "columns": {}},
    ]
    paths = []
    for i, c in enumerate(chunks):
        p = tmp_path / f"{i}.json"
        p.write_text(json.dumps(c))
        paths.append(str(p))
    build = tmp_path / "A.build.json"
    build.write_text(json.dumps({"filename": "A.hard-filtered.vcf.gz", "expected": "auto", "build": "hg38",
                                 "source": "contig_lengths", "warnings": ["mind the gap"],
                                 "evidence": {"contig_length_votes": {"hg19": 0, "hg38": 24},
                                              "header_hints": ["hg38"], "filename_tokens": []}}))
    databases = tmp_path / "databases.tsv"
    databases.write_text("hg38\trefGene\trefGene\tg\thg38_refGene.txt\t2019-12-05\t23068672\n"
                         "hg38\tclinvar_latest\tclinvar_20261025\tf\thg38_clinvar_20261025.txt\t2026-10-28\t920000000\n"
                         "hg38\tgnomad41_genome\tgnomad41_genome\tf\thg38_gnomad41_genome.txt\t2024-05-01\t-1\n"
                         "hg19\tdbnsfp47a\tdbnsfp47a\tf\thg19_dbnsfp47a.txt\t2024-06-01\t1\n")
    sample_qc = tmp_path / "A.sample_qc.json"
    sample_qc.write_text(json.dumps({"id": "A", "samples": [
        {"sample": "A1", "sex_inferred": "female", "x_het_hom_ratio": 1.3, "sex_check": "mismatch",
         "mean_depth": 35.2, "snv_het": 100, "snv_hom": 60, "het_hom_ratio": 1.667, "ts_tv": 2.1,
         "ped": {"family": "F1", "father": None, "mother": None, "sex": "male", "affected": True}}],
        "pairs": [], "trios": []}))
    out = tmp_path / "summary"
    assert summarize.main(["--stats", *paths, "--builds", str(build), "--out", str(out),
                           "--databases", str(databases), "--annovar-version", "2025-03-02",
                           "--sample-qc", str(sample_qc)]) == 0

    # combined report (more than one input file)
    totals = (tmp_path / "summary_totals.tsv").read_text().splitlines()
    header = totals[0].split("\t")
    a = dict(zip(header, totals[1].split("\t")))
    b = dict(zip(header, totals[2].split("\t")))
    assert a["total"] == "4" and a["ti_tv"] == "2.000" and a["clinvar_plp"] == "1" and a["exonic"] == "2"
    assert b["total"] == "0" and b["clinvar_plp"] == "NA"
    assert a["genome"] == "hg38" and b["genome"] == "NA"

    matrix = (tmp_path / "summary.tsv").read_text().splitlines()
    assert matrix[0] == "Func.refGene\tExonicFunc.refGene\tA\tB"
    assert matrix[1].startswith("exonic\tnonsynonymous SNV\t2\t0")
    page = (tmp_path / "summary.html").read_text(encoding="utf-8")
    assert "<table" in page and "TTN" in page and "★★★☆" in page
    plp = (tmp_path / "summary_plp.tsv").read_text(encoding="utf-8").splitlines()
    cells = dict(zip(plp[0].split("\t"), plp[1].split("\t")))
    assert len(plp) == 2 and (cells["file"], cells["sample"], cells["gene"], cells["stars"], cells["zygosity"],
                              cells["vaf"], cells["confidence"], cells["source"]) == \
        ("A", "A1", "TTN", "3", "het", "0.48", "high", "ClinVar")
    vus = (tmp_path / "summary_vus.tsv").read_text(encoding="utf-8").splitlines()
    vus_cells = dict(zip(vus[0].split("\t"), vus[1].split("\t")))
    assert len(vus) == 2 and (vus_cells["gene"], vus_cells["flags"], vus_cells["confidence"]) == \
        ("MYH7", "LOW_GQ", "caution")
    # single affected sample without parents: homozygous rare stopgain is an AR candidate
    cand = (tmp_path / "A.summary_candidates.tsv").read_text(encoding="utf-8").splitlines()
    assert len(cand) == 2 and dict(zip(cand[0].split("\t"), cand[1].split("\t")))["model"] == "homozygous (AR)"
    assert "homozygous (AR)" in (tmp_path / "A.summary_vus.tsv").read_text(encoding="utf-8")
    samples = (tmp_path / "A.summary_samples.tsv").read_text(encoding="utf-8").splitlines()
    srow = dict(zip(samples[0].split("\t"), samples[1].split("\t")))
    assert (srow["sample"], srow["sex_declared"], srow["sex_inferred"], srow["sex_check"], srow["affected"]) == \
        ("A1", "male", "female", "mismatch", "affected")
    report = json.loads((tmp_path / "A.summary.json").read_text(encoding="utf-8"))
    assert report["schema"] == "annovar-nf.report" and report["annovar"]["version"] == "2025-03-02"
    assert report["pathogenic"][0]["genotypes"][0]["sample"] == "A1"
    assert report["inheritance_candidates"][0]["model"] == "homozygous (AR)"
    assert not (tmp_path / "summary.json").exists()           # no JSON for the combined report
    # samples, P/LP, candidates and VUS come before the functional sections; reference calls are explained
    assert page.index("<h2>Samples") < page.index("Pathogenic / likely pathogenic") \
        < page.index("Inheritance candidates") < page.index("Variants of uncertain significance") \
        < page.index("Genomic regions")
    assert "✘ mismatch" in page and "LOW_GQ" in page and "call caution" in page
    assert "173 ClinVar P/LP record(s) with a reference genotype" in page
    assert "clinvar/?term=12345[alleleid]" in page and "gnomad.broadinstitute.org/variant/1-100-A-G" in page

    # one report per input file, named after it, with the genome-build evidence
    page_a = (tmp_path / "A.summary.html").read_text(encoding="utf-8")
    assert "<title>A · annovar-nf annotation report</title>" in page_a
    assert "##contig lengths" in page_a and "hg38 24 · hg19 0" in page_a and "mind the gap" in page_a
    assert "TTN" in page_a
    # database versions of the file's build: in the header and as a table
    assert "<p>ClinVar 2026-10-25 · gnomAD 4.1 · ANNOVAR 2025-03-02</p>" in page_a
    assert page_a.index("Exonic consequences") < page_a.index("Annotation sources")
    assert "clinvar_20261025" in page_a and "(clinvar_latest)" in page_a and "gnomAD genomes" in page_a
    assert "missing" in page_a and "dbnsfp47a" not in page_a        # hg19 database: other build
    # numeric columns: header right-aligned like the cells
    assert "<th class='num'>file date</th><th class='num'>size</th>" in page_a
    assert "<th class='num'>mean depth</th>" in page_a and "<th>sex from chrX</th>" in page_a
    assert "<th class='num'>total</th>" in page_a and "<th>genome</th>" in page_a
    assert re.search(r"<th class='num'><span title='gnomAD[^>]*>gnomAD</span></th><th class='num'>REVEL</th>", page_a)
    # statistics without per-sample counts (an older version): NA, not 0
    srow_a = dict(zip(samples[0].split("\t"), samples[1].split("\t")))
    assert (srow_a["coding_rare"], srow_a["clinvar_plp"]) == ("NA", "NA")
    assert (tmp_path / "A.summary_totals.tsv").read_text().splitlines()[1].startswith("A\thg38\t")
    assert len((tmp_path / "B.summary_plp.tsv").read_text().splitlines()) == 1
    assert "TTN" not in (tmp_path / "B.summary.html").read_text(encoding="utf-8")


@pytest.mark.parametrize("status,stars", [
    ("practice_guideline", 4), ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
    ("criteria_provided,_conflicting_classifications", 1), ("criteria_provided,_single_submitter", 1),
    ("no_assertion_criteria_provided", 0), (".", 0),
])
def test_review_stars(summarize, status, stars):
    assert summarize.review_stars(status) == stars


@pytest.mark.parametrize("name,expected", [
    ("clinvar_20261025", ("ClinVar", "2026-10-25")), ("gnomad41_genome", ("gnomAD genomes", "4.1")),
    ("gnomad211_exome", ("gnomAD exomes", "2.1.1")), ("dbnsfp54a", ("dbNSFP", "5.4a")),
    ("dbnsfp47a_interpro", ("dbNSFP InterPro domains", "4.7a")), ("intervar_20250721", ("InterVar", "2025-07-21")),
    ("avsnp151", ("dbSNP", "151")), ("ALL.sites.2015_08", ("1000 Genomes (ALL)", "2015-08")),
    ("refGene", ("RefSeq genes (UCSC refGene)", "")), ("my_custom_db", ("my_custom_db", "")),
    ("pubmind_db_v111", ("PubMind", "1.1.1")),
])
def test_describe_database(summarize, name, expected):
    assert summarize.describe_database(name) == expected


@pytest.mark.parametrize("annovar,expected", [
    ("c.G590A", "c.590G>A"), ("c.4062_4063insCT", "c.4062_4063insCT"), ("c.590_592delCTT", "c.590_592del"),
    ("c.590delG", "c.590del"), ("c.123_124delinsAT", "c.123_124delinsAT"), ("c.194+2T>C", "c.194+2T>C"),
])
def test_hgvs_c(summarize, annovar, expected):
    assert summarize.hgvs_c(annovar) == expected


@pytest.mark.parametrize("annovar,expected", [
    ("p.R197H", "p.(Arg197His)"), ("p.R12X", "p.(Arg12Ter)"), ("p.G1355Lfs*167", "p.(Gly1355LeufsTer167)"),
    ("p.L120L", "p.(Leu120=)"), ("p.M1?", "p.(Met1?)"), ("p.120_121del", "p.(120_121del)"),
    ("p.L120delinsLR", "p.(Leu120delinsLeuArg)"),
])
def test_hgvs_p(summarize, annovar, expected):
    assert summarize.hgvs_p(annovar) == expected


def test_hgvs_from_aachange_or_genedetail(summarize):
    assert summarize.hgvs({"aachange": "LPL:NM_000237:exon5:c.G590A:p.R197H"}) == \
        {"transcript": "NM_000237", "hgvs_c": "c.590G>A", "hgvs_p": "p.(Arg197His)"}
    detail = "NM_003122:exon4:c.194+2T>C;NM_001354966:exon4:c.194+2T>C"     # GeneDetail: ';' between transcripts
    assert summarize.hgvs({"aachange": ".", "genedetail": detail}) == \
        {"transcript": "NM_003122", "hgvs_c": "c.194+2T>C", "hgvs_p": "."}
    assert summarize.hgvs({"aachange": "."})["hgvs_c"] == "."


def family_stats(tmp_path, samples, rare, plp=(), sample_counts=None, columns=None):
    stats = tmp_path / "F.json"
    stats.write_text(json.dumps({"sample": "F", "chunk": "all", "genome": "hg38", "counts": {"total": 2},
                                 "func": {}, "columns": columns or {}, "samples": samples, "qc": {"min_dp": 10},
                                 "plp": list(plp), "rare": rare, "sample_counts": sample_counts or {}}))
    return stats


def call(sample, zyg, vaf=0.5, dp=40, conf="high", flags=()):
    return {"sample": sample, "gt": "0/1", "zyg": zyg, "vaf": vaf, "dp": dp, "gq": 50, "alt_reads": 20,
            "flags": list(flags), "conf": conf}


def test_family_report_layout(summarize, tmp_path):
    """Family: one genotype column per member, one row per variant, doubtful rows hidden behind a switch."""
    base = {"func": "exonic", "exonic": "nonsynonymous SNV", "gnomad_af": ".", "revel": "0.9", "filter": "PASS",
            "aachange": "G1:NM_1:exon2:c.G5A:p.R2H"}
    rare = [{**base, "chr": "chr2", "start": "10", "vcf": "chr2-10-C-T", "gene": "G1",
             "genotypes": [call("DAD", "het"), call("K1", "hom", 1.0), call("K2", "hom", 1.0)]},
            {**base, "chr": "chr2", "start": "20", "vcf": "chr2-20-C-T", "gene": "G2",       # father 0/0
             "genotypes": [call("DAD", "ref", 0.0), call("K1", "hom", 1.0), call("K2", "het")]}]
    plp = [{**base, "source": "ClinVar", "chr": "chrM", "start": "3578", "vcf": "chrM-3578-T-G", "gene": "RNR2",
            "clnsig": "Likely_pathogenic", "genotypes": [call("DAD", "het", 0.04, 1026), call("K1", "ref"),
                                                         call("K2", "ref")]}]
    carried = {"DAD": {"coding_rare": 1, "clinvar_plp": 1}, "K1": {"coding_rare": 2}, "K2": {"coding_rare": 2}}
    stats = family_stats(tmp_path, ["DAD", "K1", "K2"], rare, plp, carried,
                         {"gnomad_af": "gnomad41_genome_AF", "clnsig": "CLNSIG"})
    qc = tmp_path / "F.sample_qc.json"
    qc.write_text(json.dumps({"id": "F", "pedigree_source": "inferred", "samples": [
        {"sample": "DAD", "sex_inferred": "male", "ped": {"father": None, "mother": None, "source": "inferred"}},
        {"sample": "K1", "sex_inferred": "female", "ped": {"father": "DAD", "mother": None, "source": "inferred"}},
        {"sample": "K2", "sex_inferred": "female", "ped": {"father": "DAD", "mother": None, "source": "inferred"}}],
        "pairs": [], "trios": []}))
    assert summarize.main(["--stats", str(stats), "--sample-qc", str(qc), "--out", str(tmp_path / "s")]) == 0
    page = (tmp_path / "F.s.html").read_text(encoding="utf-8")
    # the same fixed column widths in every variant table, one column per family member
    cols = re.findall(r"<colgroup>(.*?)</colgroup>", page)
    assert len(cols) == 2 and cols[0] == cols[1]                       # P/LP and candidates: same widths
    assert cols[0].count("<col style='width:112px'>") == 3 + 2         # 3 members, ClinVar, InterVar
    assert "♂ father" in page and "♀ child" in page
    # one row per variant: both children in one row, the Mendelian inconsistency hidden by default
    assert "Inheritance candidates (1)" in page and "doubtful (hidden): 1" in page
    assert "class='conf-high doubtful'" in page and "show 1 doubtful" in page
    assert "homozygous (AR)</span></div><div class='muted small' title='K1, K2 · father het'>" in page
    # HGVS and Sequence Ontology words; heteroplasmy on chrM
    assert "c.5G&gt;A p.(Arg2His)" in page and "missense · NM_1 · chr2-10-C-T" in page
    assert "heteroplasmy" in page and "4% · 1026×" in page
    cand = (tmp_path / "F.s_candidates.tsv").read_text(encoding="utf-8").splitlines()
    head = cand[0].split("\t")
    assert len(cand) == 4 and {"doubtful", "hgvs_c", "hgvs_p", "transcript"} <= set(head)
    rows = [dict(zip(head, line.split("\t"))) for line in cand[1:]]
    assert sorted((r["sample"], r["vcf"], r["doubtful"][:9]) for r in rows) == [
        ("K1", "chr2-10-C-T", "."), ("K1", "chr2-20-C-T", "Mendelian"), ("K2", "chr2-10-C-T", ".")]
    # what every sample carries: a table under the sample QC, the TSV and the JSON
    assert "<h3>Findings per sample</h3>" in page
    assert "<td class='num flag'>1</td>" in page                         # DAD: ClinVar P/LP
    assert "<td class='num' title='homozygous (AR) 1'>1</td>" in page     # K1: the doubtful one not counted
    assert "<th class='num'>InterVar-only P/LP</th>" in page
    samples = (tmp_path / "F.s_samples.tsv").read_text(encoding="utf-8").splitlines()
    srows = {r["sample"]: r for r in (dict(zip(samples[0].split("\t"), line.split("\t"))) for line in samples[1:])}
    assert [(srows[s]["coding_rare"], srows[s]["clinvar_plp"], srows[s]["intervar_plp"],
             srows[s]["inheritance_candidates"]) for s in ("DAD", "K1", "K2")] == \
        [("1", "1", "NA", "0"), ("2", "0", "NA", "1"), ("2", "0", "NA", "1")]
    report = json.loads((tmp_path / "F.s.json").read_text(encoding="utf-8"))
    assert report["sample_findings"][1] == {"sample": "K1", "coding_rare": 2, "clinvar_plp": 0, "intervar_plp": "NA",
                                            "vus_candidates": 0, "inheritance_candidates": 1,
                                            "models": {"homozygous (AR)": 1}}


def test_single_sample_layout(summarize, tmp_path):
    rare = [{"func": "exonic", "exonic": "stopgain", "gnomad_af": ".", "chr": "chr1", "start": "5", "gene": "G",
             "vcf": "chr1-5-C-T", "aachange": "G:NM_9:exon1:c.C5T:p.R2X", "genotypes": [call("S", "hom", 1.0)]}]
    stats = family_stats(tmp_path, ["S"], rare)
    assert summarize.main(["--stats", str(stats), "--out", str(tmp_path / "s")]) == 0
    page = (tmp_path / "F.s.html").read_text(encoding="utf-8")
    assert "<th>genotype</th>" in page and "<col style='width:168px'>" in page
    assert "stop gained" in page and "p.(Arg2Ter)" in page and 'class="doubtful-toggle"' not in page
    assert "Findings per sample" not in page


def test_summary_single_file_has_no_combined_report(summarize, tmp_path):
    p = tmp_path / "0.json"
    p.write_text(json.dumps({"sample": "S1", "chunk": "all", "counts": {"total": 1}, "func": {}, "columns": {}}))
    assert summarize.main(["--stats", str(p), "--out", str(tmp_path / "annotation_summary")]) == 0
    assert (tmp_path / "S1.annotation_summary.html").exists()
    assert not (tmp_path / "annotation_summary.html").exists()


def test_condition_kept_in_full(fmt):
    # The report shortens long text on screen and shows all of it on hover.
    conditions = "|".join(f"Condition_number_{i}\\x2c_type_{i}" for i in range(12))
    txt = ANNO + ["CLNDN", "AAChange.refGene"] + [f"Otherinfo{i}" for i in range(1, 15)]
    stats = fmt.Stats(fmt.build_header(txt, VCF_COLS, len(txt)))
    f = ["chr1", "100", "100", "A", "G", "exonic", "G", "nonsynonymous SNV", "Pathogenic", "0.9", ".",
         conditions, "G:NM_1:exon2:c.A1G:p.M1V,G:NM_2:exon3:c.A9G:p.M3V"]
    v = stats.record(f + ["."] * (len(txt) - len(f)), {}, None)
    assert v["clndn"] == conditions.replace("\\x2c", ",").replace("_", " ") and len(v["clndn"]) > 80
    assert v["aachange"] == "G:NM_1:exon2:c.A1G:p.M1V"
    assert fmt.shorten("abcdef", 4) == "abc…"


PUBMIND = ["PVID", "pathogenicity_sum", "paper_level_pathogenicity_score", "confidence"]


def test_pubmind_fields(fmt):
    txt = ANNO + PUBMIND + [f"Otherinfo{i}" for i in range(1, 15)]
    stats = fmt.Stats(fmt.build_header(txt, VCF_COLS, len(txt)))
    f = ["chr17", "43106487", "43106487", "A", "C", "exonic", "BRCA1", "nonsynonymous SNV", "Pathogenic", "0.9",
         ".", "PVID120370|PVID120371", "Conflicting", "0.9318379141365718", "2"]
    v = stats.record(f + ["."] * (len(txt) - len(f)), {}, None)
    assert (v["pubmind"], v["pubmind_score"], v["pubmind_confidence"], v["pubmind_pvid"]) == \
        ("Conflicting", "0.93", "2", "PVID120370|PVID120371")
    assert stats.columns["pubmind"] == "pathogenicity_sum"


def test_pubmind_columns_must_be_adjacent(fmt):
    # 'confidence' of another database is not PubMind's
    txt = ANNO + ["pathogenicity_sum", "confidence"] + [f"Otherinfo{i}" for i in range(1, 15)]
    stats = fmt.Stats(fmt.build_header(txt, VCF_COLS, len(txt)))
    assert stats.idx["pubmind"] is not None
    assert stats.idx["pubmind_confidence"] is None and stats.idx["pubmind_pvid"] is None


def test_pubmind_report_column(summarize, tmp_path):
    rare = [{"func": "exonic", "exonic": "stopgain", "gnomad_af": ".", "chr": "chr1", "start": "5", "gene": "G",
             "vcf": "chr1-5-C-T", "aachange": "G:NM_9:exon1:c.C5T:p.R2X", "genotypes": [call("S", "hom", 1.0)],
             "pubmind": "Pathogenic/Likely pathogenic", "pubmind_score": "0.94", "pubmind_confidence": "2",
             "pubmind_pvid": "PVID1|PVID2"}]
    stats = family_stats(tmp_path, ["S"], rare)
    assert summarize.main(["--stats", str(stats), "--out", str(tmp_path / "s")]) == 0
    assert ">PubMind</span></th>" not in (tmp_path / "F.s.html").read_text(encoding="utf-8")   # not in protocol
    data = json.loads(stats.read_text())
    data["columns"] = {"pubmind": "pathogenicity_sum"}
    stats.write_text(json.dumps(data))
    assert summarize.main(["--stats", str(stats), "--out", str(tmp_path / "s")]) == 0
    page = (tmp_path / "F.s.html").read_text(encoding="utf-8")
    assert ">PubMind</span></th>" in page and "<col style='width:104px'>" in page
    assert ">P/LP</span>" in page and ">0.94 · conf 2/3</div>" in page and "PVID1, PVID2" in page
    cand = (tmp_path / "F.s_candidates.tsv").read_text(encoding="utf-8").splitlines()
    row = dict(zip(cand[0].split("\t"), cand[1].split("\t")))
    assert (row["pubmind"], row["pubmind_score"], row["pubmind_pvid"]) == \
        ("Pathogenic/Likely pathogenic", "0.94", "PVID1|PVID2")


def test_report_tooltips_only_when_useful(summarize):
    # the script shows a title only when the cell is cut or the title adds detail
    assert "el.removeAttribute('title')" in summarize.JS and "scrollWidth > el.clientWidth" in summarize.JS
