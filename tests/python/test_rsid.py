"""dbSNP rsIDs of the listed variants (format_annovar_table.py) and their links in the report."""

from test_format_and_summary import ANNO, VCF_COLS, row

HEADER = ANNO + ["avsnp151", "Otherinfo1", "Otherinfo2", "Otherinfo3"] + VCF_COLS
ID_COL = HEADER.index("ID")


def fields(rs=".", vcf_id="."):
    f = row("chr1", 100, "A", "G", "exonic", "nonsynonymous SNV", "Pathogenic", "0.9").split("\t")
    f = f[:len(ANNO)] + [rs] + f[len(ANNO):]
    f[ID_COL] = vcf_id
    return f


def test_rsid_from_dbsnp_then_vcf_id(fmt):
    s = fmt.Stats(HEADER)
    assert s.rsid(fields("rs123", "rs999")) == "rs123"          # the protocol's dbSNP first
    assert s.rsid(fields(".", "rs9;rs10")) == "rs9;rs10"        # else the VCF ID column
    assert s.rsid(fields(".", "chr1_100_A_G")) == "."           # caller IDs are not rsIDs
    assert s.rsid(fields()) == "."


def test_rsid_without_dbsnp_column(fmt):
    header = ANNO + ["Otherinfo1", "Otherinfo2", "Otherinfo3"] + VCF_COLS
    f = row("chr1", 100, "A", "G", "exonic", "nonsynonymous SNV", "Pathogenic", "0.9").split("\t")
    f[header.index("ID")] = "rs42"
    assert fmt.Stats(header).rsid(f) == "rs42"


def test_listed_variant_carries_rsid(fmt):
    s = fmt.Stats(HEADER)
    s.add("\t".join(fields("rs123")) + "\n")
    assert [v["rsid"] for v in s.plp] == ["rs123"]


def test_report_links_rsid(summarize):
    v = {"rsid": "rs1;rs2", "vcf": "chr1-100-A-G", "gene": "G"}
    assert ["dbSNP", "https://www.ncbi.nlm.nih.gov/snp/rs1"] in summarize.variant_links(v, "hg38")
    html = summarize.rsid_html(v)
    assert "https://www.ncbi.nlm.nih.gov/snp/rs1" in html and ">rs2</a>" in html
    # stats written before the rsid field: no dbSNP entry, nothing shown
    old = {"vcf": "chr1-100-A-G", "gene": "G"}
    assert all(label != "dbSNP" for label, _ in summarize.variant_links(old, "hg38"))
    assert summarize.rsid_html(old) == "" and "rsid" in summarize.VARIANT_COLUMNS
