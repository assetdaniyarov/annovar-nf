"""Tests for bin/gene_info.py, `annovar_db.py genes` and the gene information in the report."""

import json
import re

import pytest
from conftest import load_script

gi = load_script("gene_info")

GENCC_HEAD = ["uuid", "gene_curie", "gene_symbol", "disease_curie", "disease_title", "classification_title",
              "moi_title", "submitter_title", "submitted_as_notes"]
GENCC = [
    ["1", "HGNC:10585", "SCN1A", "MONDO:1", "Dravet syndrome", "Definitive", "Autosomal dominant", "ClinGen", ""],
    ["2", "HGNC:10585", "SCN1A", "MONDO:1", "Dravet syndrome", "Strong", "Autosomal dominant", "Ambry Genetics",
     '"a note over\ntwo lines with ""quotes"""'],
    ["3", "HGNC:10585", "SCN1A", "MONDO:2", "familial febrile seizures", "Limited", "Autosomal recessive",
     "Ambry Genetics", ""],
    ["4", "HGNC:20", "AARS1", "MONDO:3", "leukoencephalopathy; early onset", "Moderate", "Autosomal recessive",
     "Genomics England PanelApp", ""],
    ["5", "HGNC:30", "NOPE", "MONDO:4", "something", "Refuted Evidence", "Autosomal dominant", "ClinGen", ""],
]
CONSTRAINT_HEAD = ["gene", "gene_id", "transcript", "canonical", "mane_select", "lof.pLI", "lof.oe_ci.upper",
                   "mis.z_score"]
CONSTRAINT = [
    ["SCN1A", "6323", "NM_001165963.4", "true", "true", "1.0000e+00", "9.1000e-02", "5.71e+00"],
    ["SCN1A", "ENSG1", "ENST1", "true", "true", "9.0000e-01", "2.0000e-01", "5.0e+00"],
    ["SCN1A", "ENSG1", "ENST2", "false", "false", "0", "1.5", "0"],
    ["AARS", "ENSG2", "ENST3", "true", "false", "2.0e-03", "6.5432e-01", "NA"],      # previous symbol
    ["LONE", "ENSG3", "ENST4", "true", "true", "NA", "1.2e+00", "-1.0"],
]
HGNC_HEAD = ["hgnc_id", "symbol", "name", "status", "prev_symbol"]
HGNC = [
    ["HGNC:10585", "SCN1A", "sodium voltage-gated channel alpha subunit 1", "Approved", ""],
    ["HGNC:20", "AARS1", "alanyl-tRNA synthetase 1", "Approved", '"AARS|AARS0"'],
    ["HGNC:40", "LONE", "lonely gene", "Approved", ""],
    ["HGNC:50", "OLD", "withdrawn gene", "Entry Withdrawn", ""],
]


def tsv(path, head, rows):
    path.write_text("\n".join("\t".join(r) for r in [head, *rows]) + "\n", encoding="utf-8", newline="\n")
    return path


@pytest.fixture
def sources(tmp_path):
    return {"gencc": tsv(tmp_path / "gencc.tsv", GENCC_HEAD, GENCC),
            "constraint": tsv(tmp_path / "constraint.tsv", CONSTRAINT_HEAD, CONSTRAINT),
            "hgnc": tsv(tmp_path / "hgnc.tsv", HGNC_HEAD, HGNC)}


def read(path, parser):
    with open(path, encoding="utf-8", newline="") as fh:
        return parser(fh)


def build(sources):
    return gi.build_rows(read(sources["gencc"], gi.read_gencc), read(sources["constraint"], gi.read_constraint),
                         read(sources["hgnc"], gi.read_hgnc))


def test_read_sources(sources):
    gencc = read(sources["gencc"], gi.read_gencc)
    assert set(gencc) == {"SCN1A", "AARS1"}                  # refuted evidence left out
    dravet = gencc["SCN1A"]["diseases"]["MONDO:1"]
    assert (dravet["moi"], dravet["rank"], len(dravet["submitters"])) == ({"AD"}, 0, 2)
    assert gencc["AARS1"]["diseases"]["MONDO:3"]["title"] == "leukoencephalopathy, early onset"
    constraint = read(sources["constraint"], gi.read_constraint)
    assert constraint["SCN1A"] == ("0.09", "1.00", "5.71")        # MANE Select, RefSeq row
    assert constraint["AARS"] == ("0.65", "0.00", ".") and constraint["LONE"] == ("1.20", ".", "-1.00")
    hgnc = read(sources["hgnc"], gi.read_hgnc)
    assert hgnc["AARS1"]["previous"] == ["AARS", "AARS0"] and "OLD" not in hgnc


def test_build_and_load(sources, tmp_path):
    rows = build(sources)
    by_gene = {r[0]: r for r in rows}
    assert set(by_gene) == {"SCN1A", "AARS1", "LONE"}             # AARS listed under its approved symbol
    assert by_gene["SCN1A"][3] == "AD"                            # the 'Limited' AR disease does not count
    assert by_gene["SCN1A"][4].split(";")[0] == "Dravet syndrome|AD|Definitive|2"
    assert by_gene["AARS1"][5] == "0.65" and by_gene["AARS1"][8] == "AARS,AARS0"
    path = tmp_path / "genes.tsv"
    gi.write_table(rows, path, ["test table", "GenCC (CC0)"])
    table = gi.load(str(path))
    assert table.sources == ["test table", "GenCC (CC0)"]
    scn1a = table.lookup("NONE;scn1a")
    assert scn1a.gene == "SCN1A" and scn1a.summary() == "AD · LOEUF 0.09"
    assert scn1a.disease_text().startswith("Dravet syndrome (AD, Definitive); familial febrile seizures (AR, Limited)")
    assert "intolerant of loss of function" in scn1a.detail() and "2 GenCC submitters" in scn1a.detail()
    assert table.get("AARS").gene == "AARS1" and table.symbol("aars0") == "AARS1"
    assert table.get("LONE").moi == [] and table.get("XYZ") is None and gi.load(None).genes == {}


@pytest.mark.parametrize("model,moi,fits", [
    ("de novo", ["AD"], True), ("de novo", ["AR"], False), ("homozygous (AR)", ["AR", "AD"], True),
    ("compound heterozygous", ["AD"], False), ("X-linked (hemizygous)", ["XLR"], True), ("de novo", [], None),
])
def test_moi_fit(model, moi, fits):
    gene = gi.Gene("G", ".", "", moi, [], None, None, None)
    assert gi.moi_fit(model, gene) is fits
    assert gi.moi_fit(model, None) is None


def test_read_panel(tmp_path):
    table = gi.GeneTable(genes={"AARS1": None, "SCN1A": None}, aliases={"AARS": "AARS1"})
    plain = tmp_path / "epilepsy.txt"
    plain.write_text("# my panel\ngene\nSCN1A\naars, alanyl-tRNA synthetase\n")
    panel = gi.read_panel(str(plain), table)
    assert (panel.name, panel.genes) == ("epilepsy", {"SCN1A", "AARS1"})     # header word and comments skipped
    panelapp = tmp_path / "Epilepsy.tsv"
    panelapp.write_text("Entity Name\tEntity type\tGene Symbol\tLevel4\tGEL_Status\tversion\r\n"
                        "AARS\tgene\tAARS\tEarly onset epilepsy\t3\t9.94\r\n"
                        "KCNQ2\tgene\tKCNQ2\tEarly onset epilepsy\t4\t9.94\r\n"
                        "ABC\tgene\tABC\tEarly onset epilepsy\t2\t9.94\r\n"
                        "region_1\tregion\t\tEarly onset epilepsy\t3\t9.94\r\n", encoding="utf-8")
    panel = gi.read_panel(str(panelapp), table)
    assert panel.name == "Early onset epilepsy (PanelApp v9.94)" and panel.genes == {"AARS1", "KCNQ2"}
    assert gi.panels_of("AARS;X", [panel], table) == ["Early onset epilepsy (PanelApp v9.94)"]
    assert gi.panels_of("ABC", [panel], table) == [] and gi.panels_of(".", [panel]) == []


def test_annovar_db_genes(sources, tmp_path):
    db = load_script("annovar_db")
    humandb = tmp_path / "humandb"
    args = ["genes", "--humandb", str(humandb), "--gencc", str(sources["gencc"]),
            "--constraint", str(sources["constraint"]), "--hgnc", str(sources["hgnc"])]
    assert db.main(args) == 0
    table = humandb / gi.TABLE_NAME
    text = table.read_text(encoding="utf-8")
    assert text.startswith("## annovar-nf gene table, built ") and "\nSCN1A\tHGNC:10585\t" in text
    assert table.parent.name == "annovar_nf" and not list(table.parent.glob(".annovar-nf-staging.*"))
    table.write_text("## old\n")
    assert db.main(args) == 0 and table.read_text() == "## old\n"          # kept: no --max-age
    assert db.main([*args, "--max-age", "0"]) == 0 and "SCN1A" in table.read_text(encoding="utf-8")
    missing = ["genes", "--humandb", str(humandb), "--max-age", "0", "--gencc", str(tmp_path / "none.tsv"),
               "--constraint", str(sources["constraint"]), "--hgnc", str(sources["hgnc"])]
    assert db.main([*missing, "--optional"]) == 0                 # a warning in the pipeline
    with pytest.raises(OSError):
        db.main(missing)


def test_report_gene_information_and_panel(summarize, tmp_path, sources):
    """Modes of inheritance, LOEUF and diseases next to the gene; link menu; panel switch; model fit."""
    def call(sample, zyg):
        return {"sample": sample, "gt": "0/1", "zyg": zyg, "vaf": 0.5, "dp": 40, "gq": 50, "alt_reads": 20,
                "flags": [], "conf": "high"}

    base = {"func": "exonic", "exonic": "stopgain", "gnomad_af": ".", "filter": "PASS", "clndn": "."}
    rare = [{**base, "chr": "chr2", "start": "10", "vcf": "chr2-10-C-T", "gene": "SCN1A", "clnalleleid": "77",
             "genotypes": [call("P", "hom")]},
            {**base, "chr": "chr3", "start": "20", "vcf": "chr3-20-G-A", "gene": "AARS",
             "genotypes": [call("P", "hom")]}]
    stats = tmp_path / "F.json"
    stats.write_text(json.dumps({"sample": "F", "chunk": "all", "genome": "hg38", "counts": {"total": 2}, "func": {},
                                 "columns": {}, "samples": ["P"], "qc": {}, "plp": [], "rare": rare}))
    genes = tmp_path / "genes.tsv"
    rows = build(sources)
    gi.write_table(rows, genes, ["annovar-nf gene table, built 2026-11-01", "GenCC gene-disease validity (CC0)"])
    panel = tmp_path / "metabolic.txt"
    panel.write_text("AARS1\n")
    assert summarize.main(["--stats", str(stats), "--out", str(tmp_path / "s"), "--genes", str(genes),
                           "--panel", str(panel)]) == 0
    page = (tmp_path / "F.s.html").read_text(encoding="utf-8")
    # gene: modes of inheritance and LOEUF (constrained in red); diseases under the ClinVar condition
    assert "<span class='moi'>AD</span><span class='loeuf lof'>LOEUF 0.09</span>" in page
    assert "<span class='moi'>AR</span><span class='loeuf'>LOEUF 0.65</span>" in page
    assert "GenCC: Dravet syndrome (AD) +1</div>" in page
    # link menu of the gene name: variant and gene pages, IGV
    links = json.loads(re.search(r"data-links='([^']*)' title='SCN1A: links'", page)[1].replace("&quot;", '"'))
    assert [label for label, _ in links] == ["ClinVar", "gnomAD", "Franklin", "VarSome", "OMIM", "GenCC", "IGV"]
    assert dict(links)["Franklin"].endswith("/variant/snp/chr2-10-C-T-hg38")
    assert dict(links)["VarSome"] == "https://varsome.com/variant/hg38/chr2-10-C-T"
    assert dict(links)["IGV"] == "http://localhost:60151/goto?locus=chr2:10"
    # homozygous: fits the AR gene, not the AD one; fitting candidates first
    assert "title='gene: AR (GenCC) - fits the model'>AR</span>" in page
    assert "class='fit no' title='gene: AD (GenCC) - does not fit the model'>AD</span>" in page
    assert page.index("fits the model'>AR") < page.index("does not fit the model'>AD")
    # panel: SCN1A outside, hidden behind a switch; counts in the headings
    assert "Inheritance candidates (2 · 1 in the gene panel)" in page and "show 1 outside the gene panel" in page
    assert "<b>metabolic</b> (1 genes)" in page and "GenCC gene-disease validity (CC0)" in page
    cand = (tmp_path / "F.s_candidates.tsv").read_text(encoding="utf-8").splitlines()
    head = cand[0].split("\t")
    got = {r["gene"]: r for r in (dict(zip(head, line.split("\t"))) for line in cand[1:])}
    assert (got["AARS"]["moi_fit"], got["AARS"]["gene_moi"], got["AARS"]["panels"], got["AARS"]["loeuf"]) == \
        ("yes", "AR", "metabolic", "0.65")
    assert (got["SCN1A"]["moi_fit"], got["SCN1A"]["panels"]) == ("no", ".")
    report = json.loads((tmp_path / "F.s.json").read_text(encoding="utf-8"))
    assert report["gene_panels"] == [{"name": "metabolic", "genes": 1}]
    assert report["gene_information"][1] == "GenCC gene-disease validity (CC0)"
    assert {c["variant"]["gene"]: c["moi_fit"] for c in report["inheritance_candidates"]} == \
        {"AARS": True, "SCN1A": False}
