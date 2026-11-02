"""Tests for bin/sample_qc.py (sex, relatedness, PED checks) and bin/inheritance.py."""

import json
import random

import pytest

PSC_HEADER = ("# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t"
              "[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t"
              "[14]nMissing\n")


def psc(rows):
    return PSC_HEADER + "".join("PSC\t0\t" + "\t".join(str(x) for x in r) + "\n" for r in rows)


def test_read_ped(sqc, tmp_path):
    ped = tmp_path / "f.ped"
    ped.write_text("# family sample father mother sex phenotype\n"
                   "F1 KID DAD MOM 1 2\nF1 DAD 0 0 1 1\nF1 MOM 0 0 2 0\n")
    p = sqc.read_ped(str(ped))
    assert p["KID"] == {"family": "F1", "father": "DAD", "mother": "MOM", "sex": "male", "affected": True}
    assert (p["DAD"]["affected"], p["MOM"]["affected"], p["MOM"]["sex"]) == (False, None, "female")
    ped.write_text("F1 KID DAD\n")
    with pytest.raises(ValueError, match="6 columns"):
        sqc.read_ped(str(ped))


def test_ped_spreadsheet_export_and_id_matching(sqc, tmp_path):
    # CSV export of a spreadsheet: header row, dots in the ids where the VCF has underscores
    ped = tmp_path / "families.csv"
    ped.write_text("﻿Family ID,Sample ID,Father,Mother,Sex(1=M,2=F),Phenotype(1=unaffected,2=affected)\n"
                   "FAM1,F1.11,0,0,1,1\nFAM1,F1.21,F1.11,0,2,2\nFAM1,F1.22,F1.11,0,2,1\n"
                   "FAM2,F2.21,0,0,1,2\n", encoding="utf-8")
    raw = sqc.read_ped(str(ped))
    assert sorted(raw) == ["F1.11", "F1.21", "F1.22", "F2.21"]
    matched = sqc.match_ped(raw, ["F1_11", "F1_21", "F1_22"])
    assert sorted(matched) == ["F1_11", "F1_21", "F1_22"]
    assert (matched["F1_21"]["father"], matched["F1_21"]["mother"], matched["F1_21"]["ped_id"]) == \
        ("F1_11", None, "F1.21")
    tsv = tmp_path / "families.tsv"
    tsv.write_text("Family ID\tSample ID\tFather\tMother\tSex\tPhenotype\nF\tA\t0\t0\t2\t2\n")
    assert sqc.read_ped(str(tsv))["A"]["affected"] is True
    # ambiguous normalized ids are not matched
    assert sqc.match_ped({"S.1": raw["F1.11"]}, ["S_1", "s-1"]) == {}


@pytest.mark.parametrize("x,sex", [
    ({"nHets": 25, "nNonRefHom": 0, "nHapAlt": 467}, "male"),
    ({"nHets": 397, "nNonRefHom": 280, "nHapAlt": 0}, "female"),
    ({"nHets": 100, "nNonRefHom": 250, "nHapAlt": 0}, "unclear"),
    ({"nHets": 3, "nNonRefHom": 5, "nHapAlt": 0}, "unknown"),          # too few chrX SNVs
    (None, "unknown"),
])
def test_infer_sex(sqc, x, sex):
    assert sqc.infer_sex(x)[0] == sex


def test_mendelian(sqc):
    assert sqc.mendelian_ok(1, 0, 1) and sqc.mendelian_ok(2, 1, 1) and sqc.mendelian_ok(0, 1, 1)
    assert not sqc.mendelian_ok(2, 0, 1) and not sqc.mendelian_ok(0, 2, 0) and not sqc.mendelian_ok(1, 0, 0)


def simulated_family(n_sites=3000, seed=1):
    """Genotype lines of a father, a mother, their two children and an unrelated person."""
    rng = random.Random(seed)
    lines = []
    for i in range(n_sites):
        af = rng.uniform(0.05, 0.6)
        dad = [rng.random() < af for _ in range(2)]
        mom = [rng.random() < af for _ in range(2)]
        kid1 = [rng.choice(dad), rng.choice(mom)]
        kid2 = [rng.choice(dad), rng.choice(mom)]
        other = [rng.random() < af for _ in range(2)]

        def gt(a):
            return f"{int(a[0])}/{int(a[1])}:60:30"
        lines.append("\t".join([f"chr{1 + i % 22}"] + [gt(x) for x in (dad, mom, kid1, kid2, other)]) + "\n")
    lines.append("chrX\t1/1:60:30\t0/0:60:30\t0/0:60:30\t0/0:60:30\t0/0:60:30\n")     # not autosomal: ignored
    return lines


def test_relatedness_and_pedigree_checks(sqc, tmp_path):
    names = ["DAD", "MOM", "KID1", "KID2", "OTHER"]
    stats = tmp_path / "all.stats"
    stats.write_text(psc([[n, 100, 50, 80, 90, 40, 10, 30.5, 2, 0, 0, 1] for n in names]))
    chrx = tmp_path / "chrx.stats"
    chrx.write_text(psc([["DAD", 0, 0, 2, 0, 0, 0, 0, 0, 0, 60, 0], ["MOM", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0],
                         ["KID1", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0], ["KID2", 0, 0, 1, 0, 0, 0, 0, 0, 0, 55, 0],
                         ["OTHER", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0]]))
    ped = tmp_path / "f.ped"
    # KID2 declared female but has a male chrX; KID1 is affected
    ped.write_text("F1 KID1 DAD MOM 2 2\nF1 KID2 DAD MOM 2 1\nF1 DAD 0 0 1 1\nF1 MOM 0 0 2 1\nF2 OTHER 0 0 2 1\n")
    result = sqc.sample_qc("FAM", str(stats), str(chrx), simulated_family(), str(ped))

    by_sample = {s["sample"]: s for s in result["samples"]}
    assert by_sample["DAD"]["sex_inferred"] == "male" and by_sample["KID1"]["sex_check"] == "ok"
    assert by_sample["KID2"]["sex_check"] == "mismatch"
    assert by_sample["DAD"]["ts_tv"] == 2.25 and by_sample["DAD"]["het_hom_ratio"] == 1.6
    pairs = {(p["a"], p["b"]): p for p in result["pairs"]}
    assert pairs[("DAD", "KID1")]["relationship"] == "parent-child" and pairs[("DAD", "KID1")]["check"] == "ok"
    assert pairs[("KID1", "KID2")]["relationship"] == "full siblings" and pairs[("KID1", "KID2")]["check"] == "ok"
    assert pairs[("DAD", "MOM")]["relationship"] == "unrelated" and pairs[("DAD", "MOM")]["expected"] == "unrelated"
    assert pairs[("KID1", "OTHER")]["check"] == "ok"
    assert 0.2 < pairs[("DAD", "KID1")]["kinship"] < 0.3 and abs(pairs[("DAD", "OTHER")]["kinship"]) < 0.06
    trio = {t["child"]: t for t in result["trios"]}
    assert trio["KID1"]["errors"] == 0 and trio["KID1"]["sites"] > 1000

    # a swapped sample: the PED says parent-child, the data says unrelated
    ped.write_text("F1 KID1 DAD OTHER 2 2\nF1 DAD 0 0 1 1\nF1 OTHER 0 0 2 1\n")
    swapped = sqc.sample_qc("FAM", str(stats), str(chrx), simulated_family(), str(ped))
    bad = {(p["a"], p["b"]): p for p in swapped["pairs"]}[("KID1", "OTHER")]
    assert (bad["expected"], bad["relationship"], bad["check"]) == ("parent-child", "unrelated", "mismatch")
    assert {t["child"]: t for t in swapped["trios"]}["KID1"]["error_rate"] > 0.02


def test_pedigree_without_ped_file(sqc, tmp_path):
    names = ["DAD", "MOM", "KID1", "KID2", "OTHER"]
    stats = tmp_path / "all.stats"
    stats.write_text(psc([[n, 100, 50, 80, 90, 40, 10, 30.5, 2, 0, 0, 1] for n in names]))
    chrx = tmp_path / "chrx.stats"
    chrx.write_text(psc([["DAD", 0, 0, 2, 0, 0, 0, 0, 0, 0, 60, 0], ["MOM", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0],
                         ["KID1", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0], ["KID2", 0, 0, 1, 0, 0, 0, 0, 0, 0, 55, 0],
                         ["OTHER", 0, 40, 50, 0, 0, 0, 0, 0, 0, 0, 0]]))
    # 1. inferred from the relatedness: both children of DAD (male) and MOM (female)
    result = sqc.sample_qc("FAM", str(stats), str(chrx), simulated_family(), None)
    ped = {s["sample"]: s["ped"] for s in result["samples"]}
    assert result["pedigree_source"] == "inferred"
    assert (ped["KID1"]["father"], ped["KID1"]["mother"], ped["KID2"]["mother"]) == ("DAD", "MOM", "MOM")
    assert ped["DAD"]["source"] == "inferred" and ped["OTHER"] is None and ped["DAD"]["affected"] is None
    assert all(p["check"] is None for p in result["pairs"])      # nothing declared to check
    # 2. the VCF header's ##PEDIGREE (DRAGEN): parents' sex from their role, trio checked
    header = tmp_path / "header.txt"
    header.write_text("##fileformat=VCFv4.2\n##PEDIGREE=<Child=KID1,Mother=MOM,Father=DAD>\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDAD\tMOM\tKID1\tKID2\tOTHER\n")
    result = sqc.sample_qc("FAM", str(stats), str(chrx), simulated_family(), None, str(header))
    ped = {s["sample"]: s["ped"] for s in result["samples"]}
    assert result["pedigree_source"] == "vcf header" and ped["KID2"] is None
    assert (ped["DAD"]["sex"], ped["MOM"]["sex"], ped["KID1"]["father"]) == ("male", "female", "DAD")
    assert result["trios"][0]["child"] == "KID1" and result["trios"][0]["errors"] == 0
    assert {(p["a"], p["b"]): p for p in result["pairs"]}[("DAD", "KID1")]["check"] == "ok"


def test_infer_parent_of_siblings(sqc):
    # a father and two daughters, mother not sequenced
    pairs = [{"a": "F", "b": "D1", "relationship": "parent-child"},
             {"a": "F", "b": "D2", "relationship": "parent-child"},
             {"a": "D1", "b": "D2", "relationship": "full siblings"}]
    ped = sqc.infer_pedigree(["F", "D1", "D2"], pairs, {"F": "male", "D1": "female", "D2": "female"})
    assert (ped["D1"]["father"], ped["D1"]["mother"], ped["D2"]["father"]) == ("F", None, "F")
    assert ped["F"]["father"] is None
    # a lone parent-child pair has no direction
    assert sqc.infer_pedigree(["A", "B"], [{"a": "A", "b": "B", "relationship": "parent-child"}], {}) == {}


def test_sample_qc_main_without_genotypes(sqc, tmp_path):
    stats = tmp_path / "all.stats"
    stats.write_text(psc([["S1", 100, 50, 80, 90, 40, 10, 30.5, 2, 0, 0, 1]]))
    out = tmp_path / "S.sample_qc.json"
    assert sqc.main(["--id", "S", "--stats", str(stats), "--out", str(out)]) == 0
    data = json.loads(out.read_text())
    assert data["samples"][0]["sex_inferred"] == "unknown" and data["pairs"] == []


# ---------------------------------------------------------------------------------- inheritance
def variant(vcf, gene, calls, af=".", exonic="nonsynonymous SNV", revel="0.8"):
    chrom, pos = vcf.split("-")[:2]
    return {"vcf": vcf, "chr": chrom, "start": pos, "gene": gene, "func": "exonic", "exonic": exonic,
            "gnomad_af": af, "revel": revel,
            "genotypes": [{"sample": s, "zyg": z} for s, z in calls.items()]}


TRIO = {"KID": {"family": "F", "father": "DAD", "mother": "MOM", "sex": "male", "affected": True},
        "DAD": {"family": "F", "father": None, "mother": None, "sex": "male", "affected": False},
        "MOM": {"family": "F", "father": None, "mother": None, "sex": "female", "affected": False}}
SAMPLES = ["KID", "DAD", "MOM"]


def models(cands):
    return {(c["variant"]["vcf"], c["model"]) for c in cands}


def test_trio_models(inherit):
    vs = [
        variant("chr1-100-A-G", "SCN1A", {"KID": "het", "DAD": "ref", "MOM": "ref"}),           # de novo
        variant("chr1-200-A-G", "SCN1A", {"KID": "het", "DAD": "ref", "MOM": "ref"}, af="0.005"),  # too common
        variant("chr2-300-C-T", "GENE2", {"KID": "hom", "DAD": "het", "MOM": "het"}),             # AR hom
        variant("chr3-400-C-T", "GENE3", {"KID": "het", "DAD": "het", "MOM": "ref"}),             # comp het (pat)
        variant("chr3-500-G-A", "GENE3", {"KID": "het", "DAD": "ref", "MOM": "het"}),             # comp het (mat)
        variant("chr4-600-G-A", "GENE4", {"KID": "het", "DAD": "het", "MOM": "ref"}),             # cis: both paternal
        variant("chr4-700-G-A", "GENE4", {"KID": "het", "DAD": "het", "MOM": "ref"}),
        variant("chrX-5000000-G-A", "GENEX", {"KID": "hemi", "DAD": "ref", "MOM": "het"}),        # X-linked
        variant("chr5-800-G-A", "GENE5", {"KID": "het", "DAD": "missing", "MOM": "ref"}),         # parent missing
    ]
    cands = inherit.analyse(vs, SAMPLES, TRIO, {}, "hg38")
    assert models(cands) == {
        ("chr1-100-A-G", "de novo"), ("chr2-300-C-T", "homozygous (AR)"),
        ("chr3-400-C-T", "compound heterozygous"), ("chr3-500-G-A", "compound heterozygous"),
        ("chrX-5000000-G-A", "X-linked (hemizygous)")}
    comp = [c for c in cands if c["model"] == "compound heterozygous"]
    assert comp[0]["partners"] and comp[0]["gene"] == "GENE3"
    assert cands[0]["model"] == "de novo"                    # report order


def test_segregation_excludes(inherit):
    ped = dict(TRIO, SIB={"family": "F", "father": "DAD", "mother": "MOM", "sex": "female", "affected": False})
    vs = [variant("chr2-300-C-T", "GENE2", {"KID": "hom", "DAD": "het", "MOM": "het", "SIB": "hom"}),
          variant("chr2-310-C-T", "GENE2B", {"KID": "hom", "DAD": "het", "MOM": "het", "SIB": "het"})]
    assert models(inherit.analyse(vs, SAMPLES + ["SIB"], ped, {}, "hg38")) == {("chr2-310-C-T", "homozygous (AR)")}


def test_singleton_and_affected_parent(inherit):
    # single sample, no PED: recessive models only, compound het unphased
    vs = [variant("chr3-400-C-T", "GENE3", {"S": "het"}), variant("chr3-500-G-A", "GENE3", {"S": "het"}),
          variant("chr1-100-A-G", "GENE1", {"S": "het"}), variant("chr2-1-A-G", "GENE2", {"S": "hom"}),
          variant("chrX-5000000-G-A", "GENEX", {"S": "hom"})]
    cands = inherit.analyse(vs, ["S"], {}, {"S": "male"}, "hg38")
    assert models(cands) == {("chr3-400-C-T", "possible compound heterozygous"),
                             ("chr3-500-G-A", "possible compound heterozygous"),
                             ("chr2-1-A-G", "homozygous (AR)"), ("chrX-5000000-G-A", "X-linked (hemizygous)")}
    # several samples without PED: no proband, no candidates
    assert inherit.analyse(vs, ["S", "T"], {}, {}, "hg38") == []
    # dominant: carried by an affected parent
    ped = {k: dict(v) for k, v in TRIO.items()}
    ped["DAD"]["affected"] = True
    vs = [variant("chr6-1-A-G", "GENE6", {"KID": "het", "DAD": "het", "MOM": "ref"})]
    assert models(inherit.analyse(vs, SAMPLES, ped, {}, "hg38")) == {("chr6-1-A-G", "dominant (affected parent)")}


def test_children_are_probands_without_affected_status(inherit):
    ped = {k: dict(v, affected=None, source="vcf header") for k, v in TRIO.items()}
    assert inherit.probands(SAMPLES, ped) == ["KID"]
    vs = [variant("chr1-100-A-G", "SCN1A", {"KID": "het", "DAD": "ref", "MOM": "ref"})]
    vs[0]["genotypes"][0].update(dn="DeNovo", dq=25.0)
    cands = inherit.analyse(vs, SAMPLES, ped, {}, "hg38")
    assert [(c["model"], c["note"]) for c in cands] == \
        [("de novo", "father and mother reference · DRAGEN DeNovo (DQ 25)")]
    # PED with affected status but nobody affected: no proband
    assert inherit.probands(SAMPLES, {k: dict(v, affected=False) for k, v in TRIO.items()}) == []


def test_de_novo_and_recessive_filters(inherit):
    ped = {k: dict(v, affected=None, source="vcf header") for k, v in TRIO.items()}
    low = variant("chr1-100-A-G", "G1", {"KID": "het", "DAD": "ref", "MOM": "ref"})
    low["genotypes"][1]["conf"] = "low"                                 # father's reference call
    lowdq = variant("chr1-200-A-G", "G2", {"KID": "het", "DAD": "ref", "MOM": "ref"})
    lowdq["genotypes"][0].update(dn="LowDQ", dq=0.0)                    # DRAGEN: not de novo
    alt_in_parent = variant("chr1-300-A-G", "G3", {"KID": "het", "DAD": "ref", "MOM": "ref"})
    alt_in_parent["genotypes"][2]["flags"] = ["ALT_READS_IN_REF"]
    good = variant("chr1-400-A-G", "G4", {"KID": "het", "DAD": "ref", "MOM": "ref"})
    hom_in_parent = variant("chr2-1-C-T", "G5", {"KID": "hom", "DAD": "het", "MOM": "hom"})
    comp = [variant("chr3-1-C-T", "G6", {"KID": "het", "DAD": "het", "MOM": "ref"}),
            variant("chr3-2-C-T", "G6", {"KID": "het", "DAD": "ref", "MOM": "hom"})]   # mother homozygous
    cands = inherit.analyse([low, lowdq, alt_in_parent, good, hom_in_parent, *comp], SAMPLES, ped, {}, "hg38")
    assert models(cands) == {("chr1-400-A-G", "de novo")}


def test_without_parents_homozygous_damaging_only(inherit):
    vs = [variant("chr2-1-A-G", "G1", {"S": "hom"}, revel="0.1"),                        # benign missense
          variant("chr2-2-A-G", "G2", {"S": "hom"}, revel="0.7"),                        # damaging missense
          variant("chr2-3-A-AT", "G3", {"S": "hom"}, exonic="frameshift insertion", revel=".")]
    ncrna = variant("chr2-4-A-G", "G4", {"S": "hom"}, exonic=".")
    ncrna["func"] = "ncRNA_splicing"                                                       # not protein-coding
    cands = inherit.analyse(vs + [ncrna], ["S"], {}, {"S": "female"}, "hg38")
    assert models(cands) == {("chr2-2-A-G", "homozygous (AR)"), ("chr2-3-A-AT", "homozygous (AR)")}
    # a sequenced parent phases the call: any rare homozygous variant is a candidate
    ped = {"KID": {"father": "DAD", "affected": True}, "DAD": {"affected": False}}
    benign = variant("chr2-1-A-G", "G1", {"KID": "hom", "DAD": "het"}, revel="0.1")
    assert models(inherit.analyse([benign], ["KID", "DAD"], ped, {}, "hg38")) == {("chr2-1-A-G", "homozygous (AR)")}


def test_doubtful_candidates(inherit):
    ped = {k: dict(v, affected=None, source="inferred") for k, v in TRIO.items()}
    mendel = variant("chr2-300-C-T", "G1", {"KID": "hom", "DAD": "ref", "MOM": "het"})        # father 0/0
    low = variant("chr2-400-C-T", "G2", {"KID": "hom", "DAD": "het", "MOM": "het"})
    low["genotypes"][0]["conf"] = "low"
    good = variant("chr2-500-C-T", "G3", {"KID": "hom", "DAD": "het", "MOM": "het"})
    cands = inherit.analyse([mendel, low, good], SAMPLES, ped, {}, "hg38")
    reasons = {c["variant"]["vcf"]: c["doubtful"] for c in cands}
    assert reasons["chr2-500-C-T"] == ""
    assert reasons["chr2-300-C-T"].startswith("Mendelian inconsistency: father homozygous reference")
    assert reasons["chr2-400-C-T"] == "low-confidence genotype call"
    assert cands[0]["variant"]["vcf"] == "chr2-500-C-T"                  # doubtful ones last


def test_x_non_par(inherit):
    assert inherit.x_non_par({"chr": "chrX", "start": "5000000"}, "hg38")
    assert not inherit.x_non_par({"chr": "chrX", "start": "100000"}, "hg38")      # PAR1
    assert not inherit.x_non_par({"chr": "7", "start": "5000000"}, "hg19")
