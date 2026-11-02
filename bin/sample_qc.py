#!/usr/bin/env python3
"""
sample_qc.py
------------
Per-sample quality control of one input VCF, for the summary report:

  * counts from `bcftools stats -s -` (PSC lines): heterozygous / homozygous
    SNVs, Ti/Tv, indels, mean depth, singletons, missing genotypes;
  * sex from chromosome X outside the pseudo-autosomal regions (a second
    `bcftools stats -s -` restricted to it): the ratio of heterozygous to
    homozygous + hemizygous SNVs is about 0 in males and around 1 or above in
    females;
  * relatedness of every sample pair (files with few samples): KING-robust
    kinship and the opposite-homozygote (IBS0) rate on high-quality autosomal
    SNV genotypes, read from
    `bcftools query -f '%CHROM[\\t%GT:%GQ:%DP]\\n'` (bi-allelic PASS SNVs);
  * the pedigree: a PED file (or a CSV / TSV export of a spreadsheet with
    the PED columns; ids matched to the VCF ignoring case and '.' / '_' /
    '-'), else the VCF header's ##PEDIGREE lines (DRAGEN joint calling), else
    inferred from the relatedness (trios, a parent with siblings). Declared
    sex and relationships are checked against the data, and every trio gets
    its Mendelian error rate.

Usage:
    sample_qc.py --id FAM1 --stats all.stats --chrx chrx.stats \\
        [--genotypes gt.tsv|-] [--ped family.ped] --out FAM1.sample_qc.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from pathlib import Path

# Sex from chrX (non-PAR) heterozygous / (homozygous + hemizygous) SNVs.
MALE_MAX_X_RATIO = 0.25
FEMALE_MIN_X_RATIO = 0.6
MIN_X_SNVS = 30
# Genotypes used for relatedness and Mendelian errors.
REL_MIN_GQ = 20
REL_MIN_DP = 10
REL_MIN_SITES = 1000
# KING kinship boundaries (Manichaikul et al. 2010) and the IBS0 rate that
# separates parent-child (about 0) from full siblings.
KINSHIP_DEGREES = ((0.354, "duplicate / identical twin"), (0.177, "first-degree"),
                   (0.0884, "second-degree"), (0.0442, "third-degree"))
PARENT_CHILD_MAX_IBS0 = 0.005
AUTOSOMES = {str(i) for i in range(1, 23)}


PED_CODES = {"-9", "0", "1", "2"}


def split_ped_line(line: str) -> list[str]:
    """Tab-, comma- or whitespace-separated columns (PED files and spreadsheet exports)."""
    line = line.rstrip("\r\n")
    for sep in ("\t", ","):
        if sep in line:
            return [x.strip().strip('"') for x in line.split(sep)]
    return line.split()


def read_ped(path: str | None) -> dict[str, dict]:
    """
    PED file (family, sample, father, mother, sex 1/2, phenotype 2 = affected)
    keyed by sample. Also a CSV / TSV export of a spreadsheet with these six
    columns first: a header row (sex / phenotype not a PED code) is skipped.
    """
    ped: dict[str, dict] = {}
    if not path:
        return ped
    first = True
    with open(path, encoding="utf-8-sig") as fh:
        for n, line in enumerate(fh, 1):
            f = split_ped_line(line)
            if not any(f) or f[0].startswith("#"):
                continue
            if len(f) < 6:
                raise ValueError(f"{path}:{n}: a PED line needs 6 columns (family sample father mother sex phenotype)")
            if first and (f[4] not in PED_CODES or f[5] not in PED_CODES):
                first = False
                continue        # header row
            first = False
            ped[f[1]] = {
                "family": f[0], "father": None if f[2] in ("", "0", ".") else f[2],
                "mother": None if f[3] in ("", "0", ".") else f[3],
                "sex": {"1": "male", "2": "female"}.get(f[4], "unknown"),
                "affected": {"2": True, "1": False}.get(f[5]),
            }
    return ped


def read_vcf_pedigree(path: str | None) -> dict[str, dict]:
    """
    '##PEDIGREE=<Child=C,Mother=M,Father=F>' (or ID=) lines of a VCF header, as
    written by DRAGEN joint calling with a pedigree; the parents' sex follows
    from their role, affected status is unknown.
    """
    out: dict[str, dict] = {}
    if not path or not Path(path).exists():
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                break
            if not line.startswith("##PEDIGREE=<"):
                continue
            fields = dict(kv.split("=", 1) for kv in line.strip()[len("##PEDIGREE=<"):-1].split(",") if "=" in kv)
            child = fields.get("Child") or fields.get("ID")
            if not child or not (fields.get("Father") or fields.get("Mother")):
                continue
            out[child] = {"family": "vcf", "father": fields.get("Father"), "mother": fields.get("Mother"),
                          "sex": out.get(child, {}).get("sex", "unknown"), "affected": None}
            for parent, sex in ((fields.get("Father"), "male"), (fields.get("Mother"), "female")):
                if parent:
                    entry = out.setdefault(parent, {"family": "vcf", "father": None, "mother": None, "affected": None})
                    entry["sex"] = sex
    return out


def infer_pedigree(names: list[str], pairs: list[dict], sex: dict[str, str]) -> dict[str, dict]:
    """
    Parents and children from the data alone: a sample with parent-child
    relationships to two samples unrelated to each other is their child; a
    sample in a parent-child relationship with two full siblings is their
    parent. Father / mother by inferred sex. A lone parent-child pair (no
    direction) gives nothing.
    """
    rel = {}
    for p in pairs:
        rel[(p["a"], p["b"])] = rel[(p["b"], p["a"])] = p["relationship"]
    pc = {s: [t for t in names if rel.get((s, t)) == "parent-child"] for s in names}
    parents_of: dict[str, tuple] = {}
    for child in names:
        for a, b in itertools.combinations(pc[child], 2):
            if rel.get((a, b)) in ("unrelated", "third-degree"):
                parents_of[child] = (a, b)
                break
    for parent in names:
        if parent in parents_of:
            continue
        for x, y in itertools.combinations(pc[parent], 2):
            if rel.get((x, y)) == "full siblings":
                for kid in (x, y):
                    parents_of.setdefault(kid, (parent,))
    out: dict[str, dict] = {}
    for kid, parents in parents_of.items():
        father = next((p for p in parents if sex.get(p) == "male"), None)
        mother = next((p for p in parents if sex.get(p) == "female"), None)
        for p in parents:           # sex unknown / unclear: fill the free role
            if p not in (father, mother):
                if father is None:
                    father = p
                elif mother is None:
                    mother = p
        out[kid] = {"family": "inferred", "father": father, "mother": mother, "sex": "unknown", "affected": None}
    for entry in list(out.values()):
        for p in (entry["father"], entry["mother"]):
            if p and p not in out:
                out[p] = {"family": "inferred", "father": None, "mother": None, "sex": "unknown", "affected": None}
    return out


def id_key(sample: str) -> str:
    """Sample id for matching PED and VCF names: case and '.' / '_' / '-' / space separators ignored."""
    return re.sub(r"[._\- ]+", "_", sample.strip()).lower()


def match_ped(ped: dict[str, dict], names: list[str]) -> dict[str, dict]:
    """
    PED entries re-keyed by the VCF sample names (F1.01 in the PED matches
    F1_01 in the VCF); parents are renamed the same way. 'ped_id' keeps the
    PED's own id. Ambiguous matches are left out.
    """
    by_key: dict[str, list[str]] = {}
    for name in names:
        by_key.setdefault(id_key(name), []).append(name)

    def vcf_name(pid: str | None) -> str | None:
        if pid is None or pid in names:
            return pid
        hits = by_key.get(id_key(pid), [])
        return hits[0] if len(hits) == 1 else pid

    out = {}
    for pid, entry in ped.items():
        name = vcf_name(pid)
        if name in names and name not in out:
            out[name] = {**entry, "father": vcf_name(entry["father"]), "mother": vcf_name(entry["mother"]),
                         "ped_id": pid}
    return out


def read_psc(path: str | None) -> dict[str, dict]:
    """Per-sample counts (PSC lines) of `bcftools stats -s -`, keyed by sample."""
    out: dict[str, dict] = {}
    if not path or not Path(path).exists():
        return out
    names: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# PSC"):
                names = [re.sub(r"^\[\d+\]", "", c) for c in line.rstrip("\n").split("\t")]
            elif line.startswith("PSC\t"):
                f = line.rstrip("\n").split("\t")
                row = dict(zip(names, f))
                skip = ("# PSC", "id", "sample")
                out[row.get("sample", f[2])] = {k: _num(v) for k, v in row.items() if k not in skip}
    return out


def _num(value: str):
    try:
        x = float(value)
    except ValueError:
        return value
    return int(x) if x.is_integer() else x


def infer_sex(x: dict | None) -> tuple[str, float | None, int]:
    """('male' | 'female' | 'unclear' | 'unknown', het/hom ratio, informative SNVs) from chrX PSC counts."""
    if not x:
        return "unknown", None, 0
    het = x.get("nHets", 0)
    hom = x.get("nNonRefHom", 0) + x.get("nHapAlt", 0)
    n = het + hom
    ratio = het / hom if hom else (float("inf") if het else None)
    if n < MIN_X_SNVS or ratio is None:
        return "unknown", ratio, n
    if ratio <= MALE_MAX_X_RATIO:
        return "male", ratio, n
    if ratio >= FEMALE_MIN_X_RATIO:
        return "female", ratio, n
    return "unclear", ratio, n


def alt_count(gt: str) -> int | None:
    """Alt alleles of a diploid bi-allelic genotype; None when missing or not diploid."""
    a = re.split(r"[/|]", gt)
    if len(a) != 2 or "." in a:
        return None
    return sum(1 for x in a if x != "0")


def good_call(field: str) -> int | None:
    """'GT:GQ:DP' -> alt allele count if the call passes the GQ / DP limits (unknown values pass)."""
    parts = (field.split(":") + [".", "."])[:3]
    c = alt_count(parts[0])
    if c is None:
        return None
    for value, limit in ((parts[1], REL_MIN_GQ), (parts[2], REL_MIN_DP)):
        try:
            if float(value) < limit:
                return None
        except ValueError:
            pass
    return c


TRANSMIT = {0: (0,), 1: (0, 1), 2: (1,)}


def mendelian_ok(child: int, father: int, mother: int) -> bool:
    return any(a + b == child for a in TRANSMIT[father] for b in TRANSMIT[mother])


def relatedness(lines, names: list[str], trios: list[tuple[str, str, str]]) -> tuple[list[dict], list[dict]]:
    """Pairwise kinship / IBS0 and trio Mendelian errors from '%CHROM[\\t%GT:%GQ:%DP]' lines."""
    pairs = list(itertools.combinations(range(len(names)), 2))
    pc = [{"n": 0, "het_a": 0, "het_b": 0, "het_het": 0, "ibs0": 0} for _ in pairs]
    pos = {s: i for i, s in enumerate(names)}
    tri = [(pos[c], pos[f], pos[m]) for c, f, m in trios]
    tc = [{"n": 0, "errors": 0} for _ in tri]
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) != len(names) + 1 or f[0].removeprefix("chr") not in AUTOSOMES:
            continue
        g = [good_call(x) for x in f[1:]]
        for k, (i, j) in enumerate(pairs):
            a, b = g[i], g[j]
            if a is None or b is None or (a == 0 and b == 0):
                continue
            c = pc[k]
            c["n"] += 1
            c["het_a"] += a == 1
            c["het_b"] += b == 1
            c["het_het"] += a == 1 and b == 1
            c["ibs0"] += abs(a - b) == 2
        for k, (ci, fi, mi) in enumerate(tri):
            child, father, mother = g[ci], g[fi], g[mi]
            if None in (child, father, mother) or child + father + mother == 0:
                continue
            tc[k]["n"] += 1
            tc[k]["errors"] += not mendelian_ok(child, father, mother)
    out_pairs = []
    for (i, j), c in zip(pairs, pc):
        het = c["het_a"] + c["het_b"]
        kinship = (c["het_het"] - 2 * c["ibs0"]) / het if het else None
        ibs0 = c["ibs0"] / c["n"] if c["n"] else None
        out_pairs.append({"a": names[i], "b": names[j], "sites": c["n"],
                          "kinship": round(kinship, 4) if kinship is not None else None,
                          "ibs0_rate": round(ibs0, 5) if ibs0 is not None else None,
                          "relationship": classify(kinship, ibs0, c["n"])})
    out_trios = [{"child": c, "father": f, "mother": m, "sites": t["n"], "errors": t["errors"],
                  "error_rate": round(t["errors"] / t["n"], 5) if t["n"] else None}
                 for (c, f, m), t in zip(trios, tc)]
    return out_pairs, out_trios


def classify(kinship: float | None, ibs0: float | None, n: int) -> str:
    if kinship is None or n < REL_MIN_SITES:
        return "insufficient data"
    for limit, label in KINSHIP_DEGREES:
        if kinship >= limit:
            if label == "first-degree":
                return "parent-child" if (ibs0 or 0) < PARENT_CHILD_MAX_IBS0 else "full siblings"
            return label
    return "unrelated"


def expected_relationship(a: str, b: str, ped: dict[str, dict]) -> str | None:
    """Relationship the PED file implies for a pair; None when it does not say."""
    pa, pb = ped.get(a), ped.get(b)
    if not pa or not pb:
        return None
    if pa["family"] != pb["family"]:
        return "unrelated"
    if a in (pb["father"], pb["mother"]) or b in (pa["father"], pa["mother"]):
        return "parent-child"
    parents_a = {p for p in (pa["father"], pa["mother"]) if p}
    parents_b = {p for p in (pb["father"], pb["mother"]) if p}
    if parents_a & parents_b:
        if len(parents_a) < 2 or len(parents_b) < 2:
            return None         # one parent not given: full or half siblings
        return "full siblings" if parents_a == parents_b else "second-degree"
    # The father and mother of a common child are taken as unrelated.
    if any(p["father"] in (a, b) and p["mother"] in (a, b) for p in ped.values()):
        return "unrelated"
    return None


def relationship_matches(expected: str, observed: str) -> bool:
    if observed == "insufficient data":
        return True
    if expected == "unrelated":
        return observed in ("unrelated", "third-degree")
    return expected == observed


def sample_qc(file_id: str, stats: str | None, chrx: str | None, genotypes, ped_path: str | None,
              vcf_header: str | None = None) -> dict:
    """
    Per-sample QC of one input file. The pedigree comes from the PED file, else
    from the VCF header (##PEDIGREE), else it is inferred from the relatedness
    of the samples; declared sex and relationships are checked for the first two.
    """
    psc = read_psc(stats)
    xpsc = read_psc(chrx)
    names = list(psc)
    ped, source = match_ped(read_ped(ped_path), names), "ped"
    if not ped:
        ped, source = match_ped(read_vcf_pedigree(vcf_header), names), "vcf header"
    sexes = {s: infer_sex(xpsc.get(s)) for s in names}
    trios = [(s, p["father"], p["mother"]) for s, p in ped.items()
             if s in psc and p["father"] in psc and p["mother"] in psc]
    pairs, trio_stats = relatedness(genotypes, names, trios) if genotypes is not None else ([], [])
    for p in pairs:
        exp = expected_relationship(p["a"], p["b"], ped) if ped else None
        p["expected"] = exp
        p["check"] = None if exp is None else "ok" if relationship_matches(exp, p["relationship"]) else "mismatch"
    if not ped:
        ped, source = infer_pedigree(names, pairs, {s: v[0] for s, v in sexes.items()}), "inferred"
    for entry in ped.values():
        entry["source"] = source

    samples = []
    for s in names:
        p = psc[s]
        sex, ratio, n_x = sexes[s]
        declared = ped.get(s, {}).get("sex", "unknown")
        samples.append({
            "sample": s,
            "snv_het": p.get("nHets"), "snv_hom": p.get("nNonRefHom"),
            "het_hom_ratio": round(p["nHets"] / p["nNonRefHom"], 3) if p.get("nNonRefHom") else None,
            "ts_tv": round(p["nTransitions"] / p["nTransversions"], 3) if p.get("nTransversions") else None,
            "indels": p.get("nIndels"), "mean_depth": p.get("average depth"),
            "singletons": p.get("nSingletons"), "missing": p.get("nMissing"),
            "x_het_hom_ratio": round(ratio, 3) if ratio is not None and ratio != float("inf") else None,
            "x_snvs": n_x, "sex_inferred": sex,
            "ped": ped.get(s),
            "sex_check": (None if declared == "unknown" or sex in ("unknown", "unclear")
                          else "ok" if declared == sex else "mismatch"),
        })
    return {"id": file_id, "samples": samples, "pairs": pairs, "trios": trio_stats,
            "pedigree_source": source if ped else None,
            "pedigree_samples": sorted(set(ped) & set(names)),
            "thresholds": {"relatedness_min_gq": REL_MIN_GQ, "relatedness_min_dp": REL_MIN_DP,
                           "male_max_x_ratio": MALE_MAX_X_RATIO, "female_min_x_ratio": FEMALE_MIN_X_RATIO}}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--id", required=True, help="input-file id")
    p.add_argument("--stats", required=True, help="bcftools stats -s - output (whole file)")
    p.add_argument("--chrx", help="bcftools stats -s - output restricted to chrX outside the PARs")
    p.add_argument("--genotypes", help="'%%CHROM[\\t%%GT:%%GQ:%%DP]' lines for relatedness ('-' for STDIN)")
    p.add_argument("--ped", help="PED file with the samples' families, sex and affection status")
    p.add_argument("--vcf-header",
                   help="VCF header lines: ##PEDIGREE is used when the PED file does not cover the samples")
    p.add_argument("--out", required=True, help="output JSON")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.genotypes == "-":
        result = sample_qc(args.id, args.stats, args.chrx, sys.stdin, args.ped, args.vcf_header)
    elif args.genotypes:
        with open(args.genotypes, encoding="utf-8") as fh:
            result = sample_qc(args.id, args.stats, args.chrx, fh, args.ped, args.vcf_header)
    else:
        result = sample_qc(args.id, args.stats, args.chrx, None, args.ped, args.vcf_header)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)
    sexes = ", ".join(f"{s['sample']}={s['sex_inferred']}" for s in result["samples"][:10])
    print(f"[sample_qc] {args.id}: {len(result['samples'])} sample(s) ({sexes}); "
          f"{len(result['pairs'])} pair(s), {len(result['trios'])} trio(s); "
          f"pedigree: {result['pedigree_source'] or 'none'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ValueError as exc:
        print(f"[sample_qc] ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
