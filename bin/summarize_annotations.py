#!/usr/bin/env python3
"""
summarize_annotations.py
------------------------
Build the annotation summary from the per-chunk statistics JSON files written
by format_annovar_table.py, plus the genome-build reports of
detect_genome_build.py and the list of ANNOVAR databases used (versions). No
annotated table is re-read, so the summary is cheap even for whole-genome
multi-sample cohorts.

Note: an input VCF may hold many samples (a family or a joint-genotyped
cohort). "file" is the input-file id, "sample" a VCF sample; counts are variant
records (sites x alt alleles), not per-genotype counts. Per-sample QC
(sample_qc.py, --sample-qc) brings sex checks, relatedness and the PED data;
the inheritance analysis (inheritance.py) runs on single samples and families.

Outputs, with --out <dir>/<name> (e.g. ./annotation_summary):
    <dir>/<id>.<name>.html            one self-contained report per input file
                                      (cards, samples, P/LP, inheritance
                                      candidates, VUS, genome build, charts,
                                      annotation sources)
    <dir>/<id>.<name>.json            the same content for machines (e.g. a LIS)
    <dir>/<id>.<name>_totals.tsv      totals (SNV/indel, Ti/Tv, PASS, rare, ClinVar ...)
    <dir>/<id>.<name>_samples.tsv     per sample: QC (sex check, depth, het/hom, Ti/Tv) and
                                      what it carries (rare coding, P/LP, VUS, candidates)
    <dir>/<id>.<name>_plp.tsv         carried pathogenic / likely pathogenic variants
                                      (ClinVar, plus InterVar P/LP not P/LP in ClinVar)
    <dir>/<id>.<name>_vus.tsv         carried rare coding/splicing ClinVar VUS / conflicting
    <dir>/<id>.<name>_candidates.tsv  inheritance candidates (de novo, recessive, ...)
    <dir>/<id>.<name>.tsv             Func.refGene x ExonicFunc.refGene counts
    <dir>/<name>.*                    the same across all files (no JSON), only when
                                      there is more than one input file
The variant TSVs have one row per variant and carrier sample, with the call's
zygosity, VAF, depth, GQ, quality flags and confidence, and the HGVS notation
of the first refGene transcript (c.590G>A, p.(Arg197His)).

The HTML report keeps one section order for every input, while its variant
tables follow the input: a single sample gets one genotype column, a family
one column per member and one row per variant, a cohort (or a report of
several files) a list of carriers. Every variant table has the same columns
of the same widths. Doubtful rows - inheritance candidates and VUS whose
calls have low confidence, Mendelian inconsistencies - are hidden behind a
switch; P/LP variants are always shown.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import gene_info
import inheritance

FUNC_PRIORITY = {"exonic": 0, "exonic;splicing": 1, "splicing": 2, "ncRNA_exonic": 3, "UTR5": 4, "UTR3": 5}
REGION_CLASSES = ("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "UTR5", "UTR3",
                  "intronic", "upstream", "downstream", "intergenic")
TOTAL_COLUMNS = ("genome", "n_genotype_columns", "total", "snv", "indel", "mnv", "ti_tv", "pass",
                 *REGION_CLASSES, "rare", "coding_rare", "clinvar_plp", "clinvar_plp_not_carried", "intervar_plp",
                 "vus_candidates", "revel_ge_0.75", "cadd_phred_ge_20")
# Listed variants, one TSV row per variant and carrier sample (file = input-file id, sample = VCF sample).
GENOTYPE_COLUMNS = ("zygosity", "vaf", "dp", "gq", "alt_reads", "flags", "confidence", "inheritance")
VARIANT_COLUMNS = ("gene", "vcf", "chr", "start", "ref", "alt", "func", "exonic", "aachange", "transcript", "hgvs_c",
                   "hgvs_p", "carriers", "filter", "clnsig", "stars", "clnrevstat", "clndn", "clnalleleid", "intervar",
                   "acmg", "gnomad_af", "revel", "cadd", "alphamissense", "splice", "pubmind", "pubmind_score",
                   "pubmind_confidence", "pubmind_pvid", "gene_moi", "gene_diseases", "loeuf", "panels")
PLP_COLUMNS = ("file", "sample", "source", *GENOTYPE_COLUMNS, *VARIANT_COLUMNS)
VUS_COLUMNS = ("file", "sample", "kind", *GENOTYPE_COLUMNS, *VARIANT_COLUMNS)
CANDIDATE_COLUMNS = ("file", "sample", "model", "note", "partners", "doubtful", "moi_fit", *GENOTYPE_COLUMNS[:-1],
                     *VARIANT_COLUMNS)
FIT_TEXT = {True: "yes", False: "no", None: "."}
SAMPLE_COLUMNS = ("file", "sample", "family", "father", "mother", "sex_declared", "affected", "pedigree", "sex_inferred",
                  "x_het_hom_ratio", "sex_check", "mean_depth", "snv_het", "snv_hom", "het_hom_ratio", "ts_tv",
                  "indels", "singletons", "missing", "coding_rare", "clinvar_plp", "intervar_plp", "vus_candidates",
                  "inheritance_candidates")
# What a sample carries, per sample like the totals per file ('NA': the database is not in the protocol).
FINDING_COLUMNS = SAMPLE_COLUMNS[-5:]
CARRIER = ("het", "hom", "hemi")
MAX_ROWS_HTML = 3000          # per table; the TSV files hold every listed variant
MAX_CALLS_HTML = 12           # genotype calls shown per variant (cohorts: carriers only)
RARE_AF = 0.01                # same threshold as format_annovar_table.py
FLAG_TEXT = {
    "LOW_DP": "read depth below {min_dp}", "LOW_GQ": "genotype quality below {min_gq}",
    "LOW_ALT_READS": "fewer than {min_alt_reads} reads support the variant allele",
    "VAF_SKEW": "variant allele fraction outside {min_het_vaf}-{max_het_vaf} (heterozygous) or below "
                "{min_hom_vaf} (homozygous / hemizygous)",
    "MOSAIC": "DRAGEN: possible mosaic variant (FORMAT/MSC)",
    "HOMOLOGY": "DRAGEN: variant in a region of homology, e.g. a paralogous gene (FORMAT/HML)",
    "ALT_READS_IN_REF": "reference call with reads supporting the variant allele",
    "FT": "the variant caller filtered this genotype (FORMAT/FT)",
}

# Region groups for the stacked bar chart: (label, ANNOVAR Func values, CSS colour variable)
REGION_GROUPS = [
    ("coding", ("exonic", "exonic;splicing"), "--c1"),
    ("splicing", ("splicing",), "--c2"),
    ("UTR", ("UTR5", "UTR3", "UTR5;UTR3"), "--c3"),
    ("ncRNA", ("ncRNA_exonic", "ncRNA_intronic", "ncRNA_splicing", "ncRNA_exonic;splicing", "ncRNA_UTR5",
               "ncRNA_UTR3"), "--c4"),
    ("intronic", ("intronic",), "--c5"),
    ("up/downstream", ("upstream", "downstream", "upstream;downstream"), "--c6"),
    ("intergenic", ("intergenic",), "--c7"),
]
CONSEQUENCES = ("nonsynonymous SNV", "synonymous SNV", "stopgain", "stoploss", "startloss",
                "frameshift insertion", "frameshift deletion", "nonframeshift insertion",
                "nonframeshift deletion", "unknown")
# ANNOVAR exonic functions in Sequence Ontology words, for the report.
CONSEQUENCE_LABELS = {
    "nonsynonymous SNV": "missense", "synonymous SNV": "synonymous", "stopgain": "stop gained",
    "stoploss": "stop lost", "startloss": "start lost", "frameshift insertion": "frameshift",
    "frameshift deletion": "frameshift", "frameshift substitution": "frameshift",
    "frameshift block substitution": "frameshift", "nonframeshift insertion": "inframe insertion",
    "nonframeshift deletion": "inframe deletion", "nonframeshift substitution": "inframe indel",
}
AMINO_ACIDS = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
               "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
               "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "U": "Sec", "X": "Ter", "*": "Ter"}
TRANSCRIPT_RE = re.compile(r"^(N[MR]|X[MR])_\d|^ENST\d|^uc\d")
# How detect_genome_build.py decided the build ("source" in its report).
BUILD_SOURCES = {"explicit": "set by user (--genome)", "contig_lengths": "##contig lengths",
                 "header_hints": "header hints", "filename": "file name"}


# Databases used (databases.tsv written by the pipeline, one row per protocol entry).
DATABASE_COLUMNS = ("genome", "protocol", "database", "operation", "file", "file_date", "bytes")
OPERATIONS = {"g": "gene-based", "gx": "gene-based", "r": "region", "f": "filter"}
# Versionless databases: readable source names.
DATABASE_SOURCES = {
    "refGene": "RefSeq genes (UCSC refGene)", "knownGene": "UCSC knownGene (GENCODE)",
    "ensGene": "Ensembl genes (UCSC ensGene)", "kgXref": "UCSC kgXref", "revel": "REVEL", "mcap": "M-CAP",
    "nci60": "NCI-60", "gme": "GME Variome", "abraom": "ABraOM", "cg69": "Complete Genomics 69",
    "regsnpintron": "regSNP-intron", "gerp++elem": "GERP++ elements", "gerp++gt2": "GERP++ (RS > 2)",
    "fathmm": "FATHMM", "eigen": "Eigen", "gwava": "GWAVA",
}


def _date(d: str) -> str:
    return f"{d[:4]}-{d[4:6]}-{d[6:]}"


# (pattern, source, version) - source / version are format strings over the match groups.
DATABASE_PATTERNS = [
    (r"clinvar_(\d{8})", "ClinVar", lambda m: _date(m[1])),
    (r"intervar_(\d{8})", "InterVar", lambda m: _date(m[1])),
    (r"dbnsfp(\d)(\d)([a-z])(_interpro)?", lambda m: "dbNSFP InterPro domains" if m[4] else "dbNSFP",
     lambda m: f"{m[1]}.{m[2]}{m[3]}"),
    (r"gnomad(\d*)_(exome|genome)", lambda m: f"gnomAD {m[2]}s", lambda m: ".".join(m[1])),
    (r"avsnp(\d+)", "dbSNP", lambda m: m[1]),
    (r"exac(\d)(\d)(nontcga|nonpsych)?", lambda m: f"ExAC ({m[3]})" if m[3] else "ExAC", lambda m: f"{m[1]}.{m[2]}"),
    (r"([A-Z]{3})\.sites\.(\d{4})_(\d{2})", lambda m: f"1000 Genomes ({m[1]})", lambda m: f"{m[2]}-{m[3]}"),
    (r"esp6500siv2_(\w+)", lambda m: f"NHLBI ESP6500 ({m[1].upper()})", lambda m: "SI-V2"),
    (r"cosmic(\d+)", "COSMIC", lambda m: f"v{m[1]}"),
    (r"cadd(\d)(\d)", "CADD", lambda m: f"{m[1]}.{m[2]}"),
    (r"dbscsnv(\d)(\d)", "dbscSNV", lambda m: f"{m[1]}.{m[2]}"),
    (r"mitimpact(\d)(\d)", "MitImpact", lambda m: f"{m[1]}.{m[2]}"),
    (r"icgc(\d+)", "ICGC", lambda m: f"release {m[1]}"),
    (r"kaviar_(\d{8})", "Kaviar", lambda m: _date(m[1])),
    (r"popfreq_all_(\d{8})", "ANNOVAR popfreq (maximum over populations)", lambda m: _date(m[1])),
    (r"hrcr(\d+)", "Haplotype Reference Consortium", lambda m: f"r{m[1]}"),
    (r"pubmind_db(?:_v(\d)(\d)(\d))?", "PubMind", lambda m: f"{m[1]}.{m[2]}.{m[3]}" if m[1] else ""),
]
# Sources shown with their version in the report header.
KEY_SOURCES = ("ClinVar", "gnomAD", "dbNSFP", "InterVar", "dbSNP")


def describe_database(name: str) -> tuple[str, str]:
    """ANNOVAR database name -> (source, version), e.g. gnomad41_genome -> ('gnomAD genomes', '4.1')."""
    for pattern, source, version in DATABASE_PATTERNS:
        m = re.fullmatch(pattern, name)
        if m:
            return (source(m) if callable(source) else source), version(m)
    return DATABASE_SOURCES.get(name, name), ""


def load_databases(path: str | None) -> list[dict]:
    """Rows of databases.tsv ('#' lines and blank lines skipped)."""
    if not path:
        return []
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                rows.append(dict(zip(DATABASE_COLUMNS, line.rstrip("\n").split("\t"))))
    return rows


def review_stars(status: str) -> int:
    """ClinVar review status -> gold stars (0-4)."""
    s = status.lower().replace("\\x2c", ",")
    if "practice_guideline" in s:
        return 4
    if "expert_panel" in s:
        return 3
    if "multiple_submitters" in s and "no_conflicts" in s:
        return 2
    if "single_submitter" in s or "conflicting" in s:
        return 1
    return 0


def hgvs_c(c: str) -> str:
    """ANNOVAR cDNA change -> HGVS: c.G590A -> c.590G>A, c.590_592delCTT -> c.590_592del."""
    m = re.fullmatch(r"c\.([ACGTN])([-*]?\d+(?:[+-]\d+)?)([ACGTN])", c)
    if m:
        return f"c.{m[2]}{m[1]}>{m[3]}"
    return re.sub(r"(del|dup)[ACGTN]+$", r"\1", c)


def hgvs_p(p: str) -> str:
    """ANNOVAR protein change -> HGVS, predicted, three-letter codes: p.R197H -> p.(Arg197His),
    p.R12X -> p.(Arg12Ter), p.G1355Lfs*167 -> p.(Gly1355LeufsTer167), p.L120L -> p.(Leu120=)."""
    body = p.removeprefix("p.")
    if body in ("", ".", "?"):
        return p
    m = re.fullmatch(r"([A-Z*])(\d+)\1", body)
    if m:
        return f"p.({AMINO_ACIDS.get(m[1], m[1])}{m[2]}=)"
    return "p.(" + "".join(AMINO_ACIDS.get(ch, ch) if ch.isupper() or ch == "*" else ch for ch in body) + ")"


def hgvs(v: dict) -> dict:
    """Transcript and HGVS c. / p. notation of a listed variant, from its first refGene
    AAChange ('GENE:NM_1:exon5:c.G590A:p.R197H') or, for splice sites, GeneDetail."""
    for value in (v.get("aachange"), v.get("genedetail")):
        tokens = re.split(r"[,;]", str(value or "."))[0].split(":")      # the first transcript
        c = next((t for t in tokens if t.startswith("c.")), "")
        p = next((t for t in tokens if t.startswith("p.")), "")
        if c or p:
            tx = next((t for t in tokens if TRANSCRIPT_RE.match(t)), ".")
            return {"transcript": tx, "hgvs_c": hgvs_c(c) if c else ".", "hgvs_p": hgvs_p(p) if p else "."}
    return {"transcript": ".", "hgvs_c": ".", "hgvs_p": "."}


def consequence(v: dict) -> str:
    """'missense', 'stop gained, splice site', ... (ANNOVAR Func / ExonicFunc in Sequence Ontology words)."""
    words = [CONSEQUENCE_LABELS.get(v.get("exonic", "."), v.get("exonic"))] if v.get("exonic") not in (None, ".", "") \
        else []
    if inheritance.splicing(v):
        words.append("splice site")
    return ", ".join(words) or str(v.get("func", "."))


def load(paths: list[str]) -> tuple[dict, dict, dict, list, list, dict]:
    """Aggregate chunk statistics per input file (the stats' 'sample' is the input-file id)."""
    counts: dict[str, Counter] = defaultdict(Counter)
    func: dict[str, Counter] = defaultdict(Counter)
    meta: dict[str, dict] = {}
    plp: list[dict] = []
    vus: list[dict] = []
    rare: dict[str, list] = defaultdict(list)
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        file_id = data["sample"]
        counts[file_id].update(data.get("counts", {}))
        func[file_id].update(data.get("func", {}))
        m = meta.setdefault(file_id, {"n_genotype_columns": data.get("n_genotype_columns", 0), "genome": "",
                                      "columns": {}, "samples": [], "qc": {}, "sample_counts": None})
        if data.get("genome"):
            m["genome"] = data["genome"]
        if data.get("samples"):
            m["samples"] = data["samples"]
        if "sample_counts" in data:       # None: statistics of an older version, without per-sample counts
            m["sample_counts"] = m["sample_counts"] or {}
            for s, n in data["sample_counts"].items():
                m["sample_counts"].setdefault(s, Counter()).update(n)
        if data.get("qc"):
            m["qc"] = data["qc"]
        for key, col in (data.get("columns") or {}).items():
            if col:
                m["columns"][key] = col
        for v in data.get("plp", []):
            plp.append({"file": file_id, **v, **hgvs(v), "stars": review_stars(v.get("clnrevstat", ""))})
        for v in data.get("vus", []):
            vus.append({"file": file_id, **v, **hgvs(v), "stars": review_stars(v.get("clnrevstat", ""))})
        for v in data.get("rare", []):
            rare[file_id].append({"file": file_id, **v, **hgvs(v), "stars": review_stars(v.get("clnrevstat", ""))})
    # P/LP: ClinVar before InterVar-only, best reviewed first; VUS: most damaging REVEL first.
    plp.sort(key=lambda v: (v["file"], v.get("source") != "ClinVar", -v["stars"], v.get("gene", ""),
                            v.get("chr", ""), int(v.get("start") or 0)))
    vus.sort(key=lambda v: (v["file"], -(to_float(v.get("revel")) or -1), v.get("gene", "")))
    return counts, func, meta, plp, vus, rare


def load_sample_qc(paths: list[str]) -> dict[str, dict]:
    """<id>.sample_qc.json reports of sample_qc.py keyed by input-file id."""
    out = {}
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        out[data.get("id") or Path(path).name.split(".sample_qc")[0]] = data
    return out


def pedigree(qc: dict) -> tuple[dict, dict]:
    """PED entries and inferred sex of a file's samples, from its sample QC report."""
    samples = qc.get("samples", [])
    return ({s["sample"]: s["ped"] for s in samples if s.get("ped")},
            {s["sample"]: s.get("sex_inferred", "unknown") for s in samples})


def family_analysis(files: list[str], meta: dict, rare: dict, sample_qc: dict) -> tuple[list, dict]:
    """Inheritance candidates of every single-sample / family file, and the models per (file, sample, variant)."""
    candidates, labels = [], defaultdict(list)
    for f in files:
        if not rare.get(f) or not meta[f]["samples"]:
            continue
        ped, sex = pedigree(sample_qc.get(f, {}))
        for c in inheritance.analyse(rare[f], meta[f]["samples"], ped, sex, meta[f]["genome"]):
            c["file"] = f
            candidates.append(c)
            key = (f, c["sample"], c["variant"].get("vcf"))
            if c["model"] not in labels[key]:
                labels[key].append(c["model"])
    return candidates, labels


def num(value, digits: int = 2) -> str:
    if value is None or value == "":
        return "."
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def call_fields(g: dict) -> dict:
    """TSV fields of one genotype call."""
    return {"sample": g.get("sample", "."), "zygosity": g.get("zyg", "."), "vaf": num(g.get("vaf")),
            "dp": num(g.get("dp")), "gq": num(g.get("gq")), "alt_reads": num(g.get("alt_reads")),
            "flags": ",".join(g.get("flags") or []) or ".", "confidence": g.get("conf", ".")}


def carrier_rows(variants: list[dict], labels: dict) -> list[dict]:
    """One row per listed variant and carrier sample (one row without sample for sites-only VCFs)."""
    rows = []
    for v in variants:
        calls = [g for g in v.get("genotypes", []) if g.get("zyg") in CARRIER] or [{}]
        for g in calls:
            row = {**v, **call_fields(g)}
            row["inheritance"] = "; ".join(labels.get((v["file"], g.get("sample"), v.get("vcf")), [])) or "."
            rows.append(row)
    return rows


def candidate_rows(candidates: list[dict]) -> list[dict]:
    rows = []
    for c in candidates:
        g = next((x for x in c["variant"].get("genotypes", []) if x.get("sample") == c["sample"]), {})
        rows.append({**c["variant"], **call_fields(g), "file": c["file"], "sample": c["sample"],
                     "model": c["model"], "note": c["note"], "partners": ",".join(p for p in c["partners"] if p) or ".",
                     "doubtful": c.get("doubtful") or ".",
                     "moi_fit": FIT_TEXT[gene_info.moi_fit(c["model"], c["variant"].get("_gene"))]})
    return rows


def annotate_genes(variants, genes: gene_info.GeneTable, panels: list) -> None:
    """Gene information (diseases, modes of inheritance, LOEUF) and the gene panels of every listed variant."""
    for v in variants:
        info = genes.lookup(v.get("gene", ""))
        names = gene_info.panels_of(v.get("gene", ""), panels, genes) if panels else []
        v.update({"_gene": info, "_panel": bool(names), "panels": ",".join(names) or ".",
                  "gene_moi": ",".join(info.moi) if info and info.moi else ".",
                  "gene_diseases": (info.disease_text() if info else "") or ".",
                  "loeuf": f"{info.loeuf:.2f}" if info and info.loeuf is not None else "."})


def sample_rows(files: list[str], sample_qc: dict) -> list[dict]:
    rows = []
    for f in files:
        for s in sample_qc.get(f, {}).get("samples", []):
            ped = s.get("ped") or {}
            affected = {True: "affected", False: "unaffected"}.get(ped.get("affected"), ".")
            rows.append({**s, "file": f, "family": ped.get("family", "."), "father": ped.get("father") or ".",
                         "mother": ped.get("mother") or ".", "sex_declared": ped.get("sex", "."), "affected": affected,
                         "pedigree": ped.get("source", ".")})
    return rows


def sample_findings(files: list[str], meta: dict, candidates: list[dict], family_max: int) -> dict:
    """(file, sample) -> what the sample carries: rare coding variants, ClinVar P/LP, InterVar-only P/LP,
    rare coding VUS and its inheritance candidates (variants, doubtful ones left out; 'models' counts
    them per model). 'NA' when the database is not in the protocol or the file is not analysed."""
    out = {}
    for f in files:
        counts = meta[f].get("sample_counts")
        cols = meta[f]["columns"] if counts is not None else {}
        analysed = 0 < len(meta[f]["samples"]) <= family_max
        for s in meta[f]["samples"]:
            n = (counts or {}).get(s, {})
            mine = [c for c in candidates if c["file"] == f and c["sample"] == s and not c.get("doubtful")]
            models = Counter(m for m, _ in {(c["model"], c["variant"].get("vcf")) for c in mine})
            out[(f, s)] = {
                "coding_rare": n.get("coding_rare", 0) if cols.get("gnomad_af") else "NA",
                "clinvar_plp": n.get("clinvar_plp", 0) if cols.get("clnsig") else "NA",
                "intervar_plp": n.get("intervar_plp", 0) if cols.get("intervar") else "NA",
                "vus_candidates": n.get("vus_candidates", 0) if cols.get("clnsig") else "NA",
                "inheritance_candidates": len({c["variant"].get("vcf") for c in mine}) if analysed else "NA",
                "models": {m: models[m] for m in inheritance.MODEL_ORDER if models.get(m)},
            }
    return out


def sample_table(files: list[str], sample_qc: dict, findings: dict) -> list[dict]:
    """_samples.tsv rows: every sample with its QC (when run) and what it carries."""
    qc = {(r["file"], r["sample"]): r for r in sample_rows(files, sample_qc)}
    keys = list(dict.fromkeys([*qc, *findings]))
    return [{"file": f, "sample": s, **qc.get((f, s), {}), **findings.get((f, s), {})} for f, s in keys]


def to_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_builds(paths: list[str]) -> dict[str, dict]:
    """Genome-build reports ('<id>.build.json') keyed by input-file id."""
    builds = {}
    for path in paths:
        name = Path(path).name
        sample = name[: -len(".build.json")] if name.endswith(".build.json") else Path(path).stem
        with open(path, encoding="utf-8") as fh:
            builds[sample] = json.load(fh)
    return builds


def wide_matrix(func: dict[str, Counter], samples: list[str]) -> list[list]:
    keys = set()
    for c in func.values():
        keys.update(c)

    def order(key):
        f, e = key.split("\t", 1)
        return (FUNC_PRIORITY.get(f, 99), f, e)

    rows = []
    for key in sorted(keys, key=order):
        f, e = key.split("\t", 1)
        rows.append([f, e] + [func[s].get(key, 0) for s in samples])
    return rows


def totals_table(counts: dict[str, Counter], func: dict[str, Counter], meta: dict, samples: list[str]) -> list[dict]:
    rows = []
    for s in samples:
        c = counts[s]
        by_region = Counter()
        for key, n in func[s].items():
            by_region[key.split("\t", 1)[0]] += n
        cols = meta[s]["columns"]
        row = {
            "file": s,
            "genome": meta[s]["genome"] or "NA",
            "n_genotype_columns": meta[s]["n_genotype_columns"],
            "total": c.get("total", 0),
            "snv": c.get("snv", 0),
            "indel": c.get("indel", 0),
            "mnv": c.get("mnv", 0),
            "ti_tv": f"{c['ti'] / c['tv']:.3f}" if c.get("tv") else "NA",
            "pass": c.get("pass", 0),
        }
        for region in REGION_CLASSES:
            row[region] = by_region.get(region, 0)
        # 'NA' when the database was not part of the protocol.
        row["rare"] = c.get("rare", 0) if cols.get("gnomad_af") else "NA"
        row["coding_rare"] = c.get("coding_rare", 0) if cols.get("gnomad_af") else "NA"
        row["clinvar_plp"] = c.get("clinvar_plp", 0) if cols.get("clnsig") else "NA"
        row["clinvar_plp_not_carried"] = c.get("clinvar_plp_not_carried", 0) if cols.get("clnsig") else "NA"
        row["intervar_plp"] = c.get("intervar_plp", 0) if cols.get("intervar") else "NA"
        row["vus_candidates"] = c.get("vus_candidates", 0) if cols.get("clnsig") else "NA"
        row["revel_ge_0.75"] = c.get("revel_ge_0.75", 0) if cols.get("revel") else "NA"
        row["cadd_phred_ge_20"] = c.get("cadd_phred_ge_20", 0) if cols.get("cadd") else "NA"
        rows.append(row)
    return rows


def write_tsv(path: str, header, rows) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("." if x is None else str(x) for x in r) + "\n")


# ======================================================================================
# HTML
# ======================================================================================
def esc(x) -> str:
    return html.escape(str(x))


def fmt(x) -> str:
    return f"{x:,}" if isinstance(x, int) else esc(x)


def thead(head, numeric=()) -> str:
    """Header row (labels are HTML); the headers of numeric columns are right-aligned like their cells."""
    return ("<thead><tr>" + "".join(f"<th class='num'>{h}</th>" if h in numeric else f"<th>{h}</th>" for h in head)
            + "</tr></thead>")


def html_table(header, rows, n_label_cols: int, table_id: str = "") -> str:
    th = "".join(f"<th class='num'>{esc(h)}</th>" if i >= n_label_cols else f"<th>{esc(h)}</th>"
                 for i, h in enumerate(header))
    body = "".join(
        "<tr>" + "".join(f'<td class="{"label" if i < n_label_cols else "num"}">{fmt(x)}</td>'
                         for i, x in enumerate(r)) + "</tr>" for r in rows)
    ident = f' id="{table_id}"' if table_id else ""
    return f"<table{ident}><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>"


def file_cards(totals: list[dict], builds: dict[str, dict]) -> str:
    cards = []
    for t in totals:
        plp = t["clinvar_plp"]
        plp_cls = "warn" if isinstance(plp, int) and plp > 0 else ""
        source = builds.get(t["file"], {}).get("source")
        build_from = f"<dt>build from</dt><dd>{esc(BUILD_SOURCES.get(source, source))}</dd>" if source else ""
        cards.append(f"""
<div class="card">
  <div class="card-head"><span class="name">{esc(t['file'])}</span><span class="badge">{esc(t['genome'])}</span></div>
  <div class="big">{fmt(t['total'])}<span> variants</span></div>
  <dl>
    {build_from}
    <dt>SNV / indel</dt><dd>{fmt(t['snv'])} / {fmt(t['indel'])}</dd>
    <dt>Ti/Tv</dt><dd>{esc(t['ti_tv'])}</dd>
    <dt>PASS</dt><dd>{fmt(t['pass'])}</dd>
    <dt>samples in file</dt><dd>{fmt(t['n_genotype_columns'])}</dd>
    <dt>rare coding</dt><dd>{fmt(t['coding_rare'])}</dd>
    <dt>ClinVar P/LP carried</dt><dd class="{plp_cls}">{fmt(plp)}</dd>
    <dt>InterVar-only P/LP</dt><dd>{fmt(t['intervar_plp'])}</dd>
    <dt>rare coding VUS</dt><dd>{fmt(t['vus_candidates'])}</dd>
  </dl>
</div>""")
    return "".join(cards)


def region_chart(func: dict[str, Counter], samples: list[str]) -> str:
    rows = []
    for s in samples:
        by_region = Counter()
        for key, n in func[s].items():
            by_region[key.split("\t", 1)[0]] += n
        total = sum(by_region.values()) or 1
        segments, used = [], 0
        for label, values, colour in REGION_GROUPS:
            n = sum(by_region.get(v, 0) for v in values)
            used += n
            if n:
                segments.append((label, n, colour))
        if total - used:
            segments.append(("other", total - used, "--c8"))
        bar = "".join(f'<span style="width:{n / total * 100:.3f}%;background:var({c})" '
                      f'title="{esc(label)}: {n:,} ({n / total:.1%})"></span>' for label, n, c in segments)
        rows.append(f'<div class="bar-row"><div class="bar-label">{esc(s)}</div><div class="stack">{bar}</div></div>')
    legend = "".join(f'<span><i style="background:var({c})"></i>{esc(label)}</span>'
                     for label, _, c in REGION_GROUPS + [("other", (), "--c8")])
    return f'<div class="legend">{legend}</div>{"".join(rows)}'


def consequence_chart(func: dict[str, Counter], samples: list[str]) -> str:
    per = {s: Counter() for s in samples}
    for s in samples:
        for key, n in func[s].items():
            f, e = key.split("\t", 1)
            if f in ("exonic", "exonic;splicing") and e not in (".", ""):
                per[s][e] += n
    labels = [c for c in CONSEQUENCES if any(per[s].get(c) for s in samples)]
    labels += sorted({e for s in samples for e in per[s]} - set(labels))
    if not labels:
        return '<p class="muted">No exonic variants.</p>'
    peak = max((per[s].get(c, 0) for s in samples for c in labels), default=1) or 1
    head = "".join(f"<th>{esc(s)}</th>" for s in samples)
    body = []
    for c in labels:
        cells = []
        for s in samples:
            n = per[s].get(c, 0)
            cells.append(f'<td><div class="hbar"><span style="width:{n / peak * 100:.2f}%"></span>'
                         f"<em>{n:,}</em></div></td>")
        body.append(f'<tr><td class="label">{esc(c)}</td>{"".join(cells)}</tr>')
    return f'<table class="conseq"><thead><tr><th>consequence</th>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'


def build_table(builds: dict[str, dict], samples: list[str]) -> str:
    """How the genome build of every input file was determined."""
    rows = []
    for s in samples:
        b = builds.get(s)
        if not b:
            continue
        ev = b.get("evidence") or {}
        votes = ev.get("contig_length_votes") or {}
        vote_text = " · ".join(f"{k} {v}" for k, v in sorted(votes.items(), key=lambda kv: -kv[1])) or "none"
        warnings = "<br>".join(esc(w) for w in b.get("warnings") or []) or '<span class="muted">none</span>'
        rows.append(
            f"<tr><td class='label'>{esc(s)}</td><td class='mono'>{esc(b.get('filename', ''))}</td>"
            f"<td><span class='badge'>{esc(b.get('build', 'NA'))}</span></td>"
            f"<td>{esc(BUILD_SOURCES.get(b.get('source'), b.get('source', '')))}</td>"
            f"<td>{esc(b.get('expected', ''))}</td><td>{esc(vote_text)}</td>"
            f"<td>{esc(', '.join(ev.get('header_hints') or []) or 'none')}</td>"
            f"<td>{esc(', '.join(ev.get('filename_tokens') or []) or 'none')}</td>"
            f"<td>{warnings}</td></tr>")
    if not rows:
        return '<p class="muted">No genome-build reports.</p>'
    head = thead(("file", "input", "build", "determined by", "requested", "##contig length votes", "header hints",
                  "file-name tokens", "warnings"))
    return f'<div class="wrap"><table>{head}<tbody>{"".join(rows)}</tbody></table></div>'


def sig_badge(label: str, title: str = "") -> str:
    """Coloured classification badge (ClinVar CLNSIG, InterVar or PubMind call): short text, the full label
    (or *title*) as tooltip."""
    if label in (".", ""):
        return "<span class='muted'>·</span>"
    full = label.replace("_", " ")
    first = full.split("|")[0]
    low = first.lower()
    cls = ("vus" if low.startswith(("uncertain", "conflicting")) else "other" if "benign" in low
           else "path" if low.startswith("pathogenic") else "lpath" if low.startswith("likely pathogenic") else "other")
    text = ("Conflicting" if low.startswith("conflicting") else "VUS" if low.startswith("uncertain")
            else "P/LP" if low.startswith("pathogenic/likely pathogenic")
            else "B/LB" if low.startswith("benign/likely benign")
            else "Likely path." if low.startswith("likely pathogenic") else first)
    return f"<span class='sig {cls}' title='{esc(title or full)}'>{esc(text)}</span>"


def pubmind_cell(v: dict) -> str:
    """PubMind call, paper-level score and evidence confidence; the PubMind variant ids in the tooltip."""
    label = v.get("pubmind", ".")
    if label in (".", "", None):
        return "<td><span class='muted'>·</span></td>"
    score, conf = v.get("pubmind_score", "."), v.get("pubmind_confidence", ".")
    pvids = (v.get("pubmind_pvid") or ".").split("|")
    title = (f"PubMind (literature, extracted by a language model): {label}"
             + (f" · paper-level pathogenicity score {score} (0 benign, 1 pathogenic)" if score != "." else "")
             + (f" · evidence confidence {conf} of 3" if conf != "." else "")
             + (f" · {', '.join(pvids)}" if pvids != ["."] else ""))
    detail = " · ".join(x for x in (score if score != "." else "", f"conf {conf}/3" if conf != "." else "") if x)
    return (f"<td>{sig_badge(label, title)}<div class='muted small' title='{esc(title)}'>{esc(detail)}</div></td>")


def variant_links(v: dict, genome: str) -> list[list[str]]:
    """[label, url] of the pages about a variant (VCF coordinates) and its gene, for the link menu.
    IGV: the locus in a running IGV desktop (batch port 60151)."""
    out = []
    if v.get("clnalleleid") not in (None, ".", ""):
        out.append(["ClinVar", f"https://www.ncbi.nlm.nih.gov/clinvar/?term={v['clnalleleid']}[alleleid]"])
    parts = (v.get("vcf") or "").split("-")
    if len(parts) == 4 and "." not in parts and genome in ("hg19", "hg38"):
        vcf_chrom, pos, ref, alt = parts
        chrom = vcf_chrom.removeprefix("chr").replace("MT", "M")
        dataset = "gnomad_r2_1" if genome == "hg19" else "gnomad_r4"
        out += [["gnomAD", f"https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}?dataset={dataset}"],
                ["Franklin", f"https://franklin.genoox.com/clinical-db/variant/snp/chr{chrom}-{pos}-{ref}-{alt}"
                             + ("-hg38" if genome == "hg38" else "")],
                ["VarSome", f"https://varsome.com/variant/{genome}/chr{chrom}-{pos}-{ref}-{alt}"]]
    genes = gene_info.genes_of(v.get("gene"))
    info = v.get("_gene")
    if genes:
        out.append(["OMIM", f"https://www.omim.org/search?search={genes[0]}"])
    if info and info.hgnc_id.startswith("HGNC:"):
        out.append(["GenCC", f"https://search.thegencc.org/genes/{info.hgnc_id}"])
    if len(parts) == 4 and parts[1].isdigit():
        out.append(["IGV", f"http://localhost:60151/goto?locus={parts[0]}:{parts[1]}"])
    return out


def gene_cell(v: dict, view: View) -> str:
    """Gene name (opens the link menu) and its modes of inheritance and LOEUF; the diseases on hover."""
    info = v.get("_gene")
    links = esc(json.dumps(variant_links(v, view.genomes.get(v["file"], ""))))
    line = ""
    if info:
        moi = "".join(f"<span class='moi'>{esc(m)}</span>" for m in info.moi)
        loeuf = ("" if info.loeuf is None else
                 f"<span class='loeuf{' lof' if info.loeuf < gene_info.LOEUF_CONSTRAINED else ''}'>"
                 f"LOEUF {info.loeuf:.2f}</span>")
        line = f"<div class='small gene-info' title='{esc(info.detail())}'>{moi}{loeuf}</div>"
    off = "<div class='small muted'>outside the panel</div>" if view.panels and not v.get("_panel") else ""
    return (f"<td><button type='button' class='gene' data-links='{links}' title='{esc(v['gene'])}: links'>"
            f"{esc(v['gene'])}</button>{line or off}</td>")


def condition_cell(v: dict) -> str:
    """ClinVar condition of the variant, then the GenCC diseases of the gene."""
    condition = v.get("clndn", ".").replace("|", "; ")
    info = v.get("_gene")
    first = (f"<div class='clamp' title='{esc(condition)}'>{esc(condition)}</div>" if not (info and info.diseases)
             else f"<div title='{esc(condition)}'>{esc(condition) if condition != '.' else '·'}</div>")
    if not (info and info.diseases):
        return f"<td>{first}</td>"
    d = info.diseases[0]
    more = f" +{len(info.diseases) - 1}" if len(info.diseases) > 1 else ""
    return (f"<td>{first}<div class='muted small' title='{esc(info.disease_text())}'>GenCC: {esc(d['title'])} "
            f"({esc(d['moi'])}){more}</div></td>")


def fit_chip(model: str, info) -> str:
    """Does the model fit the gene's modes of inheritance (GenCC)? Nothing when the gene has none."""
    fit = gene_info.moi_fit(model, info)
    if fit is None:
        return ""
    moi = "/".join(info.moi)
    return (f" <span class='fit {'ok' if fit else 'no'}' title='gene: {esc(moi)} (GenCC) - "
            f"{'fits' if fit else 'does not fit'} the model'>{esc(moi)}</span>")


CONF_RANK = {"low": 2, "caution": 1, "high": 0}


def flag_text(flag: str, qc: dict) -> str:
    key = "FT" if flag.startswith("FT:") else flag
    try:
        return FLAG_TEXT.get(key, flag).format(**qc)
    except (KeyError, IndexError):
        return FLAG_TEXT.get(key, flag)


def roles(files: list[str], sample_qc: dict) -> dict[tuple[str, str], str]:
    """Short role of every sample with PED data: father / mother / affected / unaffected."""
    out = {}
    for f in files:
        ped, _ = pedigree(sample_qc.get(f, {}))
        parents = {p["father"]: "father" for p in ped.values() if p.get("father")}
        parents.update({p["mother"]: "mother" for p in ped.values() if p.get("mother")})
        for s, p in ped.items():
            status = {True: "affected", False: "unaffected"}.get(p.get("affected"), "")
            if not status and (p.get("father") or p.get("mother")):
                status = "child"
            source = {"vcf header": "VCF", "inferred": "inferred"}.get(p.get("source"), "")
            out[(f, s)] = " · ".join(x for x in (parents.get(s, ""), status) if x) + (f" ({source})" if source else "")
    return out


PEDIGREE_SOURCES = {
    "ped": "from the PED file (--ped)",
    "vcf header": "from the VCF header (##PEDIGREE, e.g. DRAGEN joint calling); affected status unknown: "
                  "the children are analysed as probands",
    "inferred": "inferred from the relatedness of the samples (no --ped, no ##PEDIGREE); affected status "
                "unknown: the children are analysed as probands",
}
# Inheritance models as the narrow model column shows them (the full name is the tooltip).
MODEL_LABELS = {"homozygous (AR)": "homozygous (AR)", "X-linked (hemizygous)": "X-linked hemizygous",
                "compound heterozygous": "compound het", "possible compound heterozygous": "compound het (possible)",
                "dominant (affected parent)": "dominant (parent)", "de novo": "de novo"}
# Variant tables: one column set for P/LP, candidates and VUS, with the same
# widths (px) in every table; 'change' and 'condition' share the rest of the
# page and get at least FLEX_MIN each (narrower windows scroll the table).
COLUMN_WIDTHS = {"file": 150, "lead": 150, "gene": 104, "genotype": 168, "member": 112, "carriers": 190,
                 "ClinVar": 112, "InterVar · ACMG": 112, "PubMind": 104, "gnomAD AF": 76, "REVEL": 58, "CADD": 54,
                 "AlphaMissense": 98}
NUMERIC_COLUMNS = ("gnomAD AF", "REVEL", "CADD")
FLEX_MIN = 150


@dataclass
class View:
    """How the variant tables of one report look: the genotype layout follows the input."""
    layout: str                   # single | family | cohort | sites (see report_layout)
    members: list                 # family: the samples, one genotype column each
    sexes: dict                   # sample -> male / female (PED, else inferred)
    role: dict                    # (file, sample) -> role
    labels: dict                  # (file, sample, vcf) -> inheritance models
    qc: dict                      # genotype-quality thresholds
    genomes: dict                 # file -> build
    show_file: bool
    pubmind: bool = False         # PubMind column: the database is in the protocol of a file
    panels: list = field(default_factory=list)      # gene panels (--gene_panel): other genes hidden


def report_layout(files: list[str], meta: dict, family_max: int) -> str:
    """'single' (one sample), 'family' (2..family_max samples), 'cohort' (more samples, or a
    report of several input files) or 'sites' (a sites-only VCF)."""
    if len(files) != 1:
        return "cohort"
    n = len(meta[files[0]]["samples"])
    return "sites" if n == 0 else "single" if n == 1 else "family" if n <= family_max else "cohort"


def make_view(files: list[str], run: Run) -> View:
    layout = report_layout(files, run.meta, run.family_max)
    members = run.meta[files[0]]["samples"] if layout == "family" else []
    ped, inferred = pedigree(run.sample_qc.get(files[0], {})) if len(files) == 1 else ({}, {})
    sexes = {s: (ped.get(s) or {}).get("sex") if (ped.get(s) or {}).get("sex") in ("male", "female")
             else inferred.get(s, "") for s in members}
    qc = next((run.meta[f]["qc"] for f in files if run.meta[f].get("qc")), {})
    genomes = {f: run.meta[f]["genome"] for f in files}
    pubmind = any(run.meta[f]["columns"].get("pubmind") for f in files)
    return View(layout, members, sexes, roles(files, run.sample_qc), run.labels, qc, genomes, len(files) > 1,
                pubmind, run.panels)


def worst_confidence(v: dict) -> str:
    confs = [g.get("conf", "high") for g in v.get("genotypes", []) if g.get("zyg") in CARRIER]
    return max(confs, key=lambda c: CONF_RANK.get(c, 0)) if confs else "high"


def best_confidence(v: dict) -> str:
    confs = [g.get("conf", "high") for g in v.get("genotypes", []) if g.get("zyg") in CARRIER]
    return min(confs, key=lambda c: CONF_RANK.get(c, 0)) if confs else "high"


def mito(v: dict) -> bool:
    return str(v.get("chr", "")).removeprefix("chr") in ("M", "MT")


def percent(x: float) -> str:
    return f"{x:.0%}" if x >= 0.01 or x == 0 else f"{x:.1%}"


def call_html(g: dict, qc: dict, on_mito: bool, name: str = "") -> str:
    """One genotype call: zygosity (heteroplasmy on chrM), VAF, depth and quality flags;
    GQ, alt reads and the flag explanations in the tooltip."""
    zyg = g.get("zyg", ".")
    carrier = zyg in CARRIER
    label = ({"het": "heteroplasmy", "hom": "homoplasmy"}.get(zyg, zyg) if on_mito
             else "no call" if zyg == "missing" else zyg)
    bits = []
    if carrier and g.get("vaf") is not None:
        bits.append(percent(g["vaf"]))
    if g.get("dp") is not None:
        bits.append(f"{g['dp']}×")
    if g.get("dn") == "DeNovo":
        bits.append("DRAGEN de novo")
    detail = [name, f"GT {g['gt']}" if g.get("gt") else "",
              f"VAF {g['vaf']:.3f}" if g.get("vaf") is not None else "",
              f"depth {g['dp']}" if g.get("dp") is not None else "", f"GQ {g['gq']}" if g.get("gq") is not None else "",
              f"alt reads {g['alt_reads']}" if g.get("alt_reads") is not None else "",
              f"DRAGEN {g['dn']}" + (f" (DQ {g['dq']:g})" if g.get("dq") is not None else "") if g.get("dn") else "",
              f"confidence {g['conf']}" if g.get("conf") not in (None, "missing") else ""]
    flags = g.get("flags") or []
    title = " · ".join(x for x in detail if x) + "".join(f"\n{fl}: {flag_text(fl, qc)}" for fl in flags)
    chips = "".join(f"<span class='qf'>{esc(fl)}</span>" for fl in flags)
    return (f"<div class='call {esc(g.get('conf', ''))}' title='{esc(title)}'>"
            + (f"<span class='who'>{esc(name)}</span> " if name else "")
            + f"<span class='zyg {esc(zyg)}'>{esc(label)}</span> "
              f"<span class='muted small'>{esc(' · '.join(bits))}</span>"
            + (f"<div class='qfs'>{chips}</div>" if chips else "") + "</div>")


def genotype_head(view: View) -> list[tuple[str, str]]:
    """(header HTML, column width key) of the genotype column(s)."""
    if view.layout == "single":
        return [("genotype", "genotype")]
    if view.layout == "family":
        heads = []
        for s in view.members:
            sex = view.sexes.get(s, "")
            role = next((r for (f, x), r in view.role.items() if x == s), "")
            short = role.replace(" (inferred)", "").replace(" (VCF)", "")
            title = " · ".join(x for x in (s, sex, role) if x)
            symbol = {"male": "♂", "female": "♀"}.get(sex, "")
            heads.append((f"<span title='{esc(title)}'>{esc(s)}<br><span class='small muted'>"
                          f"{symbol} {esc(short)}</span></span>", "member"))
        return heads
    if view.layout == "cohort":
        return [("carriers", "carriers")]
    return []


def genotype_cells(v: dict, view: View) -> str:
    calls = v.get("genotypes") or []
    on_mito = mito(v)
    if view.layout == "single":
        return f"<td>{call_html(calls[0], view.qc, on_mito) if calls else esc(v.get('carriers', '.'))}</td>"
    if view.layout == "family":
        by_sample = {g.get("sample"): g for g in calls}
        return "".join(f"<td>{call_html(by_sample[s], view.qc, on_mito)}</td>" if s in by_sample
                       else "<td><span class='muted'>·</span></td>" for s in view.members)
    if view.layout == "cohort":
        carriers = [g for g in calls if g.get("zyg") in CARRIER]
        shown = carriers[:MAX_CALLS_HTML]
        more = (f"<div class='muted small'>+{len(carriers) - len(shown)} more</div>"
                if len(carriers) > len(shown) else "")
        summary = f"<div class='muted small'>{esc(v.get('carriers', '.'))}</div>"
        return f"<td>{summary}{''.join(call_html(g, view.qc, on_mito, g.get('sample', '')) for g in shown)}{more}</td>"
    return ""


def models_html(v: dict, view: View) -> str:
    """Inheritance models of a P/LP or VUS row (with the sample, unless the file has one)."""
    found = [(g["sample"], m) for g in v.get("genotypes", []) if g.get("zyg") in CARRIER
             for m in view.labels.get((v["file"], g["sample"], v.get("vcf")), [])]
    return "".join(f"<div title='{esc(m)} · {esc(s)}'><span class='model'>{esc(MODEL_LABELS.get(m, m))}</span>"
                   + (f" <span class='muted small'>{esc(s)}</span>" if view.layout != "single" else "") + "</div>"
                   for s, m in found)


def candidate_lead(row: dict, view: View) -> str:
    """Model column of a candidate row: every model with its probands, note and partner variants."""
    groups: dict[tuple, list] = {}
    for c in row["_candidates"]:
        groups.setdefault((c["model"], c["note"], tuple(p for p in c["partners"] if p), c.get("doubtful", "")),
                          []).append(c["sample"])
    out = []
    for (model, note, partners, doubtful), probands in groups.items():
        # A single sample is its own proband; in a family the probands come first on the second line.
        details = (([] if view.layout == "single" else [", ".join(probands)]) + [note]
                   + ([f"with {', '.join(partners)}"] if partners else []))
        out.append(f"<div title='{esc(model)}'><span class='model'>{esc(MODEL_LABELS.get(model, model))}</span>"
                   f"{fit_chip(model, row.get('_gene'))}</div>"
                   f"<div class='muted small' title='{esc(' · '.join(details))}'>{esc(' · '.join(details))}</div>"
                   + (f"<div class='small flag' title='{esc(doubtful)}'>{esc(doubtful)}</div>" if doubtful else ""))
    return "".join(out)


def candidate_table_rows(candidates: list[dict]) -> list[dict]:
    """One row per (file, variant) with the candidates of every proband; doubtful when all of them are."""
    rows: dict[tuple, dict] = {}
    for c in candidates:
        key = (c["file"], c["variant"].get("vcf"))
        rows.setdefault(key, {**c["variant"], "file": c["file"], "_candidates": []})["_candidates"].append(c)
    for r in rows.values():
        r["_doubtful"] = all(c.get("doubtful") for c in r["_candidates"])

    def order(r: dict) -> tuple:
        # Within a model: the gene's mode of inheritance fits (GenCC) before unknown before does not fit.
        return min((inheritance.MODEL_ORDER.index(c["model"]), bool(c.get("doubtful")),
                    {True: 0, None: 1, False: 2}[gene_info.moi_fit(c["model"], r.get("_gene"))],
                    inheritance.CONFIDENCE_ORDER.get(inheritance.call_conf(c), 3), -inheritance.severity(r), r["gene"])
                   for c in r["_candidates"])
    return sorted(rows.values(), key=lambda r: (r["file"], order(r)))


def variant_table(variants: list[dict], table_id: str, kind: str, view: View) -> str:
    """
    Sortable, filterable variant table (kind: plp, vus or candidates). Fixed
    column widths, one line per field (the full text in the tooltip); rows
    marked _doubtful are hidden until the switch above the table is on.
    """
    lead_head = {"plp": "source · model", "vus": "inheritance", "candidates": "model"}[kind]
    heads = (([("file", "file")] if view.show_file else []) + [(lead_head, "lead"), ("gene", "gene"),
             ("change · consequence", None)] + genotype_head(view)
             + [("ClinVar", "ClinVar"), ("condition", None), ("InterVar · ACMG", "InterVar · ACMG")]
             + ([("<span title='PubMind: pathogenicity reported in the literature, extracted by a language "
                  "model; paper-level score 0 (benign) to 1 (pathogenic), evidence confidence 0 to 3'>PubMind</span>",
                  "PubMind")] if view.pubmind else [])
             + [("<span title='gnomAD maximum allele frequency (any population)'>gnomAD</span>", "gnomAD AF"),
                ("REVEL", "REVEL"), ("CADD", "CADD"), ("AlphaMissense", "AlphaMissense")])
    fixed = sum(COLUMN_WIDTHS[w] for _, w in heads if w)
    cols = "".join(f"<col style='width:{COLUMN_WIDTHS[w]}px'>" if w else "<col>" for _, w in heads)

    def num_cell(value) -> str:
        return f"<td class='num'>{esc(value) if value not in (None, '') else '.'}</td>"

    rows = []
    for v in variants[:MAX_ROWS_HTML]:
        stars = "★" * v["stars"] + "☆" * (4 - v["stars"]) if v.get("clnsig") not in (".", "") else ""
        filt = v.get("filter", ".")
        hgvs_text = " ".join(x for x in (v.get("hgvs_c"), v.get("hgvs_p")) if x not in (None, ".", ""))
        main = hgvs_text or v.get("vcf", "")
        where = [consequence(v)] + [x for x in (v.get("transcript"),) if x not in (None, ".", "")] \
            + ([v.get("vcf", "")] if hgvs_text else [])
        filter_chip = f"<span class='flag'>{esc(filt)}</span> " if filt not in ("PASS", ".") else ""
        acmg = v.get("acmg", ".")
        if kind == "candidates":
            lead = candidate_lead(v, view)
        elif kind == "plp":
            lead = f"<span class='src'>{esc(v['source'])}</span>" + models_html(v, view)
        else:
            lead = models_html(v, view) or "<span class='muted'>·</span>"
        tx_title = ":".join(x for x in (v.get("transcript"), hgvs_text) if x not in (None, ".", ""))
        classes = (f"conf-{worst_confidence(v)}" + (" doubtful" if v.get("_doubtful") else "")
                   + (" offpanel" if view.panels and not v.get("_panel") else ""))
        rows.append(
            f"<tr class='{classes}'>"
            + (f"<td title='{esc(v['file'])}'>{esc(v['file'])}</td>" if view.show_file else "")
            + f"<td>{lead}</td>"
            + gene_cell(v, view)
            + f"<td><div class='mono' title='{esc(tx_title or main)}'>{esc(main)}</div>"
              f"<div class='muted small' title='{esc(' · '.join(where + [f'FILTER {filt}']))}'>{filter_chip}"
              f"{esc(' · '.join(where))}</div></td>"
            + genotype_cells(v, view)
            + f"<td>{sig_badge(v.get('clnsig', '.'))}<div class='stars' title='{esc(v.get('clnrevstat', ''))}'>"
              f"{stars}</div></td>"
            + condition_cell(v)
            + f"<td>{sig_badge(v.get('intervar', '.'))}<div class='muted small' title='{esc(acmg)}'>"
              f"{esc(acmg.replace(',', ' ') if acmg != '.' else '')}</div></td>"
            + (pubmind_cell(v) if view.pubmind else "")
            + num_cell(v.get("gnomad_af", ".")) + num_cell(v.get("revel", ".")) + num_cell(v.get("cadd", "."))
            + f"<td title='{esc(v.get('alphamissense', '.'))}'>{esc(v.get('alphamissense', '.'))}</td></tr>")
    more = (f'<p class="muted">Showing {MAX_ROWS_HTML:,} of {len(variants):,}; all are in the TSV file.</p>'
            if len(variants) > MAX_ROWS_HTML else "")
    n_doubtful = sum(1 for v in variants[:MAX_ROWS_HTML] if v.get("_doubtful"))
    toggle = (f'<label class="toggle"><input type="checkbox" class="doubtful-toggle" data-table="{table_id}"> '
              f'show {n_doubtful:,} doubtful (low-confidence calls, Mendelian inconsistencies)</label>'
              if n_doubtful else "")
    n_off = sum(1 for v in variants[:MAX_ROWS_HTML] if view.panels and not v.get("_panel"))
    toggle += (f'<label class="toggle"><input type="checkbox" class="offpanel-toggle" data-table="{table_id}"> '
               f'show {n_off:,} outside the gene panel</label>' if n_off else "")
    return (f'<div class="toolbar"><input class="filter" data-table="{table_id}" type="search" '
            f'placeholder="Filter by gene, sample, condition, flag ..." aria-label="Filter variants">{toggle}</div>'
            f'<div class="wrap"><table id="{table_id}" class="sortable variants hide-doubtful hide-offpanel" '
            f'style="min-width:{fixed + 2 * FLEX_MIN}px"><colgroup>{cols}</colgroup>'
            f'{thead([h for h, _ in heads], [h for h, w in heads if w in NUMERIC_COLUMNS])}'
            f'<tbody>{"".join(rows)}</tbody></table></div>{more}')


def quality_legend(qc: dict) -> str:
    flags = " · ".join(f"<b>{esc(k if k != 'FT' else 'FT:…')}</b> {esc(flag_text(k, qc))}" for k in FLAG_TEXT)
    return (f"<details class='legend-box'><summary>Genotype quality: <span class='call high'>high</span> "
            f"<span class='call caution'>caution: one flag</span> <span class='call low'>low: filtered by the "
            f"caller (FILTER / FT) or two or more flags</span></summary><p class='small'>{flags}. Low-confidence "
            f"calls need confirmation (e.g. Sanger or IGV) before they are reported. Hover a genotype for its "
            f"GQ, alt reads and flags. On chrM the variant allele fraction is the heteroplasmy level.</p></details>")


def plp_section(plp: list[dict], totals: list[dict], view: View) -> str:
    not_carried = sum(t["clinvar_plp_not_carried"] for t in totals if isinstance(t["clinvar_plp_not_carried"], int))
    note = (f'<p class="note">{not_carried:,} ClinVar P/LP record(s) with a reference genotype (0/0 or ./., e.g. '
            f'DRAGEN targeted-caller calls at known pathogenic sites) are not carried and not listed.</p>'
            if not_carried else "")
    if not plp:
        return note + ('<p class="muted">No carried pathogenic / likely pathogenic variants '
                       '(or neither ClinVar nor InterVar is in the protocol).</p>')
    low = sum(1 for v in plp if worst_confidence(v) != "high")
    low_note = (f'<p class="note">{low:,} of {len(plp):,} have a genotype call of reduced confidence '
                f'(marked in the genotype column).</p>' if low else "")
    return (note + low_note + '<p class="muted">ClinVar pathogenic / likely pathogenic, plus InterVar (automated '
            'ACMG) P/LP calls that ClinVar does not classify as P/LP. All are listed, whatever the call quality. '
            'Review every call before reporting.</p>'
            + quality_legend(view.qc) + variant_table(plp, "plp", "plp", view))


def vus_rows(vus: list[dict]) -> list[dict]:
    """VUS whose every carrier call has low confidence are doubtful."""
    return [{**v, "_doubtful": best_confidence(v) == "low"} for v in vus]


def vus_section(rows: list[dict], view: View) -> str:
    if not rows:
        return '<p class="muted">No carried rare coding/splicing ClinVar VUS (or ClinVar not in the protocol).</p>'
    return ('<p class="muted">ClinVar uncertain significance and conflicting classifications that the sample '
            f'carries, exonic or splicing, with a gnomAD frequency below {RARE_AF:.0%} in every population '
            '(or absent). Sorted by REVEL.</p>'
            + variant_table(rows, "vus", "vus", view))


def candidates_section(rows: list[dict], files: list[str], meta: dict, view: View, family_max: int) -> str:
    family_files = [f for f in files if 0 < len(meta[f]["samples"]) <= family_max]
    if not family_files:
        return '<p class="muted">Inheritance analysis runs for single samples and families (few samples per file).</p>'
    intro = (f'<p class="muted">Rare (gnomAD max AF below {inheritance.RECESSIVE_AF:.0%}; de novo and dominant below '
             f'{inheritance.DOMINANT_AF:.1%}) protein-altering or splicing variants of the affected sample(s) '
             'that fit an inheritance model and segregate with the unaffected relatives. The pedigree comes '
             'from --ped, else from the VCF header (##PEDIGREE), else from the relatedness of the samples; '
             'without affected status the children are the probands, and a single sample is its own proband '
             '(homozygous, hemizygous, possible compound heterozygous). Without parents, homozygous and '
             'unphased compound heterozygous candidates are damaging variants only (loss of function, splicing, '
             'REVEL ≥ 0.5 or CADD ≥ 20). One row per variant; confident calls first, then loss of function, '
             'then REVEL ≥ 0.5 / CADD ≥ 20.</p>')
    if not rows:
        return intro + '<p class="muted">No candidates (or no affected sample in the PED file).</p>'
    counts = Counter(m for r in rows if not r["_doubtful"] for m in dict.fromkeys(c["model"] for c in r["_candidates"]))
    summary = " · ".join(f"{m}: {counts[m]}" for m in inheritance.MODEL_ORDER if counts.get(m))
    doubtful = sum(1 for r in rows if r["_doubtful"])
    summary += f" · doubtful (hidden): {doubtful}" if doubtful else ""
    return intro + f'<p>{esc(summary)}</p>' + variant_table(rows, "candidates", "candidates", view)


def findings_table(findings: dict, role: dict, show_file: bool) -> str:
    """What every sample carries: one row per sample, the counts of the file cards."""
    show_cand = any(x["inheritance_candidates"] != "NA" for x in findings.values())
    labels = {"coding_rare": "rare coding", "clinvar_plp": "ClinVar P/LP", "intervar_plp": "InterVar-only P/LP",
              "vus_candidates": "rare coding VUS", "inheritance_candidates": "inheritance candidates"}
    keys = [k for k in FINDING_COLUMNS if show_cand or k != "inheritance_candidates"]
    body = []
    for (f, s), x in findings.items():
        models = " · ".join(f"{m} {n}" for m, n in x["models"].items())
        cells = []
        for k in keys:
            cls = "num flag" if k == "clinvar_plp" and isinstance(x[k], int) and x[k] > 0 else "num"
            tip = f" title='{esc(models)}'" if k == "inheritance_candidates" and models else ""
            cells.append(f"<td class='{cls}'{tip}>{fmt(x[k])}</td>")
        body.append("<tr>" + (f"<td>{esc(f)}</td>" if show_file else "") + f"<td class='label'>{esc(s)}</td>"
                    + f"<td>{esc(role.get((f, s), ''))}</td>" + "".join(cells) + "</tr>")
    head = (["file"] if show_file else []) + ["sample", "role"] + [labels[k] for k in keys]
    return ('<h3>Findings per sample</h3><p class="muted small">Variant records each sample carries (het, hom or '
            'hemizygous genotype), counted like the totals of the file: rare coding (protein-altering or splicing, '
            f'gnomAD max AF below {RARE_AF:.0%}), ClinVar P/LP, InterVar P/LP not P/LP in ClinVar, rare coding '
            'ClinVar VUS / conflicting'
            + (', and the inheritance candidates of the sample as proband (doubtful ones not counted; the models '
               'on hover)' if show_cand else '') + '.</p>'
            f'<div class="wrap"><table>{thead(head, [labels[k] for k in keys])}'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def samples_section(files: list[str], sample_qc: dict, role: dict, show_file: bool, layout: str,
                    findings: dict | None = None) -> str:
    """Per-sample QC (sex check, depth, het/hom, Ti/Tv), what every sample carries (more than one sample),
    pairwise relatedness and trio Mendelian errors. A single sample shows no role and no counts that only
    make sense between samples; its findings are the file card."""
    rows = sample_rows(files, sample_qc)
    single = layout == "single"
    carried = findings_table(findings, role, show_file) if findings and not single else ""
    if not rows:
        return '<p class="muted">No sample QC (skipped with --skip_sample_qc, or a sites-only VCF).</p>' + carried

    def check(value: str | None) -> str:
        return {"ok": "<span class='ok'>✔</span>", "mismatch": "<span class='flag'>✘ mismatch</span>"}.get(
            value or "", "")

    show_ped = not single or any(r["sex_declared"] not in (".", "", None) for r in rows)
    counts = ("mean_depth", "snv_het", "snv_hom", "het_hom_ratio", "ts_tv", "indels") \
        + (() if single else ("singletons", "missing"))
    body = []
    for r in rows:
        ratio = r.get("x_het_hom_ratio")
        body.append(
            "<tr>" + (f"<td>{esc(r['file'])}</td>" if show_file else "")
            + f"<td class='label'>{esc(r['sample'])}</td>"
            + ("" if single else f"<td>{esc(role.get((r['file'], r['sample']), ''))}</td>")
            + (f"<td>{esc(r['sex_declared'])}</td>" if show_ped else "")
            + f"<td>{esc(r.get('sex_inferred', '.'))} "
              f"<span class='muted small'>X het/hom {esc(ratio if ratio is not None else '.')}</span> "
              f"{check(r.get('sex_check'))}</td>"
            + "".join(f"<td class='num'>{fmt(r.get(k)) if r.get(k) is not None else '.'}</td>" for k in counts)
            + "</tr>")
    numeric = ["mean depth", "SNV het", "SNV hom", "het/hom", "Ti/Tv", "indels"] \
        + ([] if single else ["singletons", "missing"])
    head = ((["file"] if show_file else []) + ["sample"] + ([] if single else ["role"])
            + (["sex (PED)"] if show_ped else []) + ["sex from chrX"] + numeric)
    sources = [(f, sample_qc.get(f, {}).get("pedigree_source")) for f in files]
    notes = "".join(f'<p class="muted small">{esc(f + ": ") if show_file else ""}Pedigree '
                    f'{esc(PEDIGREE_SOURCES.get(src, src))}.</p>' for f, src in sources if src)
    out = [notes + f'<div class="wrap"><table>{thead(head, numeric)}<tbody>{"".join(body)}</tbody></table></div>',
           carried]

    pairs = [(f, p) for f in files for p in sample_qc.get(f, {}).get("pairs", [])]
    if pairs:
        prow = "".join(
            "<tr>" + (f"<td>{esc(f)}</td>" if show_file else "")
            + f"<td>{esc(p['a'])} – {esc(p['b'])}</td><td class='num'>{fmt(p['sites'])}</td>"
            + f"<td class='num'>{esc(p['kinship'] if p['kinship'] is not None else '.')}</td>"
            + f"<td class='num'>{esc(p['ibs0_rate'] if p['ibs0_rate'] is not None else '.')}</td>"
            + f"<td>{esc(p['relationship'])}</td><td>{esc(p.get('expected') or '.')}</td>"
            + f"<td>{check(p.get('check'))}</td></tr>" for f, p in pairs)
        phead = (["file"] if show_file else []) + ["pair", "sites", "kinship", "IBS0 rate", "from the data",
                                                   "PED", "check"]
        out.append('<h3>Relatedness</h3><p class="muted small">KING-robust kinship and opposite-homozygote (IBS0) '
                   'rate on autosomal PASS SNVs with GQ ≥ 20 and depth ≥ 10: about 0.5 duplicate, 0.25 '
                   'first-degree (IBS0 near 0: parent–child, otherwise full siblings), 0.125 second-degree, '
                   '0 unrelated.</p>'
                   f'<div class="wrap"><table>{thead(phead, ("sites", "kinship", "IBS0 rate"))}'
                   f'<tbody>{prow}</tbody></table></div>')
    trios = [(f, t) for f in files for t in sample_qc.get(f, {}).get("trios", [])]
    if trios:
        trow = ""
        for f, t in trios:
            rate = f"{t['error_rate']:.2%}" if t.get("error_rate") is not None else "."
            trow += ("<tr>" + (f"<td>{esc(f)}</td>" if show_file else "")
                     + f"<td>{esc(t['child'])}</td><td>{esc(t['father'])}</td><td>{esc(t['mother'])}</td>"
                     + f"<td class='num'>{fmt(t['sites'])}</td><td class='num'>{fmt(t['errors'])}</td>"
                     + f"<td class='num'>{esc(rate)}</td></tr>")
        trio_head = (["file"] if show_file else []) + ["child", "father", "mother", "sites", "Mendelian errors", "rate"]
        out.append('<h3>Trios</h3><p class="muted small">Genotypes of the child that neither parent can explain, '
                   'on the same high-quality SNVs; a few tenths of a percent is usual, several percent point to a '
                   'sample swap or a wrong pedigree.</p>'
                   f'<div class="wrap"><table>{thead(trio_head, ("sites", "Mendelian errors", "rate"))}'
                   f'<tbody>{trow}</tbody></table></div>')
    return "".join(out)


def human_size(n: int) -> str:
    if n < 0:
        return "missing"
    for unit, div in (("TB", 1024 ** 4), ("GB", 1024 ** 3), ("MB", 1024 ** 2), ("KB", 1024)):
        if n >= div:
            return f"{n / div:.1f} {unit}"
    return f"{n} B"


def key_versions(databases: list[dict], annovar_version: str) -> str:
    """Header line, e.g. 'ClinVar 2026-10-25 · gnomAD 4.1 · dbNSFP 5.4a · ANNOVAR 2025-03-02'."""
    labels: list[str] = []
    for key in KEY_SOURCES:
        for d in databases:
            source, version = describe_database(d["database"])
            label = f"{key} {version}"
            if source.split(" ")[0] == key and version and label not in labels:
                labels.append(label)
    if annovar_version:
        labels.append(f"ANNOVAR {annovar_version}")
    return " · ".join(labels)


def sources_section(databases: list[dict], genomes: list[str], annovar_version: str) -> str:
    """Table of the ANNOVAR databases used for the given genome builds."""
    rows = [d for d in databases if d.get("genome") in genomes]
    if not rows:
        return '<p class="muted">The list of databases used is not available.</p>'
    show_build = len({d["genome"] for d in rows}) > 1
    body = []
    for d in rows:
        source, version = describe_database(d["database"])
        alias = f' <span class="muted small">({esc(d["protocol"])})</span>' if d["protocol"] != d["database"] else ""
        size = int(d.get("bytes") or -1)
        date = d.get("file_date", "")
        cells = ([f'<td><span class="badge">{esc(d["genome"])}</span></td>'] if show_build else []) + [
            f'<td class="mono">{esc(d["database"])}{alias}</td>', f"<td>{esc(source)}</td>",
            f"<td>{esc(version) or '–'}</td>", f'<td>{esc(OPERATIONS.get(d["operation"], d["operation"]))}</td>',
            f'<td class="num{" flag" if date == "missing" else ""}">{esc(date)}</td>',
            f'<td class="num">{human_size(size)}</td>']
        body.append("<tr>" + "".join(cells) + "</tr>")
    head = (["build"] if show_build else []) + ["database", "source", "version", "type", "file date", "size"]
    note = ("File date: the ANNOVAR release date for files downloaded by annotate_variation.pl, otherwise the "
            "date the file was built or installed.")
    return (f'<p class="muted small">ANNOVAR {esc(annovar_version or "unknown version")} · {note}</p>'
            f'<div class="wrap"><table>{thead(head, ("file date", "size"))}'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def panel_note(panels: list) -> str:
    if not panels:
        return ""
    names = "; ".join(f"<b>{esc(p.name)}</b> ({len(p.genes):,} genes)" for p in panels)
    return (f'<p class="note">Gene panel{"s" if len(panels) > 1 else ""}: {names}. Variants in other genes are '
            'hidden in the tables below (switch above each table); the TSV and JSON files list every variant.</p>')


def gene_sources(genes: gene_info.GeneTable) -> str:
    """Where the gene information comes from (the header lines of the gene table)."""
    if not genes.genes:
        return ""
    return ('<p class="muted small">Gene information (modes of inheritance, diseases, LOEUF): '
            + " · ".join(esc(s) for s in genes.sources[1:] or genes.sources) + '.</p>')


CSS = """
:root { --bg:#f7f8fa; --panel:#fff; --fg:#1d2433; --muted:#667085; --line:#e4e7ec; --head:#1d2433;
        --accent:#2563eb; --warn:#b42318; --warnbg:#fef3f2; --lp:#b54708; --lpbg:#fffaeb; --vus:#175cd3; --vusbg:#eff8ff;
        --c1:#2563eb; --c2:#7c3aed; --c3:#0891b2; --c4:#16a34a; --c5:#94a3b8; --c6:#f59e0b; --c7:#cbd5e1; --c8:#e2e8f0; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#0f1420; --panel:#171e2e; --fg:#e6e9ef; --muted:#98a2b3; --line:#2a3347; --head:#0b101a;
          --warnbg:#3a1714; --warn:#fda29b; --lpbg:#3a2a10; --lp:#fec84b; --vusbg:#102a43; --vus:#84caff;
          --c5:#64748b; --c7:#475569; --c8:#334155; }
}
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,Arial,sans-serif; }
header { background:var(--head); color:#fff; padding:22px 28px; }
header h1 { margin:0; font-size:21px; font-weight:650; }
header p { margin:4px 0 0; color:#c7cdd8; font-size:13px; }
main { max-width:1760px; margin:0 auto; padding:20px 28px 48px; }
h2 { font-size:16px; margin:30px 0 10px; }
.muted { color:var(--muted); }
.overview { display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px; }
.kpi { background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:14px 16px; }
.kpi b { display:block; font-size:24px; font-variant-numeric:tabular-nums; }
.kpi span { color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.04em; }
.cards { display:grid; grid-template-columns:repeat(auto-fill,minmax(300px,1fr)); gap:12px; }
.card { background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:14px 16px; }
.card-head { display:flex; justify-content:space-between; gap:8px; align-items:center; }
.card .name { font-weight:650; word-break:break-all; }
.badge { background:var(--accent); color:#fff; border-radius:999px; padding:1px 9px; font-size:12px; }
.big { font-size:24px; font-weight:650; margin:6px 0 8px; font-variant-numeric:tabular-nums; }
.big span { font-size:13px; font-weight:400; color:var(--muted); }
dl { display:grid; grid-template-columns:auto 1fr; gap:3px 12px; margin:0; font-size:13px; }
dt { color:var(--muted); } dd { margin:0; text-align:right; font-variant-numeric:tabular-nums; }
dd.warn { color:var(--warn); font-weight:650; }
.panel { background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:16px; }
.legend { display:flex; flex-wrap:wrap; gap:6px 16px; margin-bottom:12px; font-size:12px; color:var(--muted); }
.legend i { display:inline-block; width:10px; height:10px; border-radius:2px; margin-right:6px; vertical-align:-1px; }
.bar-row { display:grid; grid-template-columns:minmax(120px,240px) 1fr; gap:12px; align-items:center; margin:6px 0; }
.bar-label { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; font-size:13px; }
.stack { display:flex; height:18px; border-radius:4px; overflow:hidden; background:var(--line); }
.stack span { display:block; height:100%; }
table { border-collapse:collapse; width:100%; font-size:13px; }
th, td { padding:5px 10px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }
th { position:sticky; top:0; background:var(--panel); color:var(--muted); font-weight:600; white-space:nowrap; cursor:default; }
table.sortable th { cursor:pointer; } table.sortable th:hover { color:var(--fg); }
td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap; }
td.label { font-weight:600; white-space:nowrap; }
td.mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:12px; }
td.gene { font-weight:650; }
.stars { color:#f5b301; white-space:nowrap; letter-spacing:1px; }
.sig { border-radius:4px; padding:1px 6px; white-space:nowrap; font-size:12px; }
.sig.path { background:var(--warnbg); color:var(--warn); } .sig.lpath { background:var(--lpbg); color:var(--lp); }
.sig.vus { background:var(--vusbg); color:var(--vus); } .sig.other { background:var(--line); color:var(--muted); }
.src { font-size:12px; font-weight:600; color:var(--muted); }
.zyg { border-radius:4px; padding:1px 6px; font-size:12px; font-weight:600; background:var(--line); }
.zyg.hom, .zyg.hemi { background:var(--accent); color:#fff; }
.flag { color:var(--warn); font-weight:600; }
.small { font-size:11px; }
.note { background:var(--lpbg); color:var(--lp); border-radius:8px; padding:8px 12px; font-size:13px; }
/* variant tables: fixed column widths (colgroup), one line per field, full text in the tooltip */
table.variants { table-layout:fixed; width:100%; }
table.variants td, table.variants th { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; padding:5px 6px; }
table.variants th { font-size:12px; }
table.variants td > div, table.variants .call { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
table.variants .gene { font-weight:650; }
table.variants .mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:12px; }
table.variants .qfs { white-space:normal; line-height:1.3; }
table.hide-doubtful tr.doubtful { display:none; }
tr.doubtful > td { opacity:.7; }
.clamp { white-space:normal !important; display:-webkit-box; -webkit-line-clamp:2; -webkit-box-orient:vertical; overflow:hidden; }
.links a { color:var(--accent); text-decoration:none; margin-right:6px; }
.toolbar { display:flex; flex-wrap:wrap; align-items:center; gap:8px 18px; margin-bottom:10px; }
input.filter { width:100%; max-width:420px; padding:7px 10px; border:1px solid var(--line);
               border-radius:8px; background:var(--panel); color:var(--fg); font:inherit; }
label.toggle { font-size:13px; color:var(--muted); cursor:pointer; user-select:none; }
.conseq td { min-width:140px; }
.hbar { position:relative; height:18px; background:var(--line); border-radius:3px; }
.hbar span { position:absolute; inset:0 auto 0 0; background:var(--c1); border-radius:3px; opacity:.85; }
.hbar em { position:relative; font-style:normal; font-size:12px; padding-left:6px; line-height:18px; mix-blend-mode:normal; }
.wrap { overflow:auto; max-height:70vh; border:1px solid var(--line); border-radius:10px; background:var(--panel);
        scrollbar-gutter:stable; }
details { margin-top:14px; } summary { cursor:pointer; font-weight:600; margin-bottom:10px; }
tbody tr:hover { background:color-mix(in srgb, var(--accent) 7%, transparent); }
h3 { font-size:14px; margin:18px 0 8px; }
.call { border-left:3px solid var(--line); padding:0 0 0 6px; margin:1px 0; white-space:nowrap; }
.call.caution { border-left-color:var(--lp); } .call.low { border-left-color:var(--warn); }
.call.high { border-left-color:#16a34a; }
summary .call { display:inline-block; font-weight:400; font-size:12px; margin-right:10px; }
.who { font-weight:600; }
.zyg.ref, .zyg.missing { background:transparent; color:var(--muted); font-weight:400; }
.qf { font-size:10.5px; border-radius:3px; padding:0 4px; margin-right:4px; background:var(--lpbg); color:var(--lp); }
.call.low .qf { background:var(--warnbg); color:var(--warn); }
tr.conf-low > td:first-child { box-shadow:inset 3px 0 var(--warn); }
tr.conf-caution > td:first-child { box-shadow:inset 3px 0 var(--lp); }
.model { border-radius:4px; padding:1px 6px; font-size:12px; font-weight:600; background:var(--vusbg); color:var(--vus); }
.ok { color:#16a34a; font-weight:600; }
.legend-box { margin:0 0 10px; font-size:13px; } .legend-box summary { font-weight:400; color:var(--muted); }
table.hide-offpanel tr.offpanel { display:none; }
tr.offpanel > td { opacity:.7; }
button.gene { all:unset; cursor:pointer; font-weight:650; display:block; overflow:hidden; text-overflow:ellipsis;
              white-space:nowrap; max-width:100%; }
button.gene:hover, button.gene:focus-visible { color:var(--accent); text-decoration:underline dotted; }
.gene-info { white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
.moi { border-radius:3px; padding:0 4px; margin-right:3px; font-size:10.5px; font-weight:600;
       background:var(--vusbg); color:var(--vus); }
.loeuf { color:var(--muted); } .loeuf.lof { color:var(--warn); font-weight:600; }
.fit { border-radius:3px; padding:0 4px; font-size:10.5px; font-weight:600; }
.fit.ok { background:color-mix(in srgb, #16a34a 15%, transparent); color:#16a34a; }
.fit.no { background:var(--warnbg); color:var(--warn); text-decoration:line-through; }
.linkmenu { position:absolute; z-index:20; display:flex; flex-direction:column; min-width:120px; padding:4px 0;
            background:var(--panel); border:1px solid var(--line); border-radius:8px;
            box-shadow:0 6px 20px rgba(0,0,0,.18); font-size:13px; }
.linkmenu[hidden] { display:none; }
.linkmenu a { padding:5px 12px; color:var(--fg); text-decoration:none; }
.linkmenu a:hover { background:color-mix(in srgb, var(--accent) 10%, transparent); color:var(--accent); }
"""

JS = """
document.querySelectorAll('input.filter').forEach(filter => filter.addEventListener('input', () => {
  const q = filter.value.toLowerCase();
  document.querySelectorAll('#' + filter.dataset.table + ' tbody tr').forEach(r => { r.style.display = r.textContent.toLowerCase().includes(q) ? '' : 'none'; });
}));
document.querySelectorAll('input.doubtful-toggle').forEach(box => box.addEventListener('change', () => {
  document.getElementById(box.dataset.table).classList.toggle('hide-doubtful', !box.checked);
}));
document.querySelectorAll('input.offpanel-toggle').forEach(box => box.addEventListener('change', () => {
  document.getElementById(box.dataset.table).classList.toggle('hide-offpanel', !box.checked);
}));
// Link menu of a variant and its gene (ClinVar, gnomAD, Franklin, VarSome, OMIM, GenCC, IGV),
// opened from the gene name. IGV: the locus in a running IGV desktop (port 60151).
const menu = document.createElement('div');
menu.className = 'linkmenu';
menu.hidden = true;
document.body.appendChild(menu);
document.addEventListener('click', e => {
  const button = e.target.closest && e.target.closest('button.gene');
  if (!button) {
    if (!menu.contains(e.target)) menu.hidden = true;
    return;
  }
  menu.replaceChildren(...JSON.parse(button.dataset.links).map(([label, url]) => {
    const a = document.createElement('a');
    a.textContent = label;
    a.href = url;
    a.target = '_blank';
    a.rel = 'noopener';
    if (label === 'IGV') a.addEventListener('click', ev => {
      ev.preventDefault();
      fetch(url, {mode: 'no-cors'}).catch(() => window.open(url, '_blank'));
    });
    return a;
  }));
  const r = button.getBoundingClientRect();
  menu.style.left = `${r.left + window.scrollX}px`;
  menu.style.top = `${r.bottom + window.scrollY + 2}px`;
  menu.hidden = false;
});
menu.addEventListener('click', e => { if (e.target.closest('a')) menu.hidden = true; });
document.addEventListener('keydown', e => { if (e.key === 'Escape') menu.hidden = true; });
document.querySelectorAll('table.sortable').forEach(t => t.querySelectorAll('th').forEach((th, i) => {
  th.addEventListener('click', () => {
    const asc = th.dataset.asc !== '1'; th.dataset.asc = asc ? '1' : '0';
    const rows = [...t.tBodies[0].rows];
    const val = r => { const s = r.cells[i].textContent.trim(); const n = parseFloat(s.replace(/,/g, '')); return isNaN(n) ? s.toLowerCase() : n; };
    rows.sort((a, b) => { const x = val(a), y = val(b); return (x > y ? 1 : x < y ? -1 : 0) * (asc ? 1 : -1); });
    rows.forEach(r => t.tBodies[0].appendChild(r));
  });
}));
// Tooltips: a cell's full text when it is cut on screen, or the extra detail a
// title carries (review status, call quality); never a copy of what is visible.
const flat = s => s.replace(/\\s+/g, ' ').trim();
document.addEventListener('mouseover', e => {
  const el = e.target.closest && e.target.closest('table.variants [title], table.variants [data-tip]');
  if (!el) return;
  if (el.dataset.tip === undefined) el.dataset.tip = el.getAttribute('title') || '';
  const cut = el.scrollWidth > el.clientWidth + 1 || el.scrollHeight > el.clientHeight + 1;
  if (cut || flat(el.dataset.tip) !== flat(el.textContent)) el.setAttribute('title', el.dataset.tip);
  else el.removeAttribute('title');
});
"""


@dataclass
class Run:
    """Everything the reports are built from."""
    counts: dict
    func: dict
    meta: dict
    plp: list
    vus: list
    rare: dict
    builds: dict
    sample_qc: dict
    databases: list
    version: str = ""
    annovar_version: str = ""
    family_max: int = 10
    candidates: list = field(default_factory=list)
    labels: dict = field(default_factory=dict)
    findings: dict = field(default_factory=dict)     # (file, sample) -> what it carries (sample_findings)
    genes: gene_info.GeneTable = field(default_factory=gene_info.GeneTable)
    panels: list = field(default_factory=list)       # gene_info.Panel per --panel file


def write_html(path, files: list[str], run: Run, totals, totals_header, totals_rows, wide_header, wide_rows,
               plp, vus, candidates):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    func, builds = run.func, run.builds

    def total(key: str) -> int:
        return sum(t[key] for t in totals if isinstance(t[key], int))

    report_genomes = sorted({t["genome"] for t in totals})
    used = [d for d in run.databases if d.get("genome") in report_genomes]
    versions = key_versions(used, run.annovar_version)
    genomes = ", ".join(report_genomes) or "NA"
    view = make_view(files, run)
    n_samples = sum(len(run.meta[f]["samples"]) for f in files)
    low = sum(1 for v in plp if worst_confidence(v) == "low")
    # Candidates and VUS: one table row per variant; doubtful rows are hidden by default.
    cand_rows = candidate_table_rows(candidates)
    vus_table = vus_rows(vus)
    n_cand = sum(1 for r in cand_rows if not r["_doubtful"])
    n_vus = sum(1 for v in vus_table if not v["_doubtful"])

    def in_panel(rows: list[dict]) -> str:
        n = sum(1 for r in rows if r.get("_panel") and not r.get("_doubtful"))
        return f" · {n:,} in the gene panel" if run.panels else ""

    kpis = "".join(f'<div class="kpi"><b>{v}</b><span>{k}</span></div>' for k, v in (
        ("variant records", f"{total('total'):,}"), ("samples", f"{n_samples:,}"),
        ("ClinVar P/LP carried", f"{total('clinvar_plp'):,}"), ("InterVar-only P/LP", f"{total('intervar_plp'):,}"),
        ("P/LP low-confidence call", f"{low:,}"), ("rare coding VUS", f"{n_vus:,}"),
        ("inheritance candidates", f"{n_cand:,}"), ("genome builds", esc(genomes))))
    # A single-file report is named after its file.
    title = f"{files[0]} · annovar-nf annotation report" if len(files) == 1 else "annovar-nf annotation report"
    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<style>{CSS}</style></head><body>
<header><h1>{esc(title)}</h1>
<p>{len(totals)} input file(s) · {n_samples} sample(s) · generated {stamp} · annovar-nf {esc(run.version)} ·
counts are variant records after splitting multi-allelic sites, per input file</p>
{f"<p>{esc(versions)}</p>" if versions else ""}</header>
<main>
<section class="overview">{kpis}</section>
{panel_note(run.panels)}

<h2>Input files</h2>
<section class="cards">{file_cards(totals, builds)}</section>

<h2>Samples</h2>
<section>{samples_section(files, run.sample_qc, view.role, view.show_file, view.layout,
                          {k: v for k, v in run.findings.items() if k[0] in files})}</section>

<h2>Pathogenic / likely pathogenic ({len(plp):,}{in_panel(plp)})</h2>
<section>{plp_section(plp, totals, view)}</section>

<h2>Inheritance candidates ({n_cand:,}{in_panel(cand_rows)})</h2>
<section>{candidates_section(cand_rows, files, run.meta, view, run.family_max)}</section>

<h2>Variants of uncertain significance ({n_vus:,}{in_panel(vus_table)})</h2>
<section>{vus_section(vus_table, view)}</section>

<h2>Genome build</h2>
<section>{build_table(builds, files)}</section>

<h2>Genomic regions</h2>
<section class="panel">{region_chart(func, files)}</section>

<h2>Exonic consequences</h2>
<section class="panel" style="overflow-x:auto">{consequence_chart(func, files)}</section>

<h2>Annotation sources</h2>
<section>{sources_section(used, report_genomes, run.annovar_version)}{gene_sources(run.genes)}</section>

<details><summary>All per-file metrics</summary>
<div class="wrap">{html_table(totals_header, totals_rows, 2)}</div></details>
<details><summary>Func × ExonicFunc matrix</summary>
<div class="wrap">{html_table(wide_header, wide_rows, 2)}</div></details>
</main>
<script>{JS}</script>
</body></html>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(page)


def json_variant(v: dict, labels: dict) -> dict:
    """A listed variant for the JSON report: annotation, genotype calls and inheritance models per sample."""
    out = {k: v.get(k) for k in VARIANT_COLUMNS if k != "carriers"}
    for k in ("source", "kind"):
        if k in v:
            out[k] = v[k]
    out["genotypes"] = v.get("genotypes", [])
    out["inheritance"] = {g["sample"]: labels[(v["file"], g["sample"], v.get("vcf"))]
                          for g in v.get("genotypes", []) if (v["file"], g.get("sample"), v.get("vcf")) in labels}
    return out


def write_json(path: str, file_id: str, run: Run, totals: dict, plp, vus, candidates) -> None:
    """Machine-readable report of one input file (e.g. for a LIS): versions, QC, samples, listed variants."""
    qc = run.sample_qc.get(file_id, {})
    build = run.builds.get(file_id, {})
    genome = run.meta[file_id]["genome"]
    databases = []
    for d in run.databases:
        if d.get("genome") == genome:
            source, version = describe_database(d["database"])
            databases.append({"database": d["database"], "protocol": d["protocol"], "source": source,
                              "version": version, "type": OPERATIONS.get(d["operation"], d["operation"]),
                              "file_date": d.get("file_date")})
    payload = {
        "schema": "annovar-nf.report", "schema_version": 1,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline": {"name": "annovar-nf", "version": run.version},
        "annovar": {"version": run.annovar_version, "databases": databases},
        "input_file": {"id": file_id, "genome": genome, "filename": build.get("filename"),
                       "build_determined_by": build.get("source"), "samples": run.meta[file_id]["samples"]},
        "genotype_quality_thresholds": run.meta[file_id].get("qc") or {},
        "totals": totals,
        "samples": qc.get("samples", []), "relatedness": qc.get("pairs", []), "trios": qc.get("trios", []),
        "sample_findings": [{"sample": s, **run.findings[(f, s)]} for f, s in run.findings if f == file_id],
        "pathogenic": [json_variant(v, run.labels) for v in plp],
        "uncertain_significance": [json_variant(v, run.labels) for v in vus],
        "inheritance_candidates": [{"sample": c["sample"], "model": c["model"], "note": c["note"],
                                    "partners": c["partners"], "doubtful": c.get("doubtful") or None,
                                    "moi_fit": gene_info.moi_fit(c["model"], c["variant"].get("_gene")),
                                    "variant": json_variant(c["variant"], {})}
                                   for c in candidates],
        "gene_information": run.genes.sources,
        "gene_panels": [{"name": p.name, "genes": len(p.genes)} for p in run.panels],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)


def write_report(prefix: str, files: list[str], run: Run) -> tuple[int, int]:
    """Write <prefix>.html, .tsv, _totals / _samples / _plp / _vus / _candidates .tsv (and .json for one file)."""
    wide_header = ["Func.refGene", "ExonicFunc.refGene"] + files
    wide_rows = wide_matrix(run.func, files)
    write_tsv(f"{prefix}.tsv", wide_header, wide_rows)

    totals = totals_table(run.counts, run.func, run.meta, files)
    totals_header = ["file", *TOTAL_COLUMNS]
    totals_rows = [[t[k] for k in totals_header] for t in totals]
    write_tsv(f"{prefix}_totals.tsv", totals_header, totals_rows)
    findings = {k: v for k, v in run.findings.items() if k[0] in files}
    write_tsv(f"{prefix}_samples.tsv", SAMPLE_COLUMNS,
              [[r.get(k) for k in SAMPLE_COLUMNS] for r in sample_table(files, run.sample_qc, findings)])

    plp = [v for v in run.plp if v["file"] in files]
    vus = [v for v in run.vus if v["file"] in files]
    candidates = [c for c in run.candidates if c["file"] in files]
    write_tsv(f"{prefix}_plp.tsv", PLP_COLUMNS,
              [[r.get(k, ".") for k in PLP_COLUMNS] for r in carrier_rows(plp, run.labels)])
    write_tsv(f"{prefix}_vus.tsv", VUS_COLUMNS,
              [[r.get(k, ".") for k in VUS_COLUMNS] for r in carrier_rows(vus, run.labels)])
    write_tsv(f"{prefix}_candidates.tsv", CANDIDATE_COLUMNS,
              [[r.get(k, ".") for k in CANDIDATE_COLUMNS] for r in candidate_rows(candidates)])

    write_html(f"{prefix}.html", files, run, totals, totals_header, totals_rows, wide_header, wide_rows,
               plp, vus, candidates)
    if len(files) == 1:
        write_json(f"{prefix}.json", files[0], run, totals[0], plp, vus, candidates)
    return len(plp), len(vus)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--stats", nargs="+", required=True, help="per-chunk *.stats.json files")
    p.add_argument("--builds", nargs="*", default=[], help="<id>.build.json reports of detect_genome_build.py")
    p.add_argument("--out", required=True,
                   help="<dir>/<name>: per-file reports are <dir>/<id>.<name>.*, the combined one <dir>/<name>.*")
    p.add_argument("--version", default="", help="pipeline version for the report")
    p.add_argument("--databases", help="databases.tsv: ANNOVAR databases used per build (versions)")
    p.add_argument("--annovar-version", default="", help="ANNOVAR release, e.g. 2025-03-02")
    p.add_argument("--sample-qc", nargs="*", default=[], help="<id>.sample_qc.json reports of sample_qc.py")
    p.add_argument("--family-max-samples", type=int, default=10,
                   help="inheritance analysis for files with at most this many samples")
    p.add_argument("--genes", help="gene table (humandb/annovar_nf/genes.tsv, see gene_info.py)")
    p.add_argument("--panel", nargs="*", default=[],
                   help="gene panels: gene lists or PanelApp TSV exports (other genes hidden in the report)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    counts, func, meta, plp, vus, rare = load(args.stats)
    run = Run(counts, func, meta, plp, vus, rare, load_builds(args.builds), load_sample_qc(args.sample_qc),
              load_databases(args.databases), args.version, args.annovar_version, args.family_max_samples)
    files = sorted(counts)
    run.genes = gene_info.load(args.genes)
    run.panels = [gene_info.read_panel(p, run.genes) for p in args.panel]
    annotate_genes([*plp, *vus, *(v for f in rare for v in rare[f])], run.genes, run.panels)
    run.candidates, run.labels = family_analysis(files, meta, rare, run.sample_qc)
    run.findings = sample_findings(files, meta, run.candidates, args.family_max_samples)
    out = Path(args.out)

    for f in files:
        n_plp, n_vus = write_report(str(out.with_name(f"{f}.{out.name}")), [f], run)
        n_cand = sum(1 for c in run.candidates if c["file"] == f)
        print(f"[summarize] {f}: {counts[f].get('total', 0):,} records, {len(meta[f]['samples'])} sample(s), "
              f"{n_plp} P/LP, {n_vus} VUS, {n_cand} inheritance candidates -> {out.with_name(f'{f}.{out.name}')}.*",
              file=sys.stderr)
    if len(files) > 1:
        write_report(str(out), files, run)
        print(f"[summarize] {len(files)} input files -> {out}.*", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
