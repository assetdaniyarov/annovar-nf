#!/usr/bin/env python3
"""
format_annovar_table.py
-----------------------
Post-process an ANNOVAR `*_multianno.txt` table produced with
`--vcfinput --otherinfo` in a single streaming pass:

  1. Replace the generic 'Otherinfo1..N' header columns with the real VCF
     column names (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples>).
     The extra columns ANNOVAR inserts before the VCF line (allele frequency /
     zygosity, QUAL, DP) keep the names Otherinfo1..Otherinfo3.
     Works with both ANNOVAR >= 2019Oct24 (Otherinfo1..N in the header) and
     older releases (a single 'Otherinfo' header column).
  2. Separate genes with ';' in the Gene.* and GeneDetail.* columns.
     table_annovar.pl releases from 2021 on write ',' between the genes of one
     function and ';' between functions; ';' throughout keeps tables from every
     ANNOVAR release alike. With --annovar-vcf the same is applied to the
     annotated VCF (INFO values Gene.*= / GeneDetail.*=, ';' written as \\x3b).
  3. Write the header line and the data rows to separate outputs so that the
     per-chunk results can be concatenated cheaply (header once, bodies as
     concatenated BGZF members).
  4. Collect summary statistics (functional classes, SNV/indel counts, Ti/Tv,
     FILTER=PASS, REVEL, CADD) and the clinically relevant variants carried by
     the sample(s) - ClinVar / InterVar pathogenic and likely pathogenic, rare
     coding ClinVar VUS / conflicting - with gnomAD frequency, in-silico scores
     and the genotype calls: zygosity, allele fraction, depth, GQ and quality
     flags (LOW_DP, LOW_GQ, LOW_ALT_READS, VAF_SKEW, MOSAIC, HOMOLOGY,
     ALT_READS_IN_REF, DRAGEN FORMAT/FT) with a confidence (high / caution /
     low). Files with at most --family-max-samples samples also list every rare
     coding/splicing variant a sample carries, for the inheritance analysis.
     Per sample: the rare coding, P/LP and VUS records it carries.

Only the first columns of each row are split, so very wide multi-sample rows
(hundreds of genotype columns) are cheap to process; genotypes are parsed for
listed variants only (rare coding variants: GT only).

Usage:
    format_annovar_table.py --txt chunk.hg38_multianno.txt --vcf chunk.input.vcf \
        --header-out chunk.header.tsv --stats chunk.stats.json \
        --sample S1 --chunk chr1_1 > chunk.body.tsv
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from typing import IO

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
NUCLEOTIDES = set("ACGT")
# Gene.refGene, GeneDetail.refGene, ... (Gene_refGene with --dot2underline).
GENE_COLUMN_RE = re.compile(r"^(Gene|GeneDetail)[._]")
VCF_GENE_INFO_RE = re.compile(r"((?:^|;)(?:Gene|GeneDetail)[._][^=;]*=)([^;]*)")


def open_text(path: str, mode: str = "rt") -> IO[str]:
    """Open plain or gzip/BGZF-compressed text ('-' means STDIN/STDOUT)."""
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if "r" in mode:
        with open(path, "rb") as fh:
            if fh.read(2) == b"\x1f\x8b":
                return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def read_vcf_columns(vcf_path: str) -> list[str]:
    """Return the '#CHROM ...' column names of a (possibly compressed) VCF."""
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\r\n").split("\t")
            if not line.startswith("#"):
                break
    raise ValueError(f"no '#CHROM' header line found in {vcf_path}")


def build_header(txt_header: list[str], vcf_cols: list[str], n_data_cols: int) -> list[str]:
    """
    Merge the ANNOVAR annotation columns with the real VCF column names.

    n_data_cols is the width of the data rows; it can be larger than the header
    width for ANNOVAR releases that print a single 'Otherinfo' header column.
    """
    first_other = next((i for i, c in enumerate(txt_header) if c.startswith("Otherinfo")), None)
    if first_other is None:
        raise ValueError("no 'Otherinfo' column in the ANNOVAR header; was --otherinfo/--vcfinput used?")
    anno_cols = txt_header[:first_other]
    n_other = n_data_cols - len(anno_cols)
    n_extra = n_other - len(vcf_cols)
    if n_extra < 0:
        raise ValueError(
            f"data rows have {n_other} Otherinfo columns but the VCF has {len(vcf_cols)} columns; "
            f"the ANNOVAR table does not match the VCF")
    return anno_cols + [f"Otherinfo{i}" for i in range(1, n_extra + 1)] + vcf_cols


def gene_columns(header: list[str]) -> list[int]:
    """Indexes of the Gene.* / GeneDetail.* columns."""
    return [i for i, name in enumerate(header) if GENE_COLUMN_RE.match(name)]


def gene_separator(line: str, cols: list[int]) -> str:
    """Use ';' between genes in the given columns of a table row."""
    if not cols:
        return line
    f = line.split("\t", cols[-1] + 1)
    changed = False
    for i in cols:
        if i < len(f) and "," in f[i]:
            f[i] = f[i].replace(",", ";")
            changed = True
    return "\t".join(f) if changed else line


def vcf_gene_separator(line: str) -> str:
    """Use ';' (written as \\x3b) between genes in the Gene.* / GeneDetail.* INFO values."""
    if line.startswith("#"):
        return line
    f = line.split("\t", 8)
    if len(f) < 8 or "Gene" not in f[7]:
        return line
    f[7] = VCF_GENE_INFO_RE.sub(lambda m: m.group(1) + m.group(2).replace(",", r"\x3b"), f[7])
    return "\t".join(f)


def find_column(header: list[str], exact: tuple[str, ...] = (), pattern: str | None = None) -> int | None:
    """Index of the first column matching one of *exact* names, else *pattern*."""
    for name in exact:
        if name in header:
            return header.index(name)
    if pattern:
        rx = re.compile(pattern, re.I)
        for i, col in enumerate(header):
            if rx.search(col):
                return i
    return None


def is_plp(clnsig: str) -> bool:
    """ClinVar significance counts as pathogenic / likely pathogenic."""
    v = clnsig.lower()
    return "pathogenic" in v and "conflicting" not in v


def vus_kind(clnsig: str) -> str | None:
    """'VUS' (Uncertain_significance), 'conflicting' (conflicting classifications of pathogenicity) or None."""
    v = clnsig.lower()
    if v.startswith("conflicting") and "pathogenicity" in v:
        return "conflicting"
    if v.startswith("uncertain_significance"):
        return "VUS"
    return None


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


# ExonicFunc values that change the protein (ANNOVAR refGene vocabulary).
PROTEIN_ALTERING = {
    "nonsynonymous SNV", "stopgain", "stoploss", "startloss", "frameshift insertion",
    "frameshift deletion", "frameshift substitution", "frameshift block substitution",
    "nonframeshift insertion", "nonframeshift deletion", "nonframeshift substitution",
}
CODING_FUNC = {"exonic", "splicing", "exonic;splicing"}
# Mitochondrial contig names (ANNOVAR Chr column / VCF CHROM).
MITO = {"chrM", "chrMT", "M", "MT"}
RARE_AF = 0.01
MAX_LISTED_PER_CHUNK = 2000
INTERVAR_PLP = {"pathogenic", "likely pathogenic"}
ACMG_CRITERIA = ("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2",
                 "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5",
                 "BP6", "BP7")
# gnomAD allele frequencies: overall and population maximum, exome and genome.
# gnomAD 2.1.1 (hg19) names both its exome and genome columns 'AF' / 'AF_popmax'.
GNOMAD_AF_RE = re.compile(r"^(AF|AF_popmax|gnomad\w*?_(exome|genome)_AF(_grpmax)?|gnomAD_(exome|genome)_ALL)$")


def splicing(func: str) -> bool:
    """Splice-site variant of a protein-coding transcript (Func 'splicing' / 'exonic;splicing', not ncRNA_splicing)."""
    return "splicing" in func.split(";")


def shorten(value: str, limit: int) -> str:
    return value if len(value) <= limit else value[: limit - 1] + "…"


def first_transcript(value: str, limit: int = 300) -> str:
    """'GENE:NM_1:exon2:c.A1G:p.M1V,GENE:NM_2:...' -> first entry (a safety limit only:
    the report shortens long text on screen and shows all of it on hover)."""
    return shorten(value.split(",")[0], limit)


@dataclass(frozen=True)
class Qc:
    """Genotype quality thresholds; a call outside them gets a flag."""
    min_dp: int = 10              # LOW_DP: read depth below
    min_gq: int = 10              # LOW_GQ: genotype quality below
    min_alt_reads: int = 4        # LOW_ALT_READS: reads supporting the alt allele (carriers) below
    min_het_vaf: float = 0.2      # VAF_SKEW: heterozygous alt-allele fraction outside [min, max]
    max_het_vaf: float = 0.8
    min_hom_vaf: float = 0.8      # VAF_SKEW: homozygous / hemizygous alt-allele fraction below
    ref_alt_reads: int = 2        # ALT_READS_IN_REF: reference call with at least this many alt reads


DEFAULT_QC = Qc()


def as_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def confidence(flags: list[str], site_pass: bool) -> str:
    """
    'low' when the caller filtered the site or the call (FILTER / FORMAT/FT) or
    the call has two or more flags, 'caution' with one flag, else 'high'.
    """
    if not site_pass or any(f.startswith("FT:") for f in flags) or len(flags) >= 2:
        return "low"
    return "caution" if flags else "high"


def parse_call(keys: dict[str, int], sample: str, qc: Qc, site_pass: bool = True, mito: bool = False) -> dict:
    """
    One genotype (bi-allelic record) with its quality: GT, zygosity (ref / het /
    hom / hemi / missing), alt-allele fraction (AD, else FORMAT/AF), depth, GQ,
    alt reads, DRAGEN FT / HML (homology region) / MSC (possible mosaic), the
    quality flags and the resulting confidence. On the mitochondrial genome
    (mito) the alt-allele fraction is the heteroplasmy level: no VAF_SKEW.
    """
    parts = sample.split(":")

    def value(key: str) -> str:
        i = keys.get(key, len(parts))
        return parts[i] if i < len(parts) else ""

    gt = value("GT") or "."
    alleles = re.split(r"[/|]", gt)
    called = [a for a in alleles if a != "."]
    alt = sum(1 for a in called if a != "0")
    zyg = ("missing" if not called else "ref" if not alt else "hemi" if len(alleles) == 1
           else "hom" if alt == len(alleles) else "het")
    ad = [as_int(x) for x in value("AD").split(",")] if value("AD") not in ("", ".") else []
    ad = ad if len(ad) >= 2 and None not in ad else []
    alt_reads = sum(ad[1:]) if ad else None
    vaf = (alt_reads / sum(ad) if ad and sum(ad) > 0
           else to_float(value("AF").split(",")[0]) if value("AF") not in ("", ".") else None)
    dp = as_int(value("DP")) if value("DP") not in ("", ".") else (sum(ad) if ad else None)
    gq = as_int(value("GQ"))
    ft = value("FT")

    flags = []
    if zyg != "missing":
        if ft not in ("", ".", "PASS"):
            flags.append(f"FT:{ft}")
        if dp is not None and dp < qc.min_dp:
            flags.append("LOW_DP")
        if gq is not None and gq < qc.min_gq:
            flags.append("LOW_GQ")
        if alt:
            if alt_reads is not None and alt_reads < qc.min_alt_reads:
                flags.append("LOW_ALT_READS")
            if not mito and vaf is not None and ((zyg == "het" and not qc.min_het_vaf <= vaf <= qc.max_het_vaf)
                                                 or (zyg != "het" and vaf < qc.min_hom_vaf)):
                flags.append("VAF_SKEW")
            if value("MSC") == "1":
                flags.append("MOSAIC")
        elif alt_reads is not None and alt_reads >= qc.ref_alt_reads:
            flags.append("ALT_READS_IN_REF")
        if value("HML") == "1":
            flags.append("HOMOLOGY")
    call = {"gt": gt, "zyg": zyg, "vaf": round(vaf, 3) if vaf is not None else None, "dp": dp, "gq": gq,
            "alt_reads": alt_reads, "ft": ft or None, "flags": flags,
            "conf": confidence(flags, site_pass) if zyg != "missing" else "missing"}
    # DRAGEN de novo calling with a pedigree: DN (DeNovo / Inherited / LowDQ) and DQ.
    if value("DN") not in ("", "."):
        call["dn"] = value("DN")
        call["dq"] = to_float(value("DQ")) if value("DQ") not in ("", ".") else None
    return call


def genotype_info(fmt: str, samples: list[str], names: list[str] | None = None, qc: Qc = DEFAULT_QC,
                  site_pass: bool = True, keep_all: bool = True, max_calls: int = 50, mito: bool = False) -> dict:
    """
    Carrier status of one bi-allelic record from its FORMAT and sample columns:
    carried (None without genotypes), het / hom / hemi carrier counts, and the
    parsed calls (see parse_call) of every sample (keep_all, families) or of the
    carriers only (cohorts; at most max_calls). Reference calls (0/0, ./.) - e.g.
    DRAGEN targeted-caller records at known pathogenic sites - are not carried.
    """
    keys = {k: i for i, k in enumerate(fmt.split(":"))}
    if "GT" not in keys or not samples:
        return {"carried": None, "calls": []}
    names = names or [f"sample{i + 1}" for i in range(len(samples))]
    n = Counter()
    calls = []
    for name, s in zip(names, samples):
        call = parse_call(keys, s, qc, site_pass, mito)
        if call["zyg"] in ("het", "hom", "hemi"):
            n[call["zyg"]] += 1
            if not keep_all and len(calls) < max_calls:
                calls.append({"sample": name, **call})
        if keep_all:
            calls.append({"sample": name, **call})
    return {"carried": bool(n), "het": n["het"], "hom": n["hom"], "hemi": n["hemi"], "calls": calls}


class Stats:
    """
    Streaming per-chunk statistics and the clinically relevant variants carried
    by the sample(s): ClinVar pathogenic / likely pathogenic (+ InterVar P/LP not
    in ClinVar) and rare coding/splicing ClinVar VUS / conflicting, each with the
    genotype calls and their quality. Files with at most family_max samples
    (single samples, families) also list every rare coding/splicing variant any
    sample carries, with all genotypes, for the inheritance analysis.
    """

    def __init__(self, header: list[str], qc: Qc = DEFAULT_QC, family_max: int = 10):
        def exact(name: str) -> int | None:
            return header.index(name) if name in header else None

        self.idx = {
            "chr": 0, "start": 1, "ref": 3, "alt": 4,
            "func": find_column(header, ("Func.refGene",), r"^Func\."),
            "gene": find_column(header, ("Gene.refGene",), r"^Gene\."),
            "exonic": find_column(header, ("ExonicFunc.refGene",), r"^ExonicFunc\."),
            "aachange": find_column(header, ("AAChange.refGene",), r"^AAChange\."),
            "genedetail": find_column(header, ("GeneDetail.refGene",), r"^GeneDetail\."),
            "clnsig": find_column(header, ("CLNSIG",), r"clnsig"),
            "clnrevstat": find_column(header, ("CLNREVSTAT",), r"clnrevstat"),
            "clndn": find_column(header, ("CLNDN",), r"^clndn$"),
            "clnalleleid": find_column(header, ("CLNALLELEID", "ALLELEID")),
            "intervar": exact("InterVar_automated"),
            "revel": find_column(header, ("REVEL", "REVEL_score"), r"^revel"),
            "cadd": find_column(header, ("CADD_phred", "CADD13_PHRED", "CADD_PHRED"), r"^cadd.*phred$"),
            "alphamissense": exact("AlphaMissense_pred"),
            "alphamissense_score": exact("AlphaMissense_score"),
            "splice": exact("dbscSNV_ADA_SCORE"),
            "vcf_chrom": exact("#CHROM"), "vcf_pos": exact("POS"), "vcf_ref": exact("REF"), "vcf_alt": exact("ALT"),
            "filter": exact("FILTER"),
            "format": exact("FORMAT"),
        }
        # PubMind (pubmind_db): four adjacent columns; 'PVID' and 'confidence' are too generic to look up alone.
        pm = exact("pathogenicity_sum")
        for key, offset, name in (("pubmind_pvid", -1, "PVID"), ("pubmind", 0, "pathogenicity_sum"),
                                  ("pubmind_score", 1, "paper_level_pathogenicity_score"),
                                  ("pubmind_confidence", 2, "confidence")):
            i = pm + offset if pm is not None else None
            self.idx[key] = i if i is not None and i < len(header) and header[i] == name else None
        self.gnomad = [i for i, name in enumerate(header) if GNOMAD_AF_RE.match(name)]
        # InterVar criteria follow InterVar_automated (PVS1 ... BP7).
        iv = self.idx["intervar"]
        self.acmg = [(i, header[i]) for i in range(iv + 1, min(iv + 1 + len(ACMG_CRITERIA), len(header)))
                     if header[i] in ACMG_CRITERIA] if iv is not None else []
        self.columns = {k: (header[v] if v is not None else None) for k, v in self.idx.items()
                        if k not in ("chr", "start", "ref", "alt")}
        self.columns["gnomad_af"] = ",".join(header[i] for i in self.gnomad) or None
        # Split up to FORMAT; the sample columns (huge for cohorts) stay one
        # field and are only parsed for the few listed variants.
        self.max_idx = max([v for v in self.idx.values() if v is not None] + self.gnomad)
        fmt_i = self.idx["format"]
        self.samples = header[fmt_i + 1:] if fmt_i is not None else []
        self.qc = qc
        self.family = 0 < len(self.samples) <= family_max
        self.counts = Counter()
        self.func = Counter()
        self.sample_counts: dict[str, Counter] = defaultdict(Counter)
        self.plp: list[dict] = []
        self.vus: list[dict] = []
        self.rare: list[dict] = []

    def col(self, f: list[str], key: str) -> str:
        i = self.idx[key]
        return f[i] if i is not None and i < len(f) else "."

    def gnomad_max(self, f: list[str]) -> float | None:
        """Highest gnomAD frequency (overall or any population, exome or genome); None if absent."""
        values = [v for v in (to_float(f[i]) for i in self.gnomad) if v is not None]
        return max(values) if values else None

    def genotypes(self, f: list[str]) -> dict:
        i = self.idx["format"]
        if i is None or i + 1 >= len(f):
            return {"carried": None, "calls": []}
        return genotype_info(f[i], f[i + 1].rstrip("\n").split("\t"), self.samples, self.qc,
                             self.col(f, "filter") in ("PASS", "."), keep_all=self.family, mito=f[0] in MITO)

    def carriers(self, f: list[str]) -> list[str]:
        """Samples whose GT has an alternative allele (het, hom, hemi); reads GT only, so it is
        cheap enough for every rare coding variant of a cohort."""
        i = self.idx["format"]
        keys = f[i].split(":") if i is not None and i + 1 < len(f) else []
        if "GT" not in keys:
            return []
        k = keys.index("GT")
        out = []
        for name, s in zip(self.samples, f[i + 1].rstrip("\n").split("\t")):
            parts = s.split(":")
            if k < len(parts) and any(a not in ("0", ".") for a in re.split(r"[/|]", parts[k])):
                out.append(name)
        return out

    def record(self, f: list[str], gt: dict, af: float | None) -> dict:
        """One listed variant with the fields the report shows."""
        acmg = [name for i, name in self.acmg if i < len(f) and f[i] == "1"]
        am = self.col(f, "alphamissense")
        am_score = self.col(f, "alphamissense_score")
        return {
            "chr": f[0], "start": f[1], "ref": f[3], "alt": f[4],
            "vcf": "-".join(self.col(f, k) for k in ("vcf_chrom", "vcf_pos", "vcf_ref", "vcf_alt")),
            "gene": self.col(f, "gene"), "func": self.col(f, "func"), "exonic": self.col(f, "exonic"),
            "aachange": first_transcript(self.col(f, "aachange")),
            # Splice sites have no AAChange; GeneDetail holds their c. notation, one transcript per
            # ';'-separated item (NM_1:exon3:c.194+2T>C;NM_2:...).
            "genedetail": (first_transcript(self.col(f, "genedetail").split(";")[0])
                           if splicing(self.col(f, "func")) else "."),
            "clnsig": self.col(f, "clnsig"), "clnrevstat": self.col(f, "clnrevstat").replace("\\x2c", ","),
            # Every condition ('|'-separated), in full: the report shows it on hover.
            "clndn": shorten(self.col(f, "clndn").replace("\\x2c", ",").replace("_", " "), 2000),
            "clnalleleid": self.col(f, "clnalleleid"),
            "intervar": self.col(f, "intervar"), "acmg": ",".join(acmg) or ".",
            "gnomad_af": f"{af:.3g}" if af is not None else ".",
            "revel": self.col(f, "revel"), "cadd": self.col(f, "cadd"),
            "alphamissense": f"{am} ({am_score})" if am not in (".", "") and am_score not in (".", "") else am,
            "splice": self.col(f, "splice"),
            "pubmind": self.col(f, "pubmind"),
            "pubmind_score": f"{score:.2f}" if (score := to_float(self.col(f, "pubmind_score"))) is not None else ".",
            "pubmind_confidence": self.col(f, "pubmind_confidence"),
            "pubmind_pvid": self.col(f, "pubmind_pvid"),
            "filter": self.col(f, "filter"),
            "carriers": (f"{gt['het']} het / {gt['hom'] + gt['hemi']} hom" if gt.get("carried") else "."),
            "genotypes": gt.get("calls", []),
        }

    def add(self, line: str) -> None:
        f = line.rstrip("\n").split("\t", self.max_idx + 1)
        c = self.counts
        c["total"] += 1
        ref, alt = f[3], f[4]
        if ref == "-" or alt == "-" or len(ref) != len(alt):
            c["indel"] += 1
        elif len(ref) == 1 and ref in NUCLEOTIDES and alt in NUCLEOTIDES:
            c["snv"] += 1
            c["ti" if (ref, alt) in TRANSITIONS else "tv"] += 1
        else:
            c["mnv"] += 1

        func, exonic = self.col(f, "func"), self.col(f, "exonic")
        if self.idx["func"] is not None:
            self.func[f"{func}\t{exonic}"] += 1
        if self.idx["filter"] is not None and f[self.idx["filter"]] in ("PASS", "."):
            c["pass"] += 1

        # Rare: the highest gnomAD frequency (any population) below 1%, or absent from gnomAD.
        af = self.gnomad_max(f)
        rare = af is None or af < RARE_AF
        coding = exonic in PROTEIN_ALTERING or splicing(func)
        coding_rare = bool(self.gnomad) and rare and coding
        if self.gnomad and rare:
            c["rare"] += 1
        if coding_rare:
            c["coding_rare"] += 1

        revel = to_float(self.col(f, "revel"))
        if revel is not None and revel >= 0.75:
            c["revel_ge_0.75"] += 1
        cadd = to_float(self.col(f, "cadd"))
        if cadd is not None and cadd >= 20:
            c["cadd_phred_ge_20"] += 1

        # Clinically relevant variants: only those the sample(s) carry are listed.
        clnsig = self.col(f, "clnsig")
        plp = self.idx["clnsig"] is not None and is_plp(clnsig)
        intervar_plp = self.col(f, "intervar").lower() in INTERVAR_PLP
        kind = vus_kind(clnsig) if self.idx["clnsig"] is not None else None
        vus = kind is not None and func in CODING_FUNC and rare
        # Rare coding/splicing variants of single samples and families: input of the inheritance analysis.
        candidate = self.family and rare and coding
        # What every sample carries, counted like the file totals.
        hits = [k for k, hit in (("coding_rare", coding_rare), ("clinvar_plp", plp),
                                 ("intervar_plp", intervar_plp and not plp), ("vus_candidates", vus)) if hit]
        if hits and self.samples:
            for s in self.carriers(f):
                self.sample_counts[s].update(hits)
        if not (plp or intervar_plp or vus or candidate):
            return
        gt = self.genotypes(f)
        if gt["carried"] is False:
            if plp:
                c["clinvar_plp_not_carried"] += 1
            return
        rec = self.record(f, gt, af)
        if plp or intervar_plp:
            c["clinvar_plp" if plp else "intervar_plp"] += 1
            if len(self.plp) < MAX_LISTED_PER_CHUNK:
                self.plp.append({"source": "ClinVar" if plp else "InterVar", **rec})
        if vus:
            c["vus_candidates"] += 1
            if len(self.vus) < MAX_LISTED_PER_CHUNK:
                self.vus.append({"kind": kind, **rec})
        if candidate and gt["carried"] and len(self.rare) < MAX_LISTED_PER_CHUNK:
            self.rare.append(rec)

    def to_dict(self) -> dict:
        return {"counts": dict(self.counts), "func": dict(self.func), "columns": self.columns,
                "samples": self.samples, "sample_counts": {s: dict(n) for s, n in self.sample_counts.items()},
                "qc": asdict(self.qc), "plp": self.plp, "vus": self.vus, "rare": self.rare}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--txt", required=True, help="ANNOVAR *_multianno.txt ('-' for STDIN)")
    p.add_argument("--vcf", required=True, help="VCF given to ANNOVAR (plain or .gz)")
    p.add_argument("--header-out", required=True, help="write the fixed header line here")
    p.add_argument("--out", default="-", help="data rows (no header) [STDOUT]")
    p.add_argument("--stats", help="write per-chunk statistics as JSON")
    p.add_argument("--sample", default="sample", help="sample / dataset id for the stats")
    p.add_argument("--chunk", default="all", help="chunk name for the stats")
    p.add_argument("--genome", default="", help="genome build, recorded in the stats")
    p.add_argument("--n-samples", type=int, default=None, help="number of genotype columns (stats only)")
    p.add_argument("--annovar-vcf", help="ANNOVAR *_multianno.vcf to write with the same gene separator")
    p.add_argument("--vcf-out", help="output for --annovar-vcf (plain VCF)")
    d = Qc()
    p.add_argument("--qc-min-dp", type=int, default=d.min_dp, help="LOW_DP below this read depth")
    p.add_argument("--qc-min-gq", type=int, default=d.min_gq, help="LOW_GQ below this genotype quality")
    p.add_argument("--qc-min-alt-reads", type=int, default=d.min_alt_reads,
                   help="LOW_ALT_READS: carriers with fewer alt-allele reads")
    p.add_argument("--qc-het-vaf", default=f"{d.min_het_vaf},{d.max_het_vaf}",
                   help="VAF_SKEW: heterozygous alt-allele fraction outside MIN,MAX")
    p.add_argument("--qc-min-hom-vaf", type=float, default=d.min_hom_vaf,
                   help="VAF_SKEW: homozygous / hemizygous alt-allele fraction below")
    p.add_argument("--family-max-samples", type=int, default=10,
                   help="list all genotypes and rare coding variants for files with at most this many samples")
    args = p.parse_args(argv)
    if bool(args.annovar_vcf) != bool(args.vcf_out):
        p.error("--annovar-vcf and --vcf-out go together")
    try:
        lo, hi = (float(x) for x in args.qc_het_vaf.split(","))
    except ValueError:
        p.error("--qc-het-vaf takes MIN,MAX, e.g. 0.2,0.8")
    args.qc = Qc(args.qc_min_dp, args.qc_min_gq, args.qc_min_alt_reads, lo, hi, args.qc_min_hom_vaf)
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    vcf_cols = read_vcf_columns(args.vcf)

    n_rows = 0
    header: list[str] | None = None
    stats: Stats | None = None
    gene_cols: list[int] = []
    src = open_text(args.txt)
    out = open_text(args.out, "w")
    try:
        txt_header = None
        for line in src:
            if txt_header is None:
                if not line.strip():
                    continue
                txt_header = line.rstrip("\r\n").split("\t")
                continue
            if header is None:
                n_data_cols = line.count("\t") + 1
                header = build_header(txt_header, vcf_cols, n_data_cols)
                stats = Stats(header, args.qc, args.family_max_samples)
                gene_cols = gene_columns(header)
            line = gene_separator(line, gene_cols)
            out.write(line)
            stats.add(line)
            n_rows += 1
        if txt_header is not None and header is None:
            # Header-only table (no variants in this chunk): assume the usual
            # 3 extra Otherinfo columns ANNOVAR writes before the VCF line.
            n_anno = next((i for i, c in enumerate(txt_header) if c.startswith("Otherinfo")), len(txt_header))
            header = build_header(txt_header, vcf_cols, max(len(txt_header), n_anno + 3 + len(vcf_cols)))
    finally:
        if out is not sys.stdout:
            out.close()
        if src is not sys.stdin:
            src.close()

    with open(args.header_out, "w", encoding="utf-8") as fh:
        if header is not None:
            fh.write("\t".join(header) + "\n")

    if args.annovar_vcf:
        with open_text(args.annovar_vcf) as vin, open(args.vcf_out, "w", encoding="utf-8") as vout:
            for line in vin:
                vout.write(vcf_gene_separator(line))

    if args.stats:
        payload = {
            "sample": args.sample,
            "chunk": args.chunk,
            "genome": args.genome,
            "n_genotype_columns": args.n_samples if args.n_samples is not None else max(len(vcf_cols) - 9, 0),
            **(stats.to_dict() if stats else {"counts": {"total": 0}, "func": {}, "columns": {}}),
        }
        with open(args.stats, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=1)

    print(f"[format_annovar_table] {args.sample}:{args.chunk} rows={n_rows} "
          f"columns={len(header) if header else 0}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ValueError as exc:
        print(f"[format_annovar_table] ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
