#!/usr/bin/env python3
"""
gene_info.py
------------
Gene-level information for the report and the gene panels of --gene_panel.

The gene table (annovar_nf/genes.tsv in humandb, built by `annovar_db.py
genes`) joins three open sources:
    GenCC        gene-disease validity: diseases, modes of inheritance and the
                 strength of the evidence, from every GenCC submitter (CC0)
    gnomAD v4.1  loss-of-function (LOEUF, pLI) and missense (Z) constraint of
                 the MANE Select transcript (CC0)
    HGNC         approved symbols, gene names and previous symbols, so genes
                 named differently in ANNOVAR's refGene or in a panel still match

Table columns: gene, hgnc_id, name, moi (modes of inheritance of the diseases
with at least supportive evidence: AD, AR, XL ...), diseases
('title|moi|classification|submitters', best evidence first, ';'-separated),
loeuf, pli, mis_z, previous (earlier symbols, ','-separated).

A gene panel is a list of gene symbols, one per line (other columns, '#'
comments and a header are ignored), or a PanelApp TSV export, of which the
green genes (GEL_Status 3 or 4) are used.

Only the Python standard library is used.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

# In a folder of its own: rebuilding it leaves the humandb directory itself unchanged, which keeps
# the -resume cache of the ANNOVAR tasks (Nextflow hashes the directory with its modification time).
TABLE_NAME = "annovar_nf/genes.tsv"
COLUMNS = ("gene", "hgnc_id", "name", "moi", "diseases", "loeuf", "pli", "mis_z", "previous")
MOI_ABBREV = {"Autosomal dominant": "AD", "Autosomal recessive": "AR", "X-linked": "XL",
              "X-linked dominant": "XLD", "X-linked recessive": "XLR", "Semidominant": "SD",
              "Mitochondrial": "MT", "Y-linked": "YL", "Somatic mosaicism": "mosaic"}
MOI_ORDER = ("AD", "AR", "SD", "XL", "XLD", "XLR", "YL", "MT", "mosaic")
# GenCC classifications kept, strongest first; disputed / refuted / no known relationship are left out.
CLASS_RANK = {"Definitive": 0, "Strong": 1, "Moderate": 2, "Supportive": 3, "Limited": 4}
SUPPORTED = 3                  # a disease counts for the gene's modes of inheritance up to 'Supportive'
LOEUF_CONSTRAINED = 0.6        # gnomAD v4: LOEUF below 0.6 - intolerant of loss of function
# Inheritance models of inheritance.py and the gene modes of inheritance they fit.
MODEL_MOI = {"de novo": {"AD", "SD", "XL", "XLD"}, "dominant (affected parent)": {"AD", "SD", "XL", "XLD"},
             "homozygous (AR)": {"AR", "SD"}, "compound heterozygous": {"AR"},
             "possible compound heterozygous": {"AR"}, "X-linked (hemizygous)": {"XL", "XLR", "XLD"}}
GENE_SPLIT = re.compile(r"[;,]")
HEADER_WORDS = ("gene symbol", "gene_symbol", "hgnc symbol", "symbol", "gene", "genes")


def clean(text: str) -> str:
    """A value safe inside the table's ';' / '|' separated fields."""
    return re.sub(r"[\t;|]+", ",", (text or "").strip())


def number(value: str) -> float | None:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if x == x else None          # NaN


# ======================================================================================
# Building the table
# ======================================================================================
def read_gencc(fh) -> dict[str, dict]:
    """GenCC submissions TSV -> {symbol: {'hgnc_id', 'diseases': {key: {title, moi, rank, submitters}}}}."""
    csv.field_size_limit(1 << 24)          # long notes fields (sys.maxsize overflows on some platforms)
    genes: dict[str, dict] = {}
    for row in csv.DictReader(fh, delimiter="\t"):
        rank = CLASS_RANK.get((row.get("classification_title") or "").strip())
        symbol = (row.get("gene_symbol") or "").strip()
        if rank is None or not symbol:
            continue
        gene = genes.setdefault(symbol, {"hgnc_id": (row.get("gene_curie") or "").strip(), "diseases": {}})
        title = clean(row.get("disease_title") or row.get("disease_original_title") or "")
        key = (row.get("disease_curie") or "").strip() or title.lower()
        d = gene["diseases"].setdefault(key, {"title": title, "moi": set(), "rank": rank, "submitters": set()})
        moi = MOI_ABBREV.get((row.get("moi_title") or "").strip())
        if moi:
            d["moi"].add(moi)
        d["rank"] = min(d["rank"], rank)
        d["submitters"].add((row.get("submitter_title") or "").strip())
    return genes


def read_constraint(fh) -> dict[str, tuple[str, str, str]]:
    """gnomAD constraint metrics -> {symbol: (LOEUF, pLI, missense Z)} of the MANE Select transcript
    (RefSeq row first), else the canonical one."""
    reader = csv.DictReader(fh, delimiter="\t")
    best: dict[str, tuple[int, tuple]] = {}
    for row in reader:
        mane, canonical = row.get("mane_select") == "true", row.get("canonical") == "true"
        if not (mane or canonical):
            continue
        priority = 0 if mane and row.get("transcript", "").startswith("NM_") else 1 if mane else 2
        symbol = row.get("gene", "")
        if symbol and (symbol not in best or priority < best[symbol][0]):
            best[symbol] = (priority, tuple(fmt_metric(row.get(k)) for k in
                                            ("lof.oe_ci.upper", "lof.pLI", "mis.z_score")))
    return {s: v for s, (_p, v) in best.items()}


def fmt_metric(value) -> str:
    x = number(value)
    return f"{x:.2f}" if x is not None else "."


def read_hgnc(fh) -> dict[str, dict]:
    """HGNC complete set -> {approved symbol: {'hgnc_id', 'name', 'previous': [...]}}."""
    out = {}
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("status") != "Approved":
            continue
        previous = [p.strip() for p in (row.get("prev_symbol") or "").strip('"').split("|") if p.strip()]
        out[row["symbol"]] = {"hgnc_id": row.get("hgnc_id", ""), "name": clean(row.get("name", "")),
                              "previous": previous}
    return out


def build_rows(gencc: dict, constraint: dict, hgnc: dict) -> list[list[str]]:
    """One row per gene with GenCC or constraint data, under its approved HGNC symbol."""
    current = {}                                # previous symbol -> approved one (unique only)
    for symbol, h in hgnc.items():
        for p in h["previous"]:
            current[p] = None if p in current else symbol
    current = {p: s for p, s in current.items() if s and p not in hgnc}

    def approved(symbol: str) -> str:
        return symbol if symbol in hgnc else current.get(symbol, symbol)

    metrics = {}
    for symbol, values in constraint.items():
        metrics.setdefault(approved(symbol), values)
    rows = []
    for symbol in sorted(set(gencc) | set(metrics)):
        if approved(symbol) != symbol and approved(symbol) in (set(gencc) | set(metrics)):
            continue                            # an outdated symbol of a gene that is listed anyway
        g = gencc.get(symbol, {"hgnc_id": "", "diseases": {}})
        h = hgnc.get(symbol, {})
        diseases = sorted(g["diseases"].values(), key=lambda d: (d["rank"], -len(d["submitters"]), d["title"]))
        moi = {m for d in diseases if d["rank"] <= SUPPORTED for m in d["moi"]}
        loeuf, pli, mis_z = metrics.get(symbol, (".", ".", "."))
        rows.append([
            symbol, h.get("hgnc_id") or g["hgnc_id"] or ".", h.get("name") or ".",
            ",".join(m for m in MOI_ORDER if m in moi) or ".",
            ";".join(f"{d['title']}|{','.join(m for m in MOI_ORDER if m in d['moi']) or '?'}|"
                     f"{next(k for k, r in CLASS_RANK.items() if r == d['rank'])}|{len(d['submitters'])}"
                     for d in diseases) or ".",
            loeuf, pli, mis_z, ",".join(h.get("previous", [])) or "."])
    return rows


def write_table(rows: list[list[str]], path: Path, sources: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        for s in sources:
            out.write(f"## {s}\n")
        out.write("#" + "\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")


# ======================================================================================
# Using the table
# ======================================================================================
@dataclass
class Gene:
    gene: str
    hgnc_id: str
    name: str
    moi: list
    diseases: list                # [{'title', 'moi', 'classification', 'submitters'}]
    loeuf: float | None
    pli: float | None
    mis_z: float | None

    def summary(self) -> str:
        """'AD · LOEUF 0.09' (what the report shows next to the gene)."""
        return " · ".join(x for x in ("/".join(self.moi), f"LOEUF {self.loeuf:.2f}" if self.loeuf is not None else "")
                          if x)

    def disease_text(self) -> str:
        """'Dravet syndrome (AD, Definitive); ...' for the TSV files and tooltips."""
        return "; ".join(f"{d['title']} ({d['moi']}, {d['classification']})" for d in self.diseases)

    def detail(self) -> str:
        """Every fact of the gene, one per line (tooltip)."""
        lines = [f"{self.gene}{f' ({self.name})' if self.name else ''}"]
        metrics = [f"{k} {v:.2f}" for k, v in (("LOEUF", self.loeuf), ("pLI", self.pli), ("missense Z", self.mis_z))
                   if v is not None]
        if metrics:
            lines.append("gnomAD constraint: " + ", ".join(metrics)
                         + (" - intolerant of loss of function" if self.loeuf is not None
                            and self.loeuf < LOEUF_CONSTRAINED else ""))
        lines += [f"{d['title']} - {d['moi']}, {d['classification']} ({d['submitters']} GenCC submitter"
                  f"{'s' if d['submitters'] != 1 else ''})" for d in self.diseases]
        return "\n".join(lines)


@dataclass
class GeneTable:
    genes: dict = field(default_factory=dict)          # approved symbol (upper case) -> Gene
    aliases: dict = field(default_factory=dict)        # previous symbol (upper case) -> approved
    sources: list = field(default_factory=list)

    def symbol(self, name: str) -> str:
        """Approved symbol (upper case) of a gene name, previous symbols resolved."""
        key = name.strip().upper()
        return key if key in self.genes else self.aliases.get(key, key)

    def get(self, name: str) -> Gene | None:
        return self.genes.get(self.symbol(name))

    def lookup(self, gene_field: str) -> Gene | None:
        """The first gene of an ANNOVAR Gene.* value ('A;B') with information."""
        return next((g for g in (self.get(x) for x in genes_of(gene_field)) if g), None)


def genes_of(gene_field: str | None) -> list[str]:
    return [g.strip() for g in GENE_SPLIT.split(gene_field or "") if g.strip() not in ("", ".", "NONE")]


def load(path: str | None) -> GeneTable:
    """The gene table (an empty one without a path)."""
    table = GeneTable()
    if not path:
        return table
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                table.sources.append(line[2:].strip())
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = dict(zip(COLUMNS, line.rstrip("\n").split("\t")))
            diseases = []
            for item in (f.get("diseases") or ".").split(";"):
                parts = item.split("|")
                if len(parts) == 4:
                    diseases.append({"title": parts[0], "moi": parts[1], "classification": parts[2],
                                     "submitters": int(parts[3]) if parts[3].isdigit() else 0})
            key = f["gene"].upper()
            table.genes[key] = Gene(f["gene"], f.get("hgnc_id", "."), "" if f.get("name") == "." else f.get("name", ""),
                                    [m for m in (f.get("moi") or ".").split(",") if m != "."], diseases,
                                    number(f.get("loeuf")), number(f.get("pli")), number(f.get("mis_z")))
            for p in (f.get("previous") or ".").split(","):
                if p not in (".", ""):
                    table.aliases.setdefault(p.upper(), key)
    table.aliases = {p: s for p, s in table.aliases.items() if p not in table.genes}
    return table


def moi_fit(model: str, gene: Gene | None) -> bool | None:
    """Does an inheritance model fit the gene's modes of inheritance? None when the gene has none."""
    if gene is None or not gene.moi or model not in MODEL_MOI:
        return None
    return bool(MODEL_MOI[model] & set(gene.moi))


# ======================================================================================
# Gene panels
# ======================================================================================
@dataclass
class Panel:
    name: str
    genes: set                    # upper-case symbols, previous symbols resolved


def read_panel(path: str, table: GeneTable | None = None) -> Panel:
    """A list of gene symbols (first column) or a PanelApp TSV export (green genes)."""
    name = Path(path).name.split(".")[0]
    with open(path, encoding="utf-8-sig") as fh:
        rows = [line.rstrip("\r\n") for line in fh if line.strip() and not line.startswith("#")]
    header = [h.strip().lower() for h in rows[0].split("\t")] if rows else []
    symbol_col = next((header.index(k) for k in HEADER_WORDS if k in header), None) if len(header) > 1 else None
    def cell(f: list[str], i: int | None) -> str:
        return f[i].strip() if i is not None and i < len(f) else ""

    genes = set()
    if symbol_col is not None:                                  # a table with a header
        col = {h: i for i, h in enumerate(header)}
        status, entity, level4, version = (col.get(k) for k in ("gel_status", "entity type", "level4", "version"))
        for row in rows[1:]:
            f = row.split("\t")
            if entity is not None and cell(f, entity).lower() not in ("gene", ""):
                continue
            if status is not None and (number(cell(f, status)) or 0) < 3:     # PanelApp: green genes only
                continue
            if cell(f, symbol_col):
                genes.add(cell(f, symbol_col))
            if cell(f, level4):
                name = cell(f, level4) + (f" (PanelApp v{cell(f, version)})" if cell(f, version) else "")
    else:                                                       # a gene list: the first word of every line
        genes = {w for w in (re.split(r"[\s,;]+", r.strip())[0] for r in rows) if w.lower() not in HEADER_WORDS}
    resolve = table.symbol if table else (lambda g: g.upper())
    return Panel(name, {resolve(g) for g in genes if g})


def panels_of(gene_field: str | None, panels: list[Panel], table: GeneTable | None = None) -> list[str]:
    """Names of the panels holding any gene of an ANNOVAR Gene.* value."""
    resolve = table.symbol if table else (lambda g: g.upper())
    symbols = {resolve(g) for g in genes_of(gene_field)}
    return [p.name for p in panels if symbols & p.genes]
