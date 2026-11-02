"""
inheritance.py
--------------
Inheritance models for the rare coding/splicing variants of a single sample or
a family (one input file), used by summarize_annotations.py.

Probands are the affected samples of the PED file; with a pedigree without
affected status (VCF header, inferred) the children; without a pedigree a
single-sample file's sample. For every proband the variants it carries are
matched against:

  de novo                 heterozygous (hemizygous on chrX in males), both
                          parents called homozygous reference without reads
                          for the variant; DRAGEN's de novo call (DN=DeNovo)
                          when the VCF has one, else no low-confidence call
                          in the trio
  homozygous (AR)         homozygous; parents present: neither homozygous;
                          neither parent sequenced: damaging variants only
                          (see below)
  X-linked (hemizygous)   male, chrX outside the PARs
  dominant (affected parent)  heterozygous, carried by an affected parent
  compound heterozygous   two or more heterozygous variants in one gene:
                          phased by the parents (one from each) or, without
                          both parents, 'possible' (phase unknown; damaging
                          variants only: loss of function, splicing,
                          REVEL >= 0.5 or CADD >= 20)

and must segregate: an unaffected family member with the same genotype
(homozygous for AR, hemizygous for X-linked, carrier for de novo / dominant,
both variants for compound het) excludes the candidate; so does a parent not
known to be affected who is homozygous for the variant (recessive models). gnomAD maximum allele
frequency limits: below DOMINANT_AF for de novo / dominant, below RECESSIVE_AF
for the recessive models.

Doubtful candidates are kept but marked ('doubtful': the reason), so that the
report can hide them by default: the proband's genotype call has low
confidence, or the genotypes break Mendelian inheritance (a homozygous child
of a parent called homozygous reference: usually a genotyping error, rarely a
deletion on the other allele or uniparental disomy).
"""

from __future__ import annotations

from collections import defaultdict

DOMINANT_AF = 0.001
RECESSIVE_AF = 0.01
X_NON_PAR = {"hg19": (2699521, 154931043), "hg38": (2781480, 155701382)}
CARRIER = ("het", "hom", "hemi")
LOF = ("stopgain", "stoploss", "startloss", "frameshift insertion", "frameshift deletion",
       "frameshift substitution", "frameshift block substitution")
PROTEIN_ALTERING = LOF + ("nonsynonymous SNV", "nonframeshift insertion", "nonframeshift deletion",
                          "nonframeshift substitution")
# Order of the models in the report.
MODEL_ORDER = ("de novo", "homozygous (AR)", "X-linked (hemizygous)", "compound heterozygous",
               "possible compound heterozygous", "dominant (affected parent)")
# Within a model: confident calls first.
CONFIDENCE_ORDER = {"high": 0, "caution": 1, "low": 2}


def to_float(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def x_non_par(v: dict, genome: str) -> bool:
    chrom = str(v.get("chr", "")).removeprefix("chr")
    if chrom != "X":
        return False
    start, end = X_NON_PAR.get(genome, X_NON_PAR["hg38"])
    pos = int(to_float(v.get("start")) or 0)
    return start <= pos <= end


def splicing(v: dict) -> bool:
    """Splice site of a protein-coding transcript (Func 'splicing' / 'exonic;splicing', not ncRNA_splicing)."""
    return "splicing" in str(v.get("func", "")).split(";")


def coding(v: dict) -> bool:
    return v.get("exonic") in PROTEIN_ALTERING or splicing(v)


def severity(v: dict) -> int:
    """Sort key: loss of function, then damaging missense (REVEL >= 0.5 or CADD >= 20), then the rest."""
    if v.get("exonic") in LOF or splicing(v):
        return 3
    revel, cadd = to_float(v.get("revel")), to_float(v.get("cadd"))
    return 2 if (revel or 0) >= 0.5 or (cadd or 0) >= 20 else 1


def genes(v: dict) -> list[str]:
    return [g for g in str(v.get("gene", "")).split(";") if g not in ("", ".", "NONE")]


def probands(samples: list[str], ped: dict[str, dict]) -> list[str]:
    """
    Affected samples of the PED file. Without affected status - a pedigree from
    the VCF header or inferred from the data - the children (a parent in the
    file) are the probands; a single sample without pedigree is its own proband.
    """
    affected = [s for s in samples if (ped.get(s) or {}).get("affected")]
    if affected:
        return affected
    if any((e or {}).get("source") in ("vcf header", "inferred") for e in ped.values()):
        return [s for s in samples if (ped.get(s) or {}).get("father") in samples
                or (ped.get(s) or {}).get("mother") in samples]
    if len(samples) == 1 and not ped.get(samples[0]):
        return samples
    return []


def dragen_de_novo(g: dict) -> str:
    """' · DRAGEN DeNovo (DQ 25)' from the call's FORMAT/DN and DQ, if present."""
    if not g.get("dn"):
        return ""
    return f" · DRAGEN {g['dn']}" + (f" (DQ {g['dq']:g})" if g.get("dq") is not None else "")


def analyse(variants: list[dict], samples: list[str], ped: dict[str, dict], sex: dict[str, str],
            genome: str) -> list[dict]:
    """
    Candidates for every proband: dicts with sample, model, gene, the variant
    (and partner variants for compound het), a note and 'doubtful' (the reason,
    or '').
    variants: rare coding/splicing records with 'genotypes' for every sample.
    ped: PED entries of the file's samples; sex: inferred sex (used when the
    PED file does not give it).
    """
    out = []
    for p in probands(samples, ped):
        info = ped.get(p) or {}
        father = info.get("father") if info.get("father") in samples else None
        mother = info.get("mother") if info.get("mother") in samples else None
        p_sex = info.get("sex") if info.get("sex") in ("male", "female") else sex.get(p, "unknown")
        unaffected = [s for s in samples if s != p and (ped.get(s) or {}).get("affected") is False]
        affected = {s for s in samples if (ped.get(s) or {}).get("affected")}

        def call(v: dict, s: str | None) -> dict:
            return next((g for g in v.get("genotypes", []) if g.get("sample") == s), {}) if s else {}

        def male(s: str) -> bool:
            s_sex = (ped.get(s) or {}).get("sex")
            return (s_sex if s_sex in ("male", "female") else sex.get(s)) == "male"

        # Parents not known to be affected count as unaffected for the genotype
        # rules (a homozygous parent rules out a recessive variant).
        healthy_parents = [par for par in (father, mother) if par and par not in affected]

        def de_novo_ok(v: dict, p: str = p, father: str | None = father, mother: str | None = mother) -> bool:
            """DRAGEN's verdict when it scored de novo calls; else a confident trio without parental alt reads."""
            child = call(v, p)
            trio = [child, call(v, father), call(v, mother)]
            if any("ALT_READS_IN_REF" in (g.get("flags") or []) for g in trio[1:]):
                return False
            if child.get("dn"):
                return child["dn"] == "DeNovo"
            return all(g.get("conf") != "low" for g in trio)

        het_by_gene = defaultdict(list)
        for v in variants:
            gp = call(v, p).get("zyg")
            if gp not in CARRIER or not coding(v):
                continue
            af = to_float(v.get("gnomad_af")) or 0.0
            on_x = x_non_par(v, genome)
            gf, gm = call(v, father).get("zyg"), call(v, mother).get("zyg")
            carriers = {s for s in samples if call(v, s).get("zyg") in CARRIER}
            hom_parents = [par for par in healthy_parents if call(v, par).get("zyg") in ("hom", "hemi")]

            def add(model: str, note: str = "", excluded_by: list[str] = (), doubtful: str = "",
                    v: dict = v, p: str = p) -> None:
                if excluded_by:
                    return
                out.append({"sample": p, "model": model, "gene": v.get("gene", "."), "variant": v,
                            "partners": [], "note": note, "doubtful": doubtful})

            dn = dragen_de_novo(call(v, p))
            if on_x and p_sex == "male":
                if gp in ("hemi", "hom"):
                    if father and mother and gm == "ref" and gf == "ref" and af < DOMINANT_AF:
                        if de_novo_ok(v):
                            add("de novo", "hemizygous, mother and father reference" + dn,
                                [s for s in unaffected if s in carriers])
                    elif af < RECESSIVE_AF:
                        add("X-linked (hemizygous)", f"mother {gm or 'not sequenced'}",
                            [s for s in unaffected if male(s) and call(v, s).get("zyg") in ("hemi", "hom")]
                            + [par for par in hom_parents if par == father])
                continue

            if gp == "het" and af < DOMINANT_AF:
                if father and mother and gf == "ref" and gm == "ref":
                    if de_novo_ok(v):
                        add("de novo", "father and mother reference" + dn, [s for s in unaffected if s in carriers])
                elif any(par in affected and call(v, par).get("zyg") in CARRIER for par in (father, mother) if par):
                    add("dominant (affected parent)", "carried by an affected parent",
                        [s for s in unaffected if s in carriers])
            # Without parents a rare homozygous call is common (population-specific
            # variants, poorly covered regions): damaging variants only.
            if gp == "hom" and af < RECESSIVE_AF and (father or mother or severity(v) >= 2):
                parents = ", ".join(f"{r} {z}" for r, z in (("father", gf), ("mother", gm)) if z)
                ref_parents = [r for r, z in (("father", gf), ("mother", gm)) if z == "ref"]
                add("homozygous (AR)", parents or "parents not sequenced",
                    [s for s in unaffected if call(v, s).get("zyg") == "hom"] + hom_parents,
                    f"Mendelian inconsistency: {' and '.join(ref_parents)} homozygous reference "
                    "(genotyping error, or a deletion on the other allele / uniparental disomy)" if ref_parents else "")
            if gp == "het" and af < RECESSIVE_AF and not hom_parents:
                for g in genes(v):
                    het_by_gene[g].append(v)

        for gene, vs in het_by_gene.items():
            if len(vs) < 2:
                continue
            if not (father and mother):
                # Unphased pairs are common (polymorphic genes): damaging variants only.
                vs = [v for v in vs if severity(v) >= 2]
                if len(vs) < 2:
                    continue
            # A pair carried together by an unaffected relative does not segregate.
            if any(all(call(v, s).get("zyg") in CARRIER for v in vs) for s in unaffected):
                continue

            def origin(v: dict, father: str | None = father, mother: str | None = mother) -> str:
                f, m = call(v, father).get("zyg"), call(v, mother).get("zyg")
                if f in CARRIER and m == "ref":
                    return "paternal"
                if m in CARRIER and f == "ref":
                    return "maternal"
                if f == "ref" and m == "ref":
                    return "de novo"
                return "unknown"

            if father and mother:
                origins = {id(v): origin(v) for v in vs}
                kinds = set(origins.values())
                if {"paternal", "maternal"} <= kinds:
                    model, note = "compound heterozygous", "one variant from each parent"
                elif "de novo" in kinds and kinds & {"paternal", "maternal"}:
                    model, note = "possible compound heterozygous", "one variant de novo, phase unknown"
                else:
                    continue
                vs = [v for v in vs if origins[id(v)] != "unknown"]
            elif father or mother:
                parent = father or mother
                in_parent = [v for v in vs if call(v, parent).get("zyg") in CARRIER]
                if not in_parent or len(in_parent) == len(vs):
                    continue        # all from the sequenced parent (cis) or none phased
                model = "possible compound heterozygous"
                note = f"one variant from the {'father' if father else 'mother'}"
            else:
                model, note = "possible compound heterozygous", "phase unknown (parents not sequenced)"
            vs = sorted(vs, key=lambda v: -severity(v))
            for v in vs:
                out.append({"sample": p, "model": model, "gene": gene, "variant": v,
                            "partners": [w.get("vcf") for w in vs if w is not v], "note": note, "doubtful": ""})

    for c in out:
        if not c["doubtful"] and call_conf(c) == "low":
            c["doubtful"] = "low-confidence genotype call"
    out.sort(key=lambda c: (c["sample"], MODEL_ORDER.index(c["model"]), bool(c["doubtful"]),
                            CONFIDENCE_ORDER.get(call_conf(c), 3), -severity(c["variant"]), c["gene"]))
    return out


def call_conf(c: dict) -> str:
    """Confidence of the proband's genotype call of a candidate."""
    return next((g.get("conf", "") for g in c["variant"].get("genotypes", []) if g.get("sample") == c["sample"]), "")
