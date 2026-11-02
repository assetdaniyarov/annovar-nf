#!/usr/bin/env python3
"""
detect_genome_build.py
----------------------
Infer the reference genome build (hg19 or hg38) of a VCF file and reconcile it
with the build requested by the user.

Evidence is collected from three sources, in decreasing order of reliability:

  1. contig  - ##contig=<ID=...,length=...> lines; the lengths of the primary
               chromosomes differ between GRCh37/hg19 and GRCh38/hg38, which
               makes this the only truly reliable signal.
  2. header  - free-text hints in ##reference=, ##contig=<assembly=...> and
               DRAGEN/GATK command lines (e.g. "/references/hg38.fa").
  3. filename- build tokens in the file name (e.g. "cohort_hg19_dragen.vcf").

Decision rules
  * --expected hg19|hg38 (explicit): the explicit build wins. If the contig
    lengths contradict it the script fails (annotating hg19 coordinates with
    hg38 databases silently produces wrong results) unless --no-strict is set.
  * --expected auto: the most reliable available source wins. A conflict
    between contig lengths and the file name is reported as an error, because
    it usually means a mislabelled or partially lifted-over file.

The resolved build is printed to STDOUT (and nothing else), diagnostics go to
STDERR, and a machine-readable JSON report is written with --report.

Usage:
    bcftools view -h input.vcf.gz > header.txt
    detect_genome_build.py --header header.txt --filename input.vcf.gz \
        --expected auto --report sample.build.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

BUILDS = ("hg19", "hg38")

# Lengths of primary chromosomes (UCSC and Ensembl/NCBI naming share lengths).
# chrM is deliberately excluded: hg19 uses chrM=16571 while GRCh37 (b37) and
# both hg38 flavours use the rCRS (16569), so it is not discriminative.
CONTIG_LENGTHS = {
    "hg19": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    },
    "hg38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    },
}

# Free-text build tokens. Bounded by non-alphanumerics so that e.g. "hg190"
# or "sample_b370" do not match.
TOKEN_PATTERNS = {
    "hg19": re.compile(r"(?<![a-z0-9])(hg19|grch37|b37|hs37d5|hs37|hg37)(?![a-z0-9])", re.I),
    "hg38": re.compile(r"(?<![a-z0-9])(hg38|grch38|b38|hs38dh|hs38|hg38dh)(?![a-z0-9])", re.I),
}

CONTIG_RE = re.compile(r"^##contig=<(.*)>\s*$")
HINT_LINE_RE = re.compile(r"^##(reference|contig|DRAGENCommandLine|GATKCommandLine|source)\b", re.I)


class BuildError(Exception):
    """Raised when the build cannot be determined or is contradictory."""


def parse_contig_attrs(body: str) -> dict[str, str]:
    """Parse 'ID=chr1,length=123,assembly=hg19' into a dict (quotes stripped)."""
    attrs = {}
    for match in re.finditer(r'(\w+)=("[^"]*"|[^,]*)', body):
        attrs[match.group(1).lower()] = match.group(2).strip('"')
    return attrs


def evidence_from_contigs(header_lines: list[str]) -> dict[str, int]:
    """Count primary contigs whose declared length matches each build."""
    votes = {b: 0 for b in BUILDS}
    for line in header_lines:
        m = CONTIG_RE.match(line)
        if not m:
            continue
        attrs = parse_contig_attrs(m.group(1))
        name = re.sub(r"^chr", "", attrs.get("id", ""), flags=re.I)
        try:
            length = int(attrs.get("length", ""))
        except ValueError:
            continue
        for build in BUILDS:
            if CONTIG_LENGTHS[build].get(name) == length:
                votes[build] += 1
    return votes


def evidence_from_text(text: str) -> set[str]:
    """Return the set of builds whose tokens appear in *text*."""
    return {b for b, pat in TOKEN_PATTERNS.items() if pat.search(text)}


def evidence_from_header_hints(header_lines: list[str]) -> set[str]:
    found = set()
    for line in header_lines:
        if HINT_LINE_RE.match(line):
            found |= evidence_from_text(line)
    return found


def single(builds: set[str]) -> str | None:
    """Return the build if exactly one was found, otherwise None."""
    return next(iter(builds)) if len(builds) == 1 else None


def resolve(header_lines: list[str], filename: str, expected: str, strict: bool) -> tuple[str, dict]:
    """Resolve the genome build. Returns (build, report)."""
    votes = evidence_from_contigs(header_lines)
    contig_builds = {b for b, n in votes.items() if n > 0}
    hint_builds = evidence_from_header_hints(header_lines)
    name_builds = evidence_from_text(Path(filename).name)

    report = {
        "filename": Path(filename).name,
        "expected": expected,
        "evidence": {
            "contig_length_votes": votes,
            "header_hints": sorted(hint_builds),
            "filename_tokens": sorted(name_builds),
        },
        "warnings": [],
    }
    warn = report["warnings"].append

    if len(contig_builds) > 1:
        msg = f"contig lengths match both builds ({votes}); header looks mixed"
        if strict:
            raise BuildError(msg)
        warn(msg)
    contig = single(contig_builds)
    hint = single(hint_builds)
    name = single(name_builds)

    # ---- explicit build requested -------------------------------------------
    if expected in BUILDS:
        if contig and contig != expected:
            msg = (f"--genome/samplesheet says '{expected}' but ##contig lengths indicate "
                   f"'{contig}'. Annotating with the wrong build gives wrong results. "
                   f"Fix the samplesheet, or use --skip_build_check if the header is known to be stale.")
            if strict:
                raise BuildError(msg)
            warn(msg)
        for source, value in (("header hints", hint), ("file name", name)):
            if value and value != expected:
                warn(f"{source} suggest '{value}' but '{expected}' was requested; using '{expected}'")
        report.update(build=expected, source="explicit")
        return expected, report

    # ---- automatic detection ------------------------------------------------
    if contig:
        for source, value in (("header hints", hint), ("file name", name)):
            if value and value != contig:
                msg = (f"##contig lengths indicate '{contig}' but {source} suggest '{value}'. "
                       f"Set the build explicitly (samplesheet 'genome' column or --genome), "
                       f"or use --skip_build_check to trust the contig lengths.")
                if strict:
                    raise BuildError(msg)
                warn(msg)
        report.update(build=contig, source="contig_lengths")
        return contig, report

    for source, value in (("header_hints", hint), ("filename", name)):
        if value:
            warn(f"no ##contig lengths available; build inferred from {source.replace('_', ' ')} only")
            report.update(build=value, source=source)
            return value, report

    raise BuildError(
        "cannot determine the genome build: no ##contig lengths, header hints or file-name "
        "tokens (hg19/GRCh37/hg38/GRCh38). Add a 'genome' column to the samplesheet or pass --genome."
    )


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--header", required=True, help="VCF header (output of 'bcftools view -h')")
    p.add_argument("--filename", required=True, help="original input file name")
    p.add_argument("--expected", default="auto", help="auto, hg19 or hg38 [auto]")
    p.add_argument("--no-strict", dest="strict", action="store_false",
                   help="downgrade contradictions to warnings")
    p.add_argument("--report", help="write a JSON report to this path")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    expected = (args.expected or "auto").strip().lower()
    expected = {"grch37": "hg19", "b37": "hg19", "grch38": "hg38", "": "auto", "null": "auto"}.get(expected, expected)
    if expected not in BUILDS + ("auto",):
        print(f"[detect_genome_build] ERROR: unsupported genome '{args.expected}'", file=sys.stderr)
        return 2

    with open(args.header, encoding="utf-8", errors="replace") as fh:
        header_lines = [line.rstrip("\n") for line in fh if line.startswith("##")]

    try:
        build, report = resolve(header_lines, args.filename, expected, args.strict)
    except BuildError as exc:
        print(f"[detect_genome_build] ERROR ({Path(args.filename).name}): {exc}", file=sys.stderr)
        if args.report:
            Path(args.report).write_text(json.dumps({"filename": args.filename, "error": str(exc)}, indent=2))
        return 1

    for w in report["warnings"]:
        print(f"[detect_genome_build] WARNING ({report['filename']}): {w}", file=sys.stderr)
    print(f"[detect_genome_build] {report['filename']}: {build} (source: {report['source']})", file=sys.stderr)
    if args.report:
        Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(build)
    return 0


if __name__ == "__main__":
    sys.exit(main())
