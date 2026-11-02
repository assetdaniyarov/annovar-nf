#!/usr/bin/env python3
"""
annovar_db.py
-------------
Check, repair and update an ANNOVAR humandb/ directory. Used by the pipeline
modes `--mode check_db` / `--mode update_db`, and usable on its own.

Sub-commands
  audit     For every database of a protocol: is it installed, is it intact
            (see 'verify'), is a newer version distributed by ANNOVAR, is a
            newer ClinVar release published by NCBI, and what should be done.
  verify    Integrity checks only: file present and complete (ends with a
            newline, .idx header matches the file size), plausible table
            layout, *Mrna.fa present for gene-based databases.
  download  Install official ANNOVAR databases. annotate_variation.pl downloads
            into a private staging directory; files are verified and only then
            moved into humandb/ - an interrupted download never leaves a
            half-written database behind.
  clinvar   Build <build>_clinvar_<YYYYMMDD>.txt + .idx straight from the NCBI
            ClinVar VCF (MD5-verified). NCBI publishes a release every week;
            ANNOVAR repackages one only now and then.
  install   Install every database of a protocol that is not installed yet
            (download / clinvar as above); installed ones are left as they are.
            Used before annotation (--auto_install_db) on a new humandb.
  genes     Build the gene table of the report (annovar_nf/genes.tsv: GenCC
            diseases and modes of inheritance, gnomAD constraint, HGNC symbols;
            see gene_info.py) into humandb/.

Protocol entries are '<database>[:<operation>]'. '<family>_latest' (e.g.
'clinvar_latest') means the newest installed '<family>_<YYYYMMDD>' database.

Examples
  annovar_db.py audit    --humandb humandb --build hg38 --entries refGene:g,gnomad41_genome:f,clinvar_latest:f
  annovar_db.py verify   --humandb humandb --build hg19 --entries refGene:g,avsnp151:f
  annovar_db.py download --humandb humandb --build hg38 --annovar-src /opt/annovar gnomad41_genome
  annovar_db.py clinvar  --humandb humandb --build hg38 --date latest
  annovar_db.py install  --humandb humandb --build hg38 --annovar-src /opt/annovar --entries refGene:g,clinvar_latest:f
  annovar_db.py genes    --humandb humandb --max-age 30

Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import http.client
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait as wait_futures
from dataclasses import asdict, dataclass, field
from pathlib import Path

import gene_info

ANNOVAR_DOWNLOAD = "http://www.openbioinformatics.org/annovar/download"
# Databases ANNOVAR announces as released: the table of the download page
# (annovar.openbioinformatics.org/en/latest/user-guide/download/), read from its
# source in the documentation repository. avdblist also lists unannounced or
# unfinished uploads (e.g. hg19 dbnsfp54a before its release, intervar_20250801
# without an index), which must not be suggested as upgrades.
ANNOVAR_PUBLISHED = "https://raw.githubusercontent.com/WGLab/doc-ANNOVAR/master/docs/user-guide/download.md"
NCBI_CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
NCBI_ASSEMBLY = {"hg19": "GRCh37", "hg38": "GRCh38"}

# Protocol names that are installed by a differently named ANNOVAR package.
THOUSAND_GENOMES = {"2015_08": "1000g2015aug", "2014_10": "1000g2014oct", "2014_09": "1000g2014sep",
                    "2014_08": "1000g2014aug", "2012_04": "1000g2012apr", "2012_02": "1000g2012feb"}

# ClinVar INFO fields -> ANNOVAR columns, in the layout of the clinvar_2024*
# databases (germline + oncogenicity + somatic clinical impact).
CLINVAR_COLUMNS = [
    ("CLNALLELEID", "ALLELEID"), ("CLNDN", "CLNDN"), ("CLNDISDB", "CLNDISDB"),
    ("CLNREVSTAT", "CLNREVSTAT"), ("CLNSIG", "CLNSIG"),
    ("ONCDN", "ONCDN"), ("ONCDISDB", "ONCDISDB"), ("ONCREVSTAT", "ONCREVSTAT"), ("ONC", "ONC"),
    ("SCIDN", "SCIDN"), ("SCIDISDB", "SCIDISDB"), ("SCIREVSTAT", "SCIREVSTAT"), ("SCI", "SCI"),
]
CHROM_ORDER = {c: i for i, c in enumerate([str(n) for n in range(1, 23)] + ["X", "Y", "M", "MT"])}
INDEX_BIN = 1000
BIG_FILE = 1 << 30   # databases above 1 GB should have an .idx

USE_COLOR = sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def paint(text: str, code: str) -> str:
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else text


def log(msg: str) -> None:
    print(f"[annovar_db] {msg}", file=sys.stderr, flush=True)


# ======================================================================================
# Networking
# ======================================================================================
def fetch(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as resp:  # noqa: S310 (fixed hosts)
        return resp.read()


def download_file(url: str, dest: Path, size: int | None = None, attempts: int = 30,
                  timeout: int = 120, offset: int | None = None, quiet: bool = False) -> None:
    """
    Stream *url* to *dest*, reporting progress every ~5%. Interrupted transfers
    are retried and resumed (HTTP Range) from what is already in *dest*, so an
    existing partial file - e.g. from a previous run - is continued, not
    restarted. *size*, when known, is the expected final size; with *offset*
    only bytes [offset, offset + size) of the remote file are fetched.
    """
    for attempt in range(1, attempts + 1):
        have = dest.stat().st_size if dest.exists() else 0
        if size and have == size:
            return
        if size and have > size:
            dest.unlink()
            have = 0
        if offset is not None:
            headers = {"Range": f"bytes={offset + have}-{offset + size - 1}"}
        else:
            headers = {"Range": f"bytes={have}-"} if have else {}
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as resp:  # noqa: S310 (fixed hosts)
                if headers and resp.status != 206:
                    if offset is not None:
                        raise OSError("the server does not support partial downloads")
                    have = 0                        # server ignored the range: start over
                total = have + int(resp.headers.get("Content-Length") or 0)
                total = size or total
                if have and not quiet:
                    log(f"  {dest.name}: resuming at {have / 1e6:,.0f} of {total / 1e6:,.0f} MB")
                done, start, t0 = have, have, time.monotonic()
                next_report = (done / total + 0.05) if total and not quiet else 2.0
                with open(dest, "ab" if have else "wb") as out:
                    while chunk := resp.read(1 << 20):
                        out.write(chunk)
                        done += len(chunk)
                        if total and done / total >= next_report:
                            rate = (done - start) / max(time.monotonic() - t0, 1e-6) / 1e6
                            log(f"  {dest.name}: {done / total:4.0%} of {total / 1e6:,.0f} MB · {rate:.1f} MB/s")
                            next_report += 0.05
            if not size or dest.stat().st_size == size:
                return
            raise OSError(f"connection closed at {dest.stat().st_size:,} of {size:,} bytes")
        except urllib.error.HTTPError as exc:
            if exc.code == 416 and have and not size:   # nothing left to fetch
                return
            if exc.code in (403, 404) or attempt == attempts:
                raise
            error = exc
        except (OSError, http.client.HTTPException) as exc:
            if attempt == attempts:
                raise OSError(str(exc)) from exc
            error = exc
        wait = min(60, 5 * attempt)
        log(f"  {dest.name}: interrupted ({error}); retry {attempt}/{attempts - 1} in {wait}s")
        time.sleep(wait)


def download_parallel(url: str, dest: Path, size: int, connections: int = 8,
                      part_size: int = 32 << 20, report_every: float = 30.0) -> None:
    """
    Download a large file as byte ranges of *part_size* over *connections*
    parallel connections (the ANNOVAR server throttles each connection to well
    under 1 MB/s). Small parts handed out as connections free up keep a slow
    connection from holding up the whole file. Every part is resumable on its
    own as '<dest>.partN'; the parts are joined into *dest* at the end.
    """
    if dest.exists() and dest.stat().st_size == size:
        return
    if size <= 2 * part_size:
        download_file(url, dest, size=size)
        return
    ranges = [(off, min(part_size, size - off)) for off in range(0, size, part_size)]
    paths = [dest.with_name(f"{dest.name}.part{k}") for k in range(len(ranges))]

    def downloaded() -> int:
        return sum(p.stat().st_size for p in paths if p.exists())

    log(f"  {dest.name}: {len(ranges)} parts over {connections} parallel connections")
    t0, start = time.monotonic(), downloaded()
    with ThreadPoolExecutor(max_workers=connections) as pool:
        jobs = [pool.submit(download_file, url, path, size=n, offset=off, quiet=True)
                for path, (off, n) in zip(paths, ranges)]
        while not all(j.done() for j in jobs):
            wait_futures(jobs, timeout=report_every)
            done = downloaded()
            rate = (done - start) / max(time.monotonic() - t0, 1e-6) / 1e6
            eta = (size - done) / (rate * 1e6) if rate > 0 else 0
            log(f"  {dest.name}: {done / size:4.0%} of {size / 1e6:,.0f} MB · {rate:.1f} MB/s"
                + (f" · ~{eta / 60:.0f} min left" if done < size and eta else ""))
        for j in jobs:
            j.result()                  # re-raise the first failure (finished parts are kept)
    with open(dest, "wb") as out:
        for p in paths:
            with open(p, "rb") as src:
                shutil.copyfileobj(src, out, 1 << 24)
    if dest.stat().st_size != size:
        dest.unlink()
        raise OSError(f"joined file has the wrong size (expected {size:,} bytes)")
    for p in paths:
        p.unlink()


def avdblist_text(build: str) -> str:
    return gzip.decompress(fetch(f"{ANNOVAR_DOWNLOAD}/{build}_avdblist.txt.gz")).decode("utf-8", "replace")


def fetch_avdblist(build: str) -> dict[str, str]:
    return parse_avdblist(avdblist_text(build), build)


def parse_published(page: str, build: str) -> set[str]:
    """
    Database names listed for *build* in the markdown tables of the ANNOVAR
    download page: rows '| build | name | explanation | date |'. Builds like
    'hg19/hg38' count for each; '1000g2015aug (6 data sets)' -> '1000g2015aug'.
    """
    names = set()
    for line in page.splitlines():
        cells = [c.strip().strip("`*") for c in line.strip().strip("|").split("|")]
        if not line.lstrip().startswith("|") or len(cells) < 2 or build not in cells[0].split("/"):
            continue
        name = cells[1].split()[0].strip("`*") if cells[1].split() else ""
        if re.fullmatch(r"[\w.+-]+", name):
            names.add(name)
    return names


def fetch_published(build: str) -> set[str]:
    return parse_published(fetch(ANNOVAR_PUBLISHED).decode("utf-8", "replace"), build)


def latest_ncbi_clinvar(build: str) -> str | None:
    """Newest clinvar_YYYYMMDD.vcf.gz date listed by NCBI for the build."""
    listing = fetch(f"{NCBI_CLINVAR}/vcf_{NCBI_ASSEMBLY[build]}/").decode("utf-8", "replace")
    dates = re.findall(r"clinvar_(\d{8})\.vcf\.gz\"", listing)
    return max(dates) if dates else None


# ======================================================================================
# Names and versions
# ======================================================================================
def parse_entries(text: str) -> list[tuple[str, str]]:
    """'refGene:g,avsnp151' -> [('refGene', 'g'), ('avsnp151', 'f')] (gene models default to g)."""
    out = []
    for item in filter(None, (e.strip() for e in text.split(","))):
        name, _, op = item.partition(":")
        if not op:
            op = "g" if re.fullmatch(r"(?i)refGene\w*|knownGene|ensGene\w*|ccdsGene|wgEncodeGencode\w+", name) else "f"
        out.append((name, op))
    return out


def family_and_version(name: str) -> tuple[str, tuple[int, ...]]:
    """'gnomad41_exome' -> ('gnomad#_exome', (41,)); 'clinvar_20240917' -> ('clinvar_#', (20240917,))."""
    return re.sub(r"\d+", "#", name), tuple(int(n) for n in re.findall(r"\d+", name))


def parse_avdblist(text: str, build: str) -> dict[str, str]:
    """avdblist lines '<build>_<name>.txt.gz|.zip<TAB>date<TAB>size' -> {name: date}."""
    out = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 2 or ".idx." in parts[0]:
            continue
        m = re.match(rf"^{build}_(.+?)\.(txt\.gz|zip)$", parts[0])
        if m:
            out[m.group(1)] = parts[1]
    return out


def package_files(text: str, build: str, package: str) -> dict[str, int]:
    """
    Files of one ANNOVAR package listed in avdblist, with their sizes:
    '<build>_<package>.txt.gz' plus, when present, '.txt.idx.gz', 'Mrna.fa.gz'
    and 'Version.txt.gz'. Empty when the package is not listed or ships as a
    .zip (unpacked by annotate_variation.pl instead).
    """
    files = {}
    pattern = re.compile(rf"^{re.escape(build)}_{re.escape(package)}"
                         r"(\.txt\.gz|\.txt\.idx\.gz|Mrna\.fa\.gz|Version\.txt\.gz|\.zip)$")
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and pattern.match(parts[0]) and parts[2].strip().isdigit():
            files[parts[0]] = int(parts[2])
    if any(name.endswith(".zip") for name in files) or f"{build}_{package}.txt.gz" not in files:
        return {}
    return files


def download_name(db: str) -> str:
    """ANNOVAR package that installs protocol database *db*."""
    m = re.fullmatch(r"(?:ALL|AFR|AMR|EAS|EUR|SAS)\.sites\.(\d{4}_\d{2})", db)
    if m and m.group(1) in THOUSAND_GENOMES:
        return THOUSAND_GENOMES[m.group(1)]
    return db


def release_key(name: str, dates: dict[str, str] | None = None) -> tuple:
    """
    Sort key for versions of one family. Version numbers are not comparable
    across schemes (gnomad211 = 2.1.1 is older than gnomad41 = 4.1), so the
    ANNOVAR release date is used first when known, then the numbers.
    """
    dates = dates or {}
    # 1000 Genomes tables (ALL.sites.2015_08 ...) carry the date of their package (1000g2015aug).
    date = dates.get(name) or dates.get(download_name(name), "")
    return (int(date) if date.isdigit() else 0, family_and_version(name)[1])


def newest_in_family(name: str, candidates, dates: dict[str, str] | None = None) -> str | None:
    family, _ = family_and_version(name)
    same = [c for c in candidates if family_and_version(c)[0] == family]
    return max(same, key=lambda c: release_key(c, dates)) if same else None


def installed_versions(humandb: Path, build: str, family: str) -> list[str]:
    """Installed '<family>_<digits>' databases, oldest first."""
    found = []
    for p in humandb.glob(f"{build}_{family}_*.txt"):
        name = p.name[len(build) + 1:-4]
        if re.fullmatch(rf"{re.escape(family)}_\d+", name):
            found.append(name)
    return sorted(found, key=lambda n: family_and_version(n)[1])


def resolve(humandb: Path, build: str, name: str) -> str | None:
    """Resolve '<family>_latest' to the newest installed version (None if none)."""
    if not name.endswith("_latest"):
        return name
    versions = installed_versions(humandb, build, name[: -len("_latest")])
    return versions[-1] if versions else None


# ======================================================================================
# verify
# ======================================================================================
@dataclass
class Check:
    level: str          # ok | warning | error
    message: str


def last_byte(path: Path) -> bytes:
    with open(path, "rb") as fh:
        fh.seek(-1, os.SEEK_END)
        return fh.read(1)


def first_data_line(path: Path) -> str | None:
    with open(path, encoding="utf-8", errors="replace") as fh:
        for _ in range(50):
            line = fh.readline()
            if not line:
                return None
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\n")
    return None


def verify_database(humandb: Path, build: str, db: str, op: str) -> list[Check]:
    """Integrity checks for one database; an empty error list means usable."""
    txt = humandb / f"{build}_{db}.txt"
    if not txt.exists():
        return [Check("error", "not installed")]
    size = txt.stat().st_size
    if size == 0:
        return [Check("error", "file is empty")]
    checks = []
    if last_byte(txt) != b"\n":
        checks.append(Check("error", "file does not end with a newline (truncated download?)"))

    line = first_data_line(txt)
    n_cols = len(line.split("\t")) if line else 0
    if op.startswith("g"):
        if n_cols < 10:
            checks.append(Check("error", f"expected a genePred table (>=10 columns), found {n_cols}"))
        if not (humandb / f"{build}_{db}Mrna.fa").exists():
            checks.append(Check("error", f"{build}_{db}Mrna.fa missing (needed for exonic consequences)"))
    elif n_cols < 5 and op.startswith("f"):
        checks.append(Check("error", f"expected >=5 tab-separated columns, found {n_cols}"))

    idx = Path(f"{txt}.idx")
    if idx.exists():
        with open(idx, encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
        m = re.match(r"#BIN\t(\d+)\t(\d+)", header)
        if not m:
            checks.append(Check("warning", ".idx header is malformed; ANNOVAR will scan the whole file"))
        elif int(m.group(2)) != size:
            checks.append(Check("error", f".idx expects {int(m.group(2)):,} bytes but the file has {size:,} "
                                         "(incomplete or modified database)"))
    elif size > BIG_FILE and op.startswith("f"):
        checks.append(Check("warning", f"no .idx for a {size / 1e9:.1f} GB database: annotation will be slow"))

    for leftover in (f"{txt}.gz", f"{txt}.idx.gz"):
        if Path(leftover).exists():
            checks.append(Check("warning", f"leftover {Path(leftover).name} (interrupted download?)"))
    return checks or [Check("ok", "")]


def summarize_checks(checks: list[Check]) -> tuple[str, str]:
    for level in ("error", "warning"):
        msgs = [c.message for c in checks if c.level == level]
        if msgs:
            return level, "; ".join(msgs)
    return "ok", ""


# ======================================================================================
# audit
# ======================================================================================
@dataclass
class AuditRow:
    entry: str
    database: str
    operation: str
    installed: bool
    size_gb: float
    integrity: str
    newest_official: str
    status: str
    action: str           # none | download | redownload | clinvar | manual
    package: str = ""
    notes: list[str] = field(default_factory=list)


def audit(humandb: Path, build: str, entries: list[tuple[str, str]], avdb: dict[str, str] | None,
          ncbi_clinvar: str | None, published: set[str] | None = None) -> list[AuditRow]:
    """
    *avdb*: downloadable packages {name: date} (None when offline).
    *published*: names on the ANNOVAR download page; when given, only these are
    suggested as newer versions (None: every avdblist package counts).
    """
    rows = []
    releases = {n: d for n, d in (avdb or {}).items() if published is None or n in published}
    for entry, op in entries:
        db = resolve(humandb, build, entry)
        is_clinvar = entry.startswith("clinvar_")
        if db is None:  # '<family>_latest' with nothing installed
            action = "clinvar" if is_clinvar else "manual"
            rows.append(AuditRow(entry, entry, op, False, 0.0, "error", "-", "missing", action,
                                 notes=["no installed version"]))
            continue
        checks = verify_database(humandb, build, db, op)
        integrity, message = summarize_checks(checks)
        txt = humandb / f"{build}_{db}.txt"
        installed = txt.exists()
        size_gb = round(txt.stat().st_size / 1e9, 2) if installed else 0.0
        package = download_name(db)
        # Offline (avdb is None): availability unknown - assume official and let the download decide.
        official = avdb is None or package in avdb
        # Tables installed from a differently named package count as officially available.
        candidates = list(releases) + ([db] if official else [])
        newest = (newest_in_family(db, candidates, avdb) or "-") if avdb else "?"
        if is_clinvar:
            # ClinVar versions are NCBI release dates and are built from NCBI, so the
            # NCBI release - not ANNOVAR's older repackaged copies - is the reference.
            newest = f"clinvar_{ncbi_clinvar}" if ncbi_clinvar else newest

        row = AuditRow(entry, db, op, installed, size_gb, integrity, newest, "ok", "none", package)
        if message:
            row.notes.append(message)
        if not installed:
            row.status = "missing"
            row.action = "clinvar" if is_clinvar else ("download" if official else "manual")
        elif integrity == "error":
            row.status = "damaged"
            row.action = "clinvar" if is_clinvar else ("redownload" if official else "manual")
        elif is_clinvar and ncbi_clinvar and family_and_version(db)[1][0] < int(ncbi_clinvar):
            row.status = f"outdated (NCBI clinvar_{ncbi_clinvar})"
            row.action = "clinvar" if entry.endswith("_latest") else "none"
            if row.action == "none":
                row.notes.append("pinned version; use 'clinvar_latest' in the protocol to follow NCBI releases")
        elif is_clinvar:
            if not official:
                row.notes.append("built from NCBI ClinVar")
        elif avdb and newest not in ("-", db) and release_key(newest, avdb) > release_key(db, avdb):
            row.status = f"newer available: {newest}"
            row.notes.append("change the protocol to switch versions")
        elif avdb is not None and not official:
            row.status = "custom build"
        if integrity == "warning" and row.status == "ok":
            row.status = "ok (warning)"
        if avdb is not None and not official and row.action == "manual":
            row.notes.append("not distributed by ANNOVAR for this build; install it manually")
        rows.append(row)
    return rows


def print_audit(rows: list[AuditRow], build: str, humandb: Path, ncbi: str | None) -> None:
    cols = [("database", lambda r: r.database), ("size", lambda r: f"{r.size_gb:.2f} GB" if r.installed else "-"),
            ("newest official", lambda r: r.newest_official), ("status", lambda r: r.status),
            ("action", lambda r: r.action)]
    widths = [max(len(name), *(len(get(r)) for r in rows)) for name, get in cols]
    colors = {"ok": "32", "missing": "31", "damaged": "31", "custom build": "36"}
    print(paint(f"\n  ANNOVAR databases · {build} · {humandb}", "1"))
    print(paint("  " + "  ".join(n.ljust(w) for (n, _), w in zip(cols, widths)), "2"))
    for r in rows:
        cells = [get(r).ljust(w) for (_, get), w in zip(cols, widths)]
        code = colors.get(r.status, "33")
        cells[3] = paint(cells[3], code)
        print("  " + "  ".join(cells))
        for note in r.notes:
            print(paint(f"      ↳ {note}", "2"))
    todo = [r for r in rows if r.action != "none"]
    print(f"\n  {len(rows)} databases · {paint(str(len(todo)) + ' need action', '33' if todo else '32')}"
          + (f" · newest ClinVar at NCBI: clinvar_{ncbi}" if ncbi else "") + "\n")


def cmd_audit(args) -> int:
    humandb = Path(args.humandb)
    entries = parse_entries(args.entries)
    avdb = None
    if not args.offline:
        try:
            log(f"fetching the ANNOVAR database list for {args.build}")
            avdb = fetch_avdblist(args.build)
        except OSError as exc:
            log(f"WARNING: ANNOVAR server unreachable ({exc}); version checks skipped")
    published = None
    if avdb is not None:
        try:
            published = fetch_published(args.build) or None
        except OSError as exc:
            log(f"WARNING: ANNOVAR download page unreachable ({exc}); newer versions taken from avdblist")
    ncbi = None
    if not args.offline and any(e.startswith("clinvar_") for e, _ in entries):
        try:
            ncbi = latest_ncbi_clinvar(args.build)
        except OSError as exc:
            log(f"WARNING: NCBI unreachable ({exc}); ClinVar release check skipped")

    rows = audit(humandb, args.build, entries, avdb, ncbi, published)
    print_audit(rows, args.build, humandb, ncbi)
    if args.tsv:
        with open(args.tsv, "w", encoding="utf-8") as fh:
            fh.write("build\tentry\tdatabase\toperation\tinstalled\tsize_gb\tintegrity\tnewest_official\t"
                     "status\taction\tpackage\tnotes\n")
            for r in rows:
                fh.write("\t".join(map(str, [args.build, r.entry, r.database, r.operation, r.installed, r.size_gb,
                                             r.integrity, r.newest_official, r.status, r.action, r.package,
                                             " | ".join(r.notes)])) + "\n")
    if args.json:
        Path(args.json).write_text(json.dumps({"build": args.build, "humandb": str(humandb),
                                               "ncbi_clinvar": ncbi, "databases": [asdict(r) for r in rows]},
                                              indent=2) + "\n")
    return 1 if args.strict and any(r.action != "none" for r in rows) else 0


def cmd_verify(args) -> int:
    humandb = Path(args.humandb)
    bad = 0
    for entry, op in parse_entries(args.entries):
        db = resolve(humandb, args.build, entry)
        checks = verify_database(humandb, args.build, db, op) if db else [Check("error", "no installed version")]
        level, msg = summarize_checks(checks)
        bad += level == "error"
        mark = {"ok": paint("✔", "32"), "warning": paint("!", "33"), "error": paint("✘", "31")}[level]
        print(f"  {mark} {args.build}_{db or entry}  {paint(msg, '2')}")
    return 1 if bad else 0


# ======================================================================================
# download
# ======================================================================================
def fetch_package(staging: Path, files: dict[str, int]) -> bool:
    """Download (resumable, size-checked) and unpack the .gz files of one package."""
    for name, size in sorted(files.items(), key=lambda item: item[1]):
        log(f"  {name}: {size / 1e6:,.0f} MB")
        try:
            download_parallel(f"{ANNOVAR_DOWNLOAD}/{name}", staging / name, size=size)
        except OSError as exc:
            log(paint(f"  {name}: download failed: {exc}", "31"))
            log("  the partial file is kept; rerun to resume the download")
            return False
    log("  unpacking")
    for name in files:
        archive = staging / name
        try:
            with gzip.open(archive, "rb") as src, open(staging / name[:-3], "wb") as out:
                shutil.copyfileobj(src, out, 1 << 24)
        except (OSError, EOFError) as exc:
            log(paint(f"  {name}: corrupt archive ({exc}); deleted, rerun to download it again", "31"))
            archive.unlink(missing_ok=True)
            return False
    return True


def run_annotate_variation(staging: Path, build: str, package: str, annovar_src: Path) -> bool:
    """Fallback for packages not listed as plain .gz files (e.g. 1000 Genomes .zip)."""
    cmd = ["perl", str(annovar_src / "annotate_variation.pl"), "-buildver", build,
           "-downdb", "-webfrom", "annovar", package, str(staging)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    output = (proc.stdout + proc.stderr).strip().splitlines()
    if proc.returncode != 0 or any(re.search(r"\bfailed\b|^error", line, re.I) for line in output):
        log(paint(f"{build}_{package}: annotate_variation.pl failed:", "31"))
        log("\n".join(output[-12:]))
        return False
    return True


def install_package(humandb: Path, build: str, package: str, annovar_src: Path, force: bool,
                    files: dict[str, int] | None = None) -> bool:
    """
    Download one ANNOVAR package into a staging dir, verify, then move it into
    humandb. The staging dir keeps downloaded archives when something fails,
    so the next run resumes instead of starting over; it is removed on success.
    """
    staging = humandb / f".annovar-nf-staging.{build}_{package}"
    staging.mkdir(exist_ok=True)
    for p in staging.iterdir():     # keep (partial) archives, drop anything unpacked
        if p.is_dir():
            shutil.rmtree(p)
        elif not (files and p.name.split(".part")[0] in files):
            p.unlink()

    log(f"{build}_{package}: downloading" + (" (via annotate_variation.pl)" if not files else ""))
    ok = fetch_package(staging, files) if files else run_annotate_variation(staging, build, package, annovar_src)
    if not ok:
        return False

    unpacked = [p for p in staging.iterdir() if p.is_file() and p.name.startswith(f"{build}_")
                and re.search(r"\.(txt|txt\.idx|fa)$", p.name)]
    tables = [p for p in unpacked if p.name.endswith(".txt")]
    if not tables:
        log(paint(f"{build}_{package}: nothing was downloaded", "31"))
        return False
    for t in tables:
        db = t.name[len(build) + 1:-4]
        if db.endswith(("Mrna", "Version")):
            continue
        level, msg = summarize_checks(verify_database(staging, build, db, "f"))
        if level == "error":
            log(paint(f"{t.name}: verification failed: {msg}", "31"))
            return False
    # .idx before the table: ANNOVAR ignores an index whose recorded size does not match.
    for p in sorted(unpacked, key=lambda q: not q.name.endswith(".idx")):
        target = humandb / p.name
        if target.exists() and not force:
            log(f"  {p.name}: already present, keeping the installed copy")
            continue
        os.replace(p, target)
        log(f"  installed {p.name} ({target.stat().st_size / 1e6:,.0f} MB)")
    shutil.rmtree(staging, ignore_errors=True)
    return True


def cmd_download(args) -> int:
    humandb, src = Path(args.humandb), Path(args.annovar_src)
    if not (src / "annotate_variation.pl").exists():
        raise SystemExit(f"annotate_variation.pl not found in {src}")
    try:
        listing = avdblist_text(args.build)
    except OSError as exc:
        log(f"cannot read the ANNOVAR file list ({exc}); falling back to annotate_variation.pl")
        listing = ""
    failed = []
    for db in args.databases:
        package = download_name(db)
        if (humandb / f"{args.build}_{db}.txt").exists() and not args.force:
            log(f"{args.build}_{db}: already installed (use --force to replace)")
            continue
        files = package_files(listing, args.build, package)
        if not install_package(humandb, args.build, package, src, args.force, files):
            failed.append(db)
    if failed:
        log(paint(f"failed: {', '.join(failed)}", "31"))
        return 1
    log(paint("all requested databases are installed and verified", "32"))
    return 0


# ======================================================================================
# clinvar
# ======================================================================================
def vcf_to_annovar(pos: int, ref: str, alt: str) -> tuple[int, int, str, str]:
    """
    Convert a VCF allele to ANNOVAR coordinates like convert2annovar.pl: strip
    the shared leading bases, '-' marks an empty allele, insertions sit on the
    base preceding the inserted sequence.
    """
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    ref2, alt2 = ref[k:], alt[k:]
    if not ref2:                                   # insertion
        start = pos + k - 1
        return start, start, "-", alt2
    start = pos + k
    return start, start + len(ref2) - 1, ref2, alt2 or "-"


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for item in info.split(";"):
        key, _, value = item.partition("=")
        out[key] = value
    return out


def clinvar_records(vcf_lines):
    """Yield (chrom, start, line) for every ClinVar allele."""
    for line in vcf_lines:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t", 8)
        if len(f) < 8:
            continue
        chrom, pos, _vid, ref, alts, _q, _flt, info = f[:8]
        chrom = re.sub(r"^chr", "", chrom)
        if chrom not in CHROM_ORDER:
            continue
        info_d = parse_info(info)
        # ANNOVAR splits filter annotations on commas: escape them like ANNOVAR does.
        values = [(info_d.get(src) or ".").replace(",", "\\x2c") for _, src in CLINVAR_COLUMNS]
        for alt in alts.split(","):
            if alt in (".", "*") or alt.startswith("<"):
                continue
            start, end, r, a = vcf_to_annovar(int(pos), ref, alt)
            fields = "\t".join([str(start), str(end), r, a, *values])
            yield chrom, start, f"{chrom}\t{fields}\n"
            # UCSC-style VCFs call the mitochondrion chrM ('M' in ANNOVAR), NCBI uses MT: provide both.
            if chrom == "MT":
                yield "M", start, f"M\t{fields}\n"


def write_database(records, out_txt: Path) -> int:
    """Sort records, write the ANNOVAR table and its .idx; return the record count."""
    records = sorted(records, key=lambda r: (CHROM_ORDER[r[0]], r[1]))
    header = "#Chr\tStart\tEnd\tRef\tAlt\t" + "\t".join(c for c, _ in CLINVAR_COLUMNS) + "\n"
    bins: list[tuple[str, int, int, int]] = []
    offset = 0
    with open(out_txt, "wb") as out:
        data = header.encode()
        out.write(data)
        offset += len(data)
        for chrom, start, line in records:
            data = line.encode()
            b = start - start % INDEX_BIN
            if bins and bins[-1][0] == chrom and bins[-1][1] == b:
                bins[-1] = (chrom, b, bins[-1][2], offset + len(data))
            else:
                bins.append((chrom, b, offset, offset + len(data)))
            out.write(data)
            offset += len(data)
    with open(f"{out_txt}.idx", "w", encoding="utf-8") as idx:
        idx.write(f"#BIN\t{INDEX_BIN}\t{offset}\n")
        for chrom, b, o0, o1 in bins:
            idx.write(f"{chrom}\t{b}\t{o0}\t{o1}\n")
    return len(records)


def clinvar_urls(build: str, date: str) -> list[str]:
    base = f"{NCBI_CLINVAR}/vcf_{NCBI_ASSEMBLY[build]}"
    return [f"{base}/clinvar_{date}.vcf.gz", f"{base}/archive_2.0/{date[:4]}/clinvar_{date}.vcf.gz"]


def build_clinvar(humandb: Path, build: str, release: str, vcf_path: str | None = None,
                  tmpdir: str | None = None, force: bool = False) -> str:
    """Build <build>_clinvar_<date> (release: YYYYMMDD or 'latest') unless installed; returns its name."""
    date = latest_ncbi_clinvar(build) if release == "latest" else release
    if not date or not re.fullmatch(r"\d{8}", date):
        raise SystemExit(f"invalid or unavailable ClinVar date: {release}")
    name = f"clinvar_{date}"
    target = humandb / f"{build}_{name}.txt"
    if target.exists() and not force:
        log(f"{target.name} already installed - nothing to do")
        return name

    with tempfile.TemporaryDirectory(prefix=".annovar-nf-staging.", dir=tmpdir or humandb) as tmp:
        tmp = Path(tmp)
        vcf = tmp / f"clinvar_{date}.vcf.gz"
        if vcf_path:
            shutil.copy(vcf_path, vcf)
        else:
            for url in clinvar_urls(build, date):
                try:
                    log(f"downloading {url}")
                    download_file(url, vcf)
                    expected = fetch(url + ".md5").decode().split()[0]
                    break
                except OSError as exc:
                    log(f"  not available here ({exc})")
            else:
                raise SystemExit(f"ClinVar {date} not found for {build} at NCBI")
            digest = hashlib.md5(vcf.read_bytes()).hexdigest()  # noqa: S324 (integrity check only)
            if digest != expected:
                raise SystemExit(f"MD5 mismatch for {vcf.name}: {digest} != {expected}")
            log("MD5 verified")

        log(f"converting to {target.name}")
        tmp_txt = tmp / target.name
        with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as fh:
            n = write_database(clinvar_records(fh), tmp_txt)
        level, msg = summarize_checks(verify_database(tmp, build, name, "f"))
        if level == "error":
            raise SystemExit(f"built database failed verification: {msg}")
        # .idx first: ANNOVAR ignores an index whose recorded size does not match.
        os.replace(f"{tmp_txt}.idx", f"{target}.idx")
        os.replace(tmp_txt, target)

    log(paint(f"installed {target.name} ({n:,} alleles) + .idx", "32"))
    return name


def cmd_clinvar(args) -> int:
    print(build_clinvar(Path(args.humandb), args.build, args.date, args.vcf, args.tmpdir, args.force))
    return 0


# ======================================================================================
# install
# ======================================================================================
def missing_entries(humandb: Path, build: str, entries: list[tuple[str, str]]) -> list[str]:
    """Protocol entries without an installed table ('<family>_latest': no version installed)."""
    missing = []
    for entry, _op in entries:
        db = resolve(humandb, build, entry)
        if db is None or not (humandb / f"{build}_{db}.txt").exists():
            missing.append(entry)
    return list(dict.fromkeys(missing))


def install_missing(humandb: Path, build: str, entries: list[tuple[str, str]], annovar_src: Path,
                    clinvar_release: str = "latest", forks: int = 2) -> list[str]:
    """
    Install the protocol entries that are not installed: ANNOVAR packages (several at
    a time) and ClinVar from NCBI. Installed databases are left as they are.
    Returns the entries that could not be installed.
    """
    missing = missing_entries(humandb, build, entries)
    if not missing:
        log(f"{build}: every protocol database is installed")
        return []
    clinvar = [e for e in missing if e.startswith("clinvar_")]
    try:
        listing = avdblist_text(build)
    except OSError as exc:
        log(f"cannot read the ANNOVAR file list ({exc}); falling back to annotate_variation.pl")
        listing = ""
    available = parse_avdblist(listing, build) if listing else None
    packages: dict[str, list[str]] = {}       # package -> the protocol entries it installs
    manual = []
    for e in missing:
        if e in clinvar:
            continue
        package = download_name(e)
        if available is not None and package not in available:
            manual.append(e)
        else:
            packages.setdefault(package, []).append(e)
    files = {p: package_files(listing, build, p) for p in packages}
    size = sum(sum(f.values()) for f in files.values())
    log(f"{build}: installing {len(missing) - len(manual)} database(s) into {humandb}"
        + (f" ({size / 1e9:,.1f} GB to download)" if size else "") + ": "
        + ", ".join(e for e in missing if e not in manual))

    failed = list(manual)
    with ThreadPoolExecutor(max_workers=max(1, forks)) as pool:
        jobs = {p: pool.submit(install_package, humandb, build, p, annovar_src, False, files[p] or None)
                for p in packages}
        for p, job in jobs.items():
            try:
                ok = job.result()
            except OSError as exc:
                log(paint(f"{build}_{p}: {exc}", "31"))
                ok = False
            if not ok:
                failed += packages[p]
    for e in clinvar:
        try:
            build_clinvar(humandb, build, clinvar_release if e.endswith("_latest") else e[len("clinvar_"):])
        except SystemExit as exc:
            log(paint(f"{e}: {exc}", "31"))
            failed.append(e)
    if manual:
        log(paint(f"not distributed by ANNOVAR for {build}; install them yourself or change the protocol: "
                  f"{', '.join(manual)}", "31"))
    return failed


def cmd_install(args) -> int:
    import fcntl  # POSIX only, like ANNOVAR itself

    humandb = Path(args.humandb)
    humandb.mkdir(parents=True, exist_ok=True)
    if not (Path(args.annovar_src) / "annotate_variation.pl").exists():
        raise SystemExit(f"annotate_variation.pl not found in {args.annovar_src}")
    # Two runs that start on the same new humandb: the second waits, then finds the databases installed.
    with open(humandb / ".annovar-nf-install.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        failed = install_missing(humandb, args.build, parse_entries(args.entries), Path(args.annovar_src),
                                 args.clinvar_release, args.forks)
    if failed:
        log(paint(f"not installed: {', '.join(failed)}", "31"))
        return 1
    log(paint(f"{args.build}: every protocol database is installed and verified", "32"))
    return 0


# ======================================================================================
# genes
# ======================================================================================
GENE_SOURCES = {
    "gencc": "https://search.thegencc.org/download/action/submissions-export-tsv",
    "constraint": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
                  "gnomad.v4.1.constraint_metrics.tsv",
    "hgnc": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
}


def remote_info(url: str) -> tuple[int | None, str]:
    """Size and last-modified date (YYYY-MM-DD) of a remote file, as far as the server tells."""
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as resp:  # noqa: S310 (fixed hosts)
        size = resp.headers.get("Content-Length")
        modified = resp.headers.get("Last-Modified") or ""
    try:
        date = time.strftime("%Y-%m-%d", time.strptime(modified, "%a, %d %b %Y %H:%M:%S GMT"))
    except ValueError:
        date = ""
    return (int(size) if size and size.isdigit() else None), date


def table_age_days(path: Path) -> float:
    return (time.time() - path.stat().st_mtime) / 86400


def build_genes(humandb: Path, files: dict[str, str | None] | None = None, max_age: float | None = None,
                tmpdir: str | None = None, attempts: int = 5) -> Path:
    """
    Build humandb/annovar_nf/genes.tsv from GenCC, the gnomAD v4.1 constraint metrics and the
    HGNC complete set (downloaded, or the local *files*), unless it exists and is younger than
    *max_age* days (None: keep any existing table). Installed atomically.
    """
    target = humandb / gene_info.TABLE_NAME
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and (max_age is None or table_age_days(target) < max_age):
        log(f"{target.name} already installed ({table_age_days(target):.0f} days old) - nothing to do")
        return target
    files = files or {}
    with tempfile.TemporaryDirectory(prefix=".annovar-nf-staging.", dir=tmpdir or target.parent) as tmp:
        local, dates = {}, {}
        for key, url in GENE_SOURCES.items():
            if files.get(key):
                local[key], dates[key] = Path(files[key]), "local file"
                continue
            local[key] = Path(tmp) / f"{key}.tsv"
            try:
                size, dates[key] = remote_info(url)
            except OSError:
                size, dates[key] = None, ""
            log(f"downloading {url}" + (f" ({size / 1e6:,.0f} MB)" if size else ""))
            if size and size > 16 << 20:            # GenCC's server is slow per connection
                download_parallel(url, local[key], size, part_size=8 << 20)
            else:
                download_file(url, local[key], size=size, attempts=attempts)

        def read(key, parser):
            with open(local[key], encoding="utf-8", errors="replace", newline="") as fh:
                data = parser(fh)
            if not data:
                raise SystemExit(f"{key}: no usable rows in {local[key]} (format changed?)")
            return data

        gencc = read("gencc", gene_info.read_gencc)
        constraint = read("constraint", gene_info.read_constraint)
        hgnc = read("hgnc", gene_info.read_hgnc)
        rows = gene_info.build_rows(gencc, constraint, hgnc)
        stamp = time.strftime("%Y-%m-%d")
        sources = [f"annovar-nf gene table, built {stamp} by annovar_db.py genes",
                   f"GenCC gene-disease validity (CC0), {dates['gencc'] or 'download ' + stamp}: https://thegencc.org",
                   "gnomAD v4.1 constraint metrics (CC0): https://gnomad.broadinstitute.org",
                   f"HGNC approved symbols, {dates['hgnc'] or 'download ' + stamp}: https://www.genenames.org"]
        staged = Path(tmp) / target.name
        gene_info.write_table(rows, staged, sources)
        os.replace(staged, target)
    with_disease = sum(1 for r in rows if r[4] != ".")
    log(paint(f"installed {target.name}: {len(rows):,} genes, {with_disease:,} with GenCC diseases", "32"))
    return target


def cmd_genes(args) -> int:
    humandb = Path(args.humandb)
    humandb.mkdir(parents=True, exist_ok=True)
    files = {"gencc": args.gencc, "constraint": args.constraint, "hgnc": args.hgnc}
    try:
        build_genes(humandb, files, args.max_age, args.tmpdir, attempts=3 if args.optional else 5)
    except (OSError, SystemExit) as exc:
        if not args.optional:
            raise
        log(paint(f"gene table not installed ({exc}); the report shows no gene information", "33"))
    return 0


# ======================================================================================
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, entries=True):
        sp.add_argument("--humandb", required=True, help="ANNOVAR humandb directory")
        sp.add_argument("--build", choices=["hg19", "hg38"], required=True)
        if entries:
            sp.add_argument("--entries", required=True, help="protocol: db[:op],db[:op],...")

    a = sub.add_parser("audit", help="status and recommended actions for a protocol")
    common(a)
    a.add_argument("--tsv", help="write the audit as TSV")
    a.add_argument("--json", help="write the audit as JSON")
    a.add_argument("--offline", action="store_true", help="skip ANNOVAR/NCBI version checks")
    a.add_argument("--strict", action="store_true", help="exit 1 when any action is needed")
    a.set_defaults(func=cmd_audit)

    v = sub.add_parser("verify", help="integrity checks only (no network)")
    common(v)
    v.set_defaults(func=cmd_verify)

    d = sub.add_parser("download", help="install official ANNOVAR databases (verified, atomic)")
    common(d, entries=False)
    d.add_argument("--annovar-src", required=True, help="directory with annotate_variation.pl")
    d.add_argument("--force", action="store_true", help="replace installed files")
    d.add_argument("databases", nargs="+", help="protocol database names (e.g. ALL.sites.2015_08, gnomad41_genome)")
    d.set_defaults(func=cmd_download)

    c = sub.add_parser("clinvar", help="build <build>_clinvar_<date> from the NCBI ClinVar VCF")
    common(c, entries=False)
    c.add_argument("--date", default="latest", help="YYYYMMDD or 'latest' [latest]")
    c.add_argument("--vcf", help="use a local ClinVar VCF instead of downloading")
    c.add_argument("--tmpdir", help="scratch directory [inside humandb]")
    c.add_argument("--force", action="store_true", help="rebuild even if installed")
    c.set_defaults(func=cmd_clinvar)

    i = sub.add_parser("install", help="install the protocol databases that are not installed yet")
    common(i)
    i.add_argument("--annovar-src", required=True, help="directory with annotate_variation.pl")
    i.add_argument("--clinvar-release", default="latest",
                   help="ClinVar built for 'clinvar_latest' when no release is installed: latest | YYYYMMDD [latest]")
    i.add_argument("--forks", type=int, default=2, help="packages downloaded at the same time [2]")
    i.set_defaults(func=cmd_install)

    g = sub.add_parser("genes", help="build the report's gene table (GenCC, gnomAD constraint, HGNC)")
    g.add_argument("--humandb", required=True, help="ANNOVAR humandb directory")
    g.add_argument("--max-age", type=float, help="rebuild a table older than this many days [keep it]")
    g.add_argument("--gencc", help="local GenCC submissions TSV instead of downloading")
    g.add_argument("--constraint", help="local gnomAD constraint metrics TSV instead of downloading")
    g.add_argument("--hgnc", help="local HGNC complete set instead of downloading")
    g.add_argument("--tmpdir", help="scratch directory [next to the table]")
    g.add_argument("--optional", action="store_true", help="a failed download is a warning, not an error")
    g.set_defaults(func=cmd_genes)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
