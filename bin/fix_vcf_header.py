#!/data/miniconda3/envs/wgs_env/bin/python3
"""
fix_vcf_header.py
-----------------
Replace ANNOVAR's generic 'Otherinfo' columns in the .txt output with
the real VCF column names (#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE…)
extracted from the companion .vcf output file.

Usage:
    fix_vcf_header.py --vcf sample.hg38_multianno.vcf \
                      --txt sample.hg38_multianno.txt \
                      --out sample.annotated.header.txt
"""

import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Replace ANNOVAR Otherinfo header with real VCF column names"
    )
    parser.add_argument("--vcf", required=True,
                        help="ANNOVAR output VCF file (*.multianno.vcf)")
    parser.add_argument("--txt", required=True,
                        help="ANNOVAR output TXT file (*.multianno.txt)")
    parser.add_argument("--out", required=True,
                        help="Output file path for fixed TXT")
    return parser.parse_args()


def extract_header(filepath: str, marker: str) -> list[str]:
    """Return the first tab-split line that contains *marker*."""
    with open(filepath) as fh:
        for line in fh:
            if marker in line:
                return line.rstrip("\n").split("\t")
    raise ValueError(f"Header with marker '{marker}' not found in {filepath}")


def count_columns(filepath: str) -> int:
    """Count columns in the first non-comment line."""
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith("##"):
                return len(line.rstrip("\n").split("\t"))
    return 0


def build_new_header(txt_header: list[str], vcf_header: list[str]) -> list[str]:
    """
    Replace the trailing Otherinfo columns in txt_header with the real VCF columns.
    The two headers should share the same total column count.
    """
    n_vcf = len(vcf_header)
    anno_cols = txt_header[:-n_vcf]   # ANNOVAR annotation columns (left side)
    return anno_cols + vcf_header


def write_fixed_file(txt_path: str, new_header: list[str], out_path: str) -> int:
    """Write fixed file; return number of data rows written."""
    rows_written = 0
    with open(out_path, "w") as out_fh:
        out_fh.write("\t".join(new_header) + "\n")
        with open(txt_path) as in_fh:
            for line in in_fh:
                if "Otherinfo" in line:    # skip old header line(s)
                    continue
                if line.startswith("##"):  # skip VCF meta-lines if present
                    continue
                out_fh.write(line)
                rows_written += 1
    return rows_written


def validate_column_counts(txt_path: str, out_path: str, expected: int):
    """Quick sanity check: first data row should have expected column count."""
    for path in (txt_path, out_path):
        n = count_columns(path)
        if n != expected:
            print(f"WARNING: {path} has {n} columns, expected {expected}",
                  file=sys.stderr)


def main():
    args = parse_args()

    # Extract headers
    vcf_header = extract_header(args.vcf, "CHROM")
    txt_header = extract_header(args.txt, "Otherinfo")

    print(f"[fix_vcf_header] VCF header columns : {len(vcf_header)}")
    print(f"[fix_vcf_header] TXT header columns : {len(txt_header)}")

    # Build merged header
    new_header = build_new_header(txt_header, vcf_header)
    print(f"[fix_vcf_header] New header columns : {len(new_header)}")

    # Write output
    n_rows = write_fixed_file(args.txt, new_header, args.out)
    print(f"[fix_vcf_header] Rows written       : {n_rows}")
    print(f"[fix_vcf_header] Output             : {args.out}")

    # Sanity check
    validate_column_counts(args.txt, args.out, len(new_header))


if __name__ == "__main__":
    main()
