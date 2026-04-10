#!/data/miniconda3/envs/wgs_env/bin/python3
"""
summarize_annotations.py
------------------------
Cross-sample summary of ANNOVAR annotated files.
Produces a wide-format matrix: rows = Func.refGene x ExonicFunc.refGene,
columns = samples (one per annotated file).

Output:
    <out>.tsv              wide matrix
    <out>.html             interactive HTML table
    <out>_totals.tsv       simple per-sample totals
"""

import argparse
import os
import re
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--out", required=True, help="Output prefix (no extension)")
    return p.parse_args()


def get_sample_name(filepath: str) -> str:
    name = os.path.basename(filepath)
    name = re.sub(r'\.annotated\.header\.txt$', '', name)
    name = re.sub(r'\.hard-filtered$', '', name)
    name = re.sub(r'\.FINAL$', '', name)
    return name


def process_file(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep='\t', low_memory=False)
    func_col   = 'Func.refGene'       if 'Func.refGene'       in df.columns else df.columns[5]
    exfunc_col = 'ExonicFunc.refGene' if 'ExonicFunc.refGene' in df.columns else None
    if exfunc_col:
        grp = df.groupby([func_col, exfunc_col]).size().reset_index(name='counts')
        grp.columns = ['Func.refGene', 'ExonicFunc.refGene', 'counts']
    else:
        grp = df.groupby([func_col]).size().reset_index(name='counts')
        grp['ExonicFunc.refGene'] = '.'
        grp = grp[['Func.refGene', 'ExonicFunc.refGene', 'counts']]
    return grp


def build_wide_matrix(file_list: list) -> pd.DataFrame:
    unique_groups = pd.DataFrame()
    for fp in file_list:
        grp = process_file(fp)
        unique_groups = pd.concat([unique_groups,
                                   grp[['Func.refGene', 'ExonicFunc.refGene']]])
    unique_groups = unique_groups.drop_duplicates().reset_index(drop=True)

    results = unique_groups.copy()
    for fp in file_list:
        sample = get_sample_name(fp)
        print(f"[summarize] {sample}", file=sys.stderr)
        grp = process_file(fp)
        results = results.merge(grp, on=['Func.refGene', 'ExonicFunc.refGene'], how='left')
        results = results.rename(columns={'counts': sample})
        results[sample] = results[sample].fillna(0).astype(int)

    priority = {'exonic': 0, 'splicing': 1, 'exonic;splicing': 2}
    results['_s'] = results['Func.refGene'].map(lambda x: priority.get(x, 99))
    results = results.sort_values(['_s', 'Func.refGene', 'ExonicFunc.refGene']).drop(columns='_s').reset_index(drop=True)
    return results


def build_totals(file_list: list) -> pd.DataFrame:
    rows = []
    for fp in file_list:
        sample = get_sample_name(fp)
        df = pd.read_csv(fp, sep='\t', low_memory=False)
        row = {'sample': sample, 'total_variants': len(df)}
        func_col = 'Func.refGene' if 'Func.refGene' in df.columns else None
        if func_col:
            for cat in ['exonic', 'splicing', 'intronic', 'intergenic', 'UTR3', 'UTR5']:
                row[cat] = int((df[func_col].str.lower() == cat).sum())
        for col in df.columns:
            if 'clnsig' in col.lower() or 'clinvar' in col.lower():
                row['clinvar_pathogenic'] = int(
                    df[col].astype(str).str.contains('Pathogenic', na=False).sum()
                )
                break
        if 'REVEL' in df.columns:
            row['revel_high_075'] = int(
                pd.to_numeric(df['REVEL'], errors='coerce').ge(0.75).sum()
            )
        rows.append(row)
    return pd.DataFrame(rows).fillna(0)


def write_html(df: pd.DataFrame, out_path: str):
    n_samples = len(df.columns) - 2
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>annovar-nf summary</title>
<style>
  body{{font-family:Arial,sans-serif;margin:30px;font-size:13px}}
  h1{{color:#2c3e50}}
  table{{border-collapse:collapse;width:100%}}
  th,td{{border:1px solid #ddd;padding:5px 10px;text-align:right;white-space:nowrap}}
  th{{background:#2c3e50;color:white;position:sticky;top:0}}
  td:first-child,td:nth-child(2){{text-align:left;font-weight:bold}}
  tr:nth-child(even){{background:#f5f5f5}}
  tr:hover{{background:#d6eaf8}}
</style></head><body>
<h1>annovar-nf — Functional variant summary</h1>
<p>Samples: <b>{n_samples}</b></p>
{df.to_html(index=False, border=0, na_rep='0')}
</body></html>"""
    with open(out_path, 'w') as fh:
        fh.write(html)


def main():
    args = parse_args()
    wide = build_wide_matrix(args.inputs)
    wide.to_csv(args.out + ".tsv", index=False, sep='\t')
    write_html(wide, args.out + ".html")
    print(f"[summarize] Wide matrix  → {args.out}.tsv  ({len(wide)} rows)")

    totals = build_totals(args.inputs)
    totals.to_csv(args.out + "_totals.tsv", index=False, sep='\t')
    print(f"[summarize] Totals table → {args.out}_totals.tsv")
    print(totals.to_string(index=False))


if __name__ == "__main__":
    main()
