#!/usr/bin/env bash
# Example launcher for the DRAGEN server.
#
#   bash examples/run_dragen.sh /path/to/vcf_dir_or_samplesheet /path/to/outdir [extra nextflow args]
#
# Uses the 'nf' conda environment for Nextflow + Java 21, the 'dragen' site
# profile (ANNOVAR paths and installed database versions) and the pinned
# conda environment for bcftools/htslib/perl/python.
set -euo pipefail

INPUT=${1:?usage: run_dragen.sh <input> <outdir> [nextflow args]}
OUTDIR=${2:?usage: run_dragen.sh <input> <outdir> [nextflow args]}
shift 2

PIPELINE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NF_ENV=/staging/anaconda3/envs/nf

export JAVA_HOME=${NF_ENV}
export NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR:-/staging/adaniyarov/annovar-nf/conda_cache}
# No "new version available" notice. NXF_ANSI_LOG is left unset so Nextflow shows
# its live progress display in a terminal and plain lines when redirected to a file.
export NXF_DISABLE_CHECK_LATEST=true

"${NF_ENV}/bin/nextflow" run "${PIPELINE_DIR}/main.nf" \
    -profile dragen,conda \
    --input  "${INPUT}" \
    --outdir "${OUTDIR}" \
    -resume \
    "$@"
