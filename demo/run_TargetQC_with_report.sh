#!/usr/bin/env bash
set -euo pipefail

# Example caller for a fresh user.
# This script clones TargetQC, lets TargetQC.sh create/activate conda environments,
# runs TargetQC, creates plots, and renders an HTML report.

REPO_URL="https://github.com/linwan404/TargetQC.git"
REPO_DIR="TargetQC"

# -------------------------
# User inputs: edit these paths
# -------------------------
SAMPLE_NAME="sample01"
SAMPLE_TYPE="clinical"     # clinical or standard
SEQ_TYPE="ES"              # ES or WGS
BAM="/path/to/sample.bam"
VCF="/path/to/sample.vcf.gz"
GVCF="/path/to/sample.g.vcf.gz"                 # required only for standard
CAPTURE="/path/to/capture.bed"                  # required for ES
GENE="/path/to/gene.txt"
REFERENCE="/path/to/reference.gtf"
PATHOGENIC="/path/to/pathogenic.tsv"
REF_BED="/path/to/benchmark.bed"                # required for WGS standard; optional for ES standard
REF_VCF="/path/to/benchmark.vcf.gz"             # required for standard
REF_FASTA="/path/to/reference.fa"               # required for standard
REGIONS=""                                      # optional additional target regions
PLATFORM="MyPlatform"
OUT="${PWD}/TargetQC_result/${SAMPLE_NAME}"

# -------------------------
# Install/fetch TargetQC
# -------------------------
if [[ ! -d "$REPO_DIR" ]]; then
    git clone "$REPO_URL" "$REPO_DIR"
fi
cd "$REPO_DIR"
chmod +x TargetQC.sh

COMMON_ARGS=(--sample-type "$SAMPLE_TYPE" --seq "$SEQ_TYPE" --bam "$BAM" --vcf "$VCF" --gene "$GENE" --reference "$REFERENCE" --pathogenic "$PATHOGENIC" --out "$OUT" --platform "$PLATFORM" --report --auto-conda)
if [[ "$SEQ_TYPE" == "ES" ]]; then COMMON_ARGS+=(--capture "$CAPTURE"); fi
if [[ -n "$REGIONS" ]]; then COMMON_ARGS+=(--regions "$REGIONS"); fi

if [[ "$SAMPLE_TYPE" == "standard" ]]; then
    COMMON_ARGS+=(--gvcf "$GVCF" --ref-vcf "$REF_VCF" --ref-fasta "$REF_FASTA")
    if [[ -n "$REF_BED" ]]; then COMMON_ARGS+=(--ref-bed "$REF_BED"); fi
fi

./TargetQC.sh "${COMMON_ARGS[@]}"

echo "Done. HTML report: ${OUT}/$(basename "$OUT")_TargetQC_report.html"
