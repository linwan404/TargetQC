#!/usr/bin/env bash
set -euo pipefail

# TargetQC.sh - Main workflow control script with portable conda setup and HTML report output.
# Put this file in the root of the TargetQC GitHub repository.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
SRC_DIR="${SCRIPT_DIR}/src"
REPORT_RMD_DEFAULT="${SCRIPT_DIR}/reports/TargetQC_report.Rmd"
PLOT_R_DEFAULT="${SRC_DIR}/targetqc_plot_flat.R"
RENDER_R_DEFAULT="${SRC_DIR}/render_TargetQC_report.R"
VERSION="1.1"

# -------------------------
# Default parameters
# -------------------------
PYTHON3_PATH=""
PYTHON2_PATH=""
BASH_PATH="bash"
HAP_PATH=""
MOSDEPTH=""
BEDTOOLS=""
RSCRIPT_PATH=""
THRESHOLDS="10,20,30,40"
UP_EXON=30
DOWN_EXON=20
BASE_EXON=5
UP_GENE=30
DOWN_GENE=20
WELL=95.0
POOR=5.0
UPSTREAM=150
DOWNSTREAM=150
GENE=""
REFERENCE=""
REF_BED=""
REF_VCF=""
REF_FASTA=""
PATHOGENIC=""
SAMPLE_TYPE="clinical"
MODE=""
SEQ_TYPE=""
USER_REGIONS=""
TARGET_REGIONS_BED=""
BAM_FILE=""
VCF_FILE=""
GVCF_FILE=""
OUT=""
CHROM=""
PLATFORM=""
REPORT=1
REPORT_RMD="$REPORT_RMD_DEFAULT"
PLOT_R="$PLOT_R_DEFAULT"
RENDER_R="$RENDER_R_DEFAULT"
REPORT_HTML=""
AUTO_CONDA=1
CONDA_ENV="targetqc"
HAPPY_CONDA_ENV="targetqc_happy"
ENV_YML="${SCRIPT_DIR}/environment.yml"
HAPPY_ENV_YML="${SCRIPT_DIR}/environment_happy.yml"

show_version() { echo "TargetQC version $VERSION"; exit 0; }

usage() {
    cat <<'EOF'
Quick Start:
Use TargetQC with user-provided files and generate an HTML report in one step.

Usage:
  TargetQC.sh --sample-type standard|clinical --seq ES|WGS --bam <file.bam> --vcf <file.vcf> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --pathogenic <pathogenic.txt> --out <outdir_or_outprefix> [options]

Required parameters:
  --sample-type <standard|clinical>       standard for benchmark sample; clinical for sample without benchmark
  --seq <ES|WGS>                         sequencing type
  --bam <file.bam>                       input BAM file
  --vcf <file.vcf>                       input VCF file
  --capture <capture.bed>                capture regions; required for ES
  --gene <gene.txt>                      target genes/transcripts file
  --reference <reference.txt>            gene annotation file
  --pathogenic <pathogenic.txt>          ClinVar/HGMD pathogenic site file
  --out <outdir_or_outprefix>            output directory/prefix

Required for standard sample:
  --gvcf <file.gvcf>                     input GVCF file
  --ref-vcf <ref_vcf>                    benchmark VCF file
  --ref-fasta <ref_fasta>                reference genome FASTA
  --ref-bed <ref_bed>                    benchmark BED file; required for WGS standard benchmark

Additional general options:
  --thresholds <string>                  coverage thresholds, default: 10,20,30,40
  --up-exon <int>                        well-covered exon threshold, default: 30
  --down-exon <int>                      poor-covered exon threshold, default: 20
  --base-exon <int>                      poor-covered base count in exon, default: 5
  --upstream <int>                       upstream extension length, default: 150
  --downstream <int>                     downstream extension length, default: 150
  --up-gene <int>                        well-covered gene threshold, default: 30
  --down-gene <int>                      poor-covered gene threshold, default: 20
  --well <float>                         percentage threshold for well-covered genes, default: 95.0
  --poor <float>                         percentage threshold for poor-covered genes, default: 5.0
  --chrom <string>                       restrict analysis to one chromosome
  --regions <regions.bed>                additional regions for variant QC/benchmark
  --platform <string>                    platform name shown in HTML report

Report options:
  --report                               generate HTML report, default: enabled
  --no-report                            do not generate HTML report
  --report-html <file.html>              custom HTML report path
  --report-rmd <file.Rmd>                custom Rmd template
  --plot-r <file.R>                      custom plot-generation R script
  --render-r <file.R>                    custom report-rendering R script

Conda options:
  --auto-conda                           create/activate conda env automatically, default: enabled
  --skip-conda                           do not touch conda; use current PATH or explicit executable paths
  --conda-env <name>                     main conda env, default: targetqc
  --happy-conda-env <name>               hap.py/Python2 env for standard mode, default: targetqc_happy
  --env-yml <environment.yml>            main env yml, default: ./environment.yml
  --happy-env-yml <environment_happy.yml> hap.py env yml, default: ./environment_happy.yml

Executable override options:
  --python3-path <path>                  Python3 executable
  --python2-path <path>                  Python2 executable for hap.py
  --hap-path <path>                      hap.py executable
  --mosdepth-path <path>                 mosdepth executable
  --bedtools-path <path>                 bedtools executable
  --rscript-path <path>                  Rscript executable

Modes:
  --mode exon
  --mode gene
  --mode variant-calling
  --mode variant-capture
  If --mode is not specified, TargetQC runs the full workflow.

  -h, --help                             show this help message
  --version                              show version
EOF
    exit 1
}

usage_exon() { echo "Usage: TargetQC.sh --mode exon --seq ES|WGS --bam <file.bam> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --out <out> [options]"; exit 1; }
usage_gene() { echo "Usage: TargetQC.sh --mode gene --seq ES|WGS --bam <file.bam> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --out <out> [options]"; exit 1; }
usage_variant_calling() { echo "Usage: TargetQC.sh --mode variant-calling --seq ES|WGS --vcf <file.vcf> --gvcf <file.gvcf> [--capture <capture.bed>|--ref-bed <ref.bed>] --ref-vcf <ref.vcf> --ref-fasta <ref.fa> --pathogenic <pathogenic.tsv> --out <out> [options]"; exit 1; }
usage_variant_capture() { echo "Usage: TargetQC.sh --mode variant-capture --seq ES|WGS --vcf <file.vcf> [--capture <capture.bed>] --pathogenic <pathogenic.tsv> --out <out> [options]"; exit 1; }

fail() { echo "Error: $*" >&2; exit 1; }
warn() { echo "[Warning] $*" >&2; }
need_file() { [[ -f "$1" ]] || fail "$2 not found: $1"; }
need_dir() { [[ -d "$1" ]] || fail "$2 not found: $1"; }
need_exe() { command -v "$1" >/dev/null 2>&1 || fail "Executable not found: $1"; }
script_path() { local f="$1"; [[ -f "${SRC_DIR}/${f}" ]] && printf '%s\n' "${SRC_DIR}/${f}" || fail "Required script not found: ${SRC_DIR}/${f}"; }
run_py3() { "$PYTHON3_PATH" "$(script_path "$1")" "${@:2}"; }
run_sh() { "$BASH_PATH" "$(script_path "$1")" "${@:2}"; }

conda_env_exists() { conda env list | awk '{print $1}' | grep -Fxq "$1"; }
conda_env_prefix() { conda env list | awk -v n="$1" '$1==n{print $NF; found=1} END{if(!found) exit 1}'; }
setup_conda() {
    [[ "$AUTO_CONDA" -eq 1 ]] || return 0
    need_exe conda
    local conda_base
    conda_base="$(conda info --base)"
    # shellcheck disable=SC1090
    source "${conda_base}/etc/profile.d/conda.sh"
    if ! conda_env_exists "$CONDA_ENV"; then
        need_file "$ENV_YML" "Main conda environment file"
        echo "[Conda] Creating main environment: $CONDA_ENV"
        conda env create -f "$ENV_YML"
    fi
    echo "[Conda] Activating main environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
    PYTHON3_PATH="${PYTHON3_PATH:-$(command -v python3 || command -v python)}"
    RSCRIPT_PATH="${RSCRIPT_PATH:-$(command -v Rscript || true)}"
    MOSDEPTH="${MOSDEPTH:-$(command -v mosdepth || true)}"
    BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || true)}"
    if [[ "$SAMPLE_TYPE" == "standard" || "$MODE" == "variant-calling" ]]; then
        if ! conda_env_exists "$HAPPY_CONDA_ENV"; then
            need_file "$HAPPY_ENV_YML" "hap.py conda environment file"
            echo "[Conda] Creating hap.py/Python2 environment: $HAPPY_CONDA_ENV"
            conda env create -f "$HAPPY_ENV_YML"
        fi
        local happy_prefix
        happy_prefix="$(conda_env_prefix "$HAPPY_CONDA_ENV")"
        PYTHON2_PATH="${PYTHON2_PATH:-${happy_prefix}/bin/python}"
        HAP_PATH="${HAP_PATH:-${happy_prefix}/bin/hap.py}"
    fi
}

resolve_executables() {
    PYTHON3_PATH="${PYTHON3_PATH:-$(command -v python3 || command -v python || true)}"
    BASH_PATH="${BASH_PATH:-$(command -v bash || true)}"
    MOSDEPTH="${MOSDEPTH:-$(command -v mosdepth || true)}"
    BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || true)}"
    RSCRIPT_PATH="${RSCRIPT_PATH:-$(command -v Rscript || true)}"
    if [[ "$SAMPLE_TYPE" == "standard" || "$MODE" == "variant-calling" ]]; then
        PYTHON2_PATH="${PYTHON2_PATH:-$(command -v python2 || true)}"
        HAP_PATH="${HAP_PATH:-$(command -v hap.py || true)}"
    fi
    [[ -n "$PYTHON3_PATH" ]] || fail "Python3 is required. Use --auto-conda or --python3-path."
    [[ -n "$BASH_PATH" ]] || fail "bash is required."
    if [[ "$MODE" == "exon" || "$MODE" == "gene" || "$MODE" == "both" ]]; then [[ -n "$MOSDEPTH" ]] || fail "mosdepth is required. Use --auto-conda or --mosdepth-path."; fi
    if [[ "$SAMPLE_TYPE" == "standard" || "$MODE" == "variant-calling" ]]; then [[ -n "$PYTHON2_PATH" ]] || fail "Python2 is required for hap.py. Use --auto-conda or --python2-path."; [[ -n "$HAP_PATH" ]] || fail "hap.py is required. Use --auto-conda or --hap-path."; fi
    if [[ "$REPORT" -eq 1 ]]; then [[ -n "$RSCRIPT_PATH" ]] || fail "Rscript is required for HTML report. Use --auto-conda or --rscript-path."; need_file "$PLOT_R" "Plot R script"; need_file "$RENDER_R" "Report render R script"; need_file "$REPORT_RMD" "Report Rmd template"; fi
}

[[ $# -gt 0 ]] || usage
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) BAM_FILE="$2"; shift 2 ;;
        --capture) TARGET_REGIONS_BED="$2"; shift 2 ;;
        --regions) USER_REGIONS="$2"; shift 2 ;;
        --seq) SEQ_TYPE="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --thresholds) THRESHOLDS="$2"; shift 2 ;;
        --chrom) CHROM="$2"; shift 2 ;;
        --up-exon) UP_EXON="$2"; shift 2 ;;
        --down-exon) DOWN_EXON="$2"; shift 2 ;;
        --base-exon) BASE_EXON="$2"; shift 2 ;;
        --upstream) UPSTREAM="$2"; shift 2 ;;
        --downstream) DOWNSTREAM="$2"; shift 2 ;;
        --up-gene) UP_GENE="$2"; shift 2 ;;
        --down-gene) DOWN_GENE="$2"; shift 2 ;;
        --well) WELL="$2"; shift 2 ;;
        --poor) POOR="$2"; shift 2 ;;
        --gene) GENE="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --sample-type) SAMPLE_TYPE="$2"; shift 2 ;;
        --vcf) VCF_FILE="$2"; shift 2 ;;
        --gvcf) GVCF_FILE="$2"; shift 2 ;;
        --ref-bed) REF_BED="$2"; shift 2 ;;
        --ref-vcf) REF_VCF="$2"; shift 2 ;;
        --ref-fasta) REF_FASTA="$2"; shift 2 ;;
        --pathogenic) PATHOGENIC="$2"; shift 2 ;;
        --platform) PLATFORM="$2"; shift 2 ;;
        --report) REPORT=1; shift ;;
        --no-report) REPORT=0; shift ;;
        --report-html) REPORT_HTML="$2"; shift 2 ;;
        --report-rmd) REPORT_RMD="$2"; shift 2 ;;
        --plot-r) PLOT_R="$2"; shift 2 ;;
        --render-r) RENDER_R="$2"; shift 2 ;;
        --auto-conda) AUTO_CONDA=1; shift ;;
        --skip-conda) AUTO_CONDA=0; shift ;;
        --conda-env) CONDA_ENV="$2"; shift 2 ;;
        --happy-conda-env) HAPPY_CONDA_ENV="$2"; shift 2 ;;
        --env-yml) ENV_YML="$2"; shift 2 ;;
        --happy-env-yml) HAPPY_ENV_YML="$2"; shift 2 ;;
        --python3-path) PYTHON3_PATH="$2"; shift 2 ;;
        --python2-path) PYTHON2_PATH="$2"; shift 2 ;;
        --hap-path) HAP_PATH="$2"; shift 2 ;;
        --mosdepth-path) MOSDEPTH="$2"; shift 2 ;;
        --bedtools-path) BEDTOOLS="$2"; shift 2 ;;
        --rscript-path) RSCRIPT_PATH="$2"; shift 2 ;;
        --version) show_version ;;
        -h|--help) case "$MODE" in exon) usage_exon ;; gene) usage_gene ;; variant-calling) usage_variant_calling ;; variant-capture) usage_variant_capture ;; *) usage ;; esac ;;
        *) fail "Unknown option: $1" ;;
    esac
done

[[ -n "$MODE" ]] || { MODE="both"; echo "--mode parameter not specified, defaulting to full workflow."; }
[[ "$SAMPLE_TYPE" == "standard" || "$SAMPLE_TYPE" == "clinical" ]] || fail "--sample-type must be standard or clinical."
[[ "$SEQ_TYPE" == "ES" || "$SEQ_TYPE" == "WGS" ]] || fail "--seq must be ES or WGS."
[[ -n "$OUT" ]] || fail "--out is required."
[[ -d "$SRC_DIR" ]] || fail "src directory not found: $SRC_DIR"

if [[ "$MODE" == "exon" || "$MODE" == "gene" || "$MODE" == "both" ]]; then
    [[ -n "$BAM_FILE" ]] || fail "--bam is required for exon/gene/full workflow."
    [[ -n "$GENE" ]] || fail "--gene is required for exon/gene/full workflow."
    [[ -n "$REFERENCE" ]] || fail "--reference is required for exon/gene/full workflow."
fi
if [[ "$MODE" == "variant-calling" || "$MODE" == "variant-capture" || "$MODE" == "both" ]]; then
    [[ -n "$VCF_FILE" ]] || fail "--vcf is required for variant/full workflow."
    [[ -n "$PATHOGENIC" ]] || fail "--pathogenic is required for variant/full workflow."
fi
if [[ "$MODE" == "variant-calling" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "standard" ) ]]; then
    [[ -n "$GVCF_FILE" ]] || fail "--gvcf is required for standard/variant-calling workflow."
    [[ -n "$REF_VCF" ]] || fail "--ref-vcf is required for standard/variant-calling workflow."
    [[ -n "$REF_FASTA" ]] || fail "--ref-fasta is required for standard/variant-calling workflow."
    if [[ "$SEQ_TYPE" == "WGS" ]]; then [[ -n "$REF_BED" ]] || fail "--ref-bed is required for WGS standard benchmark."; fi
fi
if [[ "$SEQ_TYPE" == "ES" ]]; then [[ -n "$TARGET_REGIONS_BED" ]] || fail "ES requires --capture."; fi

[[ -n "$BAM_FILE" ]] && need_file "$BAM_FILE" "BAM file"
[[ -n "$VCF_FILE" ]] && need_file "$VCF_FILE" "VCF file"
[[ -n "$GVCF_FILE" ]] && need_file "$GVCF_FILE" "GVCF file"
[[ -n "$TARGET_REGIONS_BED" ]] && need_file "$TARGET_REGIONS_BED" "Capture BED"
[[ -n "$USER_REGIONS" ]] && need_file "$USER_REGIONS" "User regions BED"
[[ -n "$GENE" ]] && need_file "$GENE" "Target gene file"
[[ -n "$REFERENCE" ]] && need_file "$REFERENCE" "Gene annotation file"
[[ -n "$PATHOGENIC" ]] && need_file "$PATHOGENIC" "Pathogenic variant file"
[[ -n "$REF_BED" ]] && need_file "$REF_BED" "Reference BED"
[[ -n "$REF_VCF" ]] && need_file "$REF_VCF" "Reference VCF"
[[ -n "$REF_FASTA" ]] && need_file "$REF_FASTA" "Reference FASTA"

setup_conda
resolve_executables

TEMP_DIR="${OUT}/temp"
OUTPUT_DIR="$OUT"
OUTPREFIX="$(basename "$OUT")"
mkdir -p "$TEMP_DIR" "$OUTPUT_DIR"

if [[ -z "$REPORT_HTML" ]]; then REPORT_HTML="${OUT}/${OUTPREFIX}_TargetQC_report.html"; fi

echo "Starting TargetQC workflow..."
echo "TargetQC version: $VERSION"
echo "Repository root: $SCRIPT_DIR"
echo "Source directory: $SRC_DIR"
echo "Analysis mode: $MODE"
echo "Sample type: $SAMPLE_TYPE"
echo "Sequencing type: $SEQ_TYPE"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $OUTPREFIX"
echo "Python3: $PYTHON3_PATH"
[[ "$SAMPLE_TYPE" == "standard" || "$MODE" == "variant-calling" ]] && echo "hap.py: $HAP_PATH"
[[ "$REPORT" -eq 1 ]] && echo "HTML report: $REPORT_HTML"

# -------------------------
# Target Exon QC
# -------------------------
if [[ "$MODE" == "exon" || "$MODE" == "both" ]]; then
    echo "Executing target exon QC analysis..."
    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next} $2=="exon" && $6 in transcript && $8==transcript[$6]{print $1 "\t" $3-1 "\t" $4}' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_exon.bed"
    run_py3 target_coverage.py --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_EXON,$UP_EXON" --out "${TEMP_DIR}/${OUTPREFIX}_exon"
    "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_exon.sh"
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        run_py3 target_coverage.py --bam "$BAM_FILE" --regions "$TARGET_REGIONS_BED" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture"
        "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_capture.sh"
    fi
    if [[ -n "$CHROM" ]]; then
        if [[ "$SEQ_TYPE" == "ES" ]]; then
            run_py3 base.py --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture.regions.bed.gz" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
        else
            run_py3 target_coverage.py --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture_exon"
            "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_capture_exon.sh"
            run_py3 base.py --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.regions.bed.gz" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
        fi
        run_py3 exon_coverage.py --thresholds-exon "${TEMP_DIR}/${OUTPREFIX}_exon.thresholds.bed.gz" --regions-exon "${TEMP_DIR}/${OUTPREFIX}_exon.regions.bed.gz" --up-exon "$UP_EXON" --down-exon "$DOWN_EXON" --base-exon "$BASE_EXON" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
    else
        if [[ "$SEQ_TYPE" == "ES" ]]; then
            run_py3 base.py --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture.regions.bed.gz" --out "${TEMP_DIR}/${OUTPREFIX}"
        else
            run_py3 target_coverage.py --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture_exon"
            "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_capture_exon.sh"
            run_py3 base.py --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.regions.bed.gz" --out "${TEMP_DIR}/${OUTPREFIX}"
        fi
        run_py3 exon_coverage.py --thresholds-exon "${TEMP_DIR}/${OUTPREFIX}_exon.thresholds.bed.gz" --regions-exon "${TEMP_DIR}/${OUTPREFIX}_exon.regions.bed.gz" --up-exon "$UP_EXON" --down-exon "$DOWN_EXON" --base-exon "$BASE_EXON" --out "${TEMP_DIR}/${OUTPREFIX}"
    fi
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        run_sh exon_capture_QC.sh --coverage-exon "${TEMP_DIR}/${OUTPREFIX}_exon_coverage.tsv" --regions "$TARGET_REGIONS_BED" --out "${TEMP_DIR}/${OUTPREFIX}"
    else
        awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "#chrom", "start", "end", "coverage", "type", "capture"; next} {print $1, $2, $3, $4, $5, "full_capture"}' "${TEMP_DIR}/${OUTPREFIX}_exon_coverage.tsv" > "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv"
    fi
    run_py3 exon_summary.py --qc-exon "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv" --out "${TEMP_DIR}/${OUTPREFIX}"
fi

# -------------------------
# Target Gene QC
# -------------------------
if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    echo "Executing target gene QC analysis..."
    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next} $2=="exon" && $6 in transcript && $8==transcript[$6]{print $1 "\t" $3-1 "\t" $4}' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_exon.bed"
    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next} $2=="CDS" && $6 in transcript && $8==transcript[$6]{print $1 "\t" $3-1 "\t" $4}' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_cds.bed"
    run_py3 target_coverage.py --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_GENE,$UP_GENE" --out "${TEMP_DIR}/${OUTPREFIX}_target_exon"
    "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_target_exon.sh"
    run_py3 target_coverage.py --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_cds.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_GENE,$UP_GENE" --out "${TEMP_DIR}/${OUTPREFIX}_target_cds"
    "$BASH_PATH" "${TEMP_DIR}/${OUTPREFIX}_target_cds.sh"
    run_py3 trans.py --coverage "${TEMP_DIR}/${OUTPREFIX}_target_exon.thresholds.bed.gz" --gene "$GENE" --reference "$REFERENCE" --up-gene "$UP_GENE" --down-gene "$DOWN_GENE" --out "${TEMP_DIR}/${OUTPREFIX}"
    run_py3 cds.py --coverage "${TEMP_DIR}/${OUTPREFIX}_target_cds.thresholds.bed.gz" --gene "$GENE" --reference "$REFERENCE" --up "$UP_GENE" --down "$DOWN_GENE" --out "${TEMP_DIR}/${OUTPREFIX}"
    awk -F'\t' '$3=="protein_coding"' "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" > "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon.tsv" || true
    run_py3 trans_proportion.py --gene "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans"
    run_py3 trans_proportion.py --gene "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon"
    run_py3 trans_proportion.py --gene "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans_cds"
    run_py3 gene_capture_QC.py --trans-exon "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" --trans-cds "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" --well "$WELL" --poor "$POOR" --out "${TEMP_DIR}/${OUTPREFIX}"
    run_py3 gene_capture_QC_summary.py --gene-qc "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.tsv" --out "${TEMP_DIR}/${OUTPREFIX}"
fi

# -------------------------
# Variant Calling QC for standard sample
# -------------------------
if [[ "$MODE" == "variant-calling" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "standard" ) ]]; then
    echo "Executing variant calling QC analysis..."
    if [[ "$SEQ_TYPE" == "ES" ]]; then REGIONS_FILE="$TARGET_REGIONS_BED"; else REGIONS_FILE="$REF_BED"; fi
    run_py3 benchmark.py --vcf "$VCF_FILE" --regions "$REGIONS_FILE" --ref-vcf "$REF_VCF" --ref-fasta "$REF_FASTA" --tools "$HAP_PATH" --python2-path "$PYTHON2_PATH" --out "${TEMP_DIR}/${OUTPREFIX}_capture"
    run_py3 variant_calling_summary.py --input "${TEMP_DIR}/${OUTPREFIX}_capture.summary.csv" --out "${TEMP_DIR}/${OUTPREFIX}_variant_calling"
    if [[ -n "$USER_REGIONS" ]]; then
        run_py3 benchmark.py --vcf "$VCF_FILE" --regions "$USER_REGIONS" --ref-vcf "$REF_VCF" --ref-fasta "$REF_FASTA" --tools "$HAP_PATH" --python2-path "$PYTHON2_PATH" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        run_py3 variant_calling_summary.py --input "${TEMP_DIR}/${OUTPREFIX}_regions.summary.csv" --out "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling"
    fi
    run_py3 variant_calling_QC.py --vcf-orig "$GVCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "FN" --out "${TEMP_DIR}/${OUTPREFIX}"
    run_py3 variant_calling_QC.py --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "FP" --out "${TEMP_DIR}/${OUTPREFIX}"
    run_py3 variant_calling_QC.py --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "TP" --out "${TEMP_DIR}/${OUTPREFIX}"
    { echo -e "#chrom\tpos\ttype\tvariant_type\tgenotype\tDP\tBAF"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"TP",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_TP_QC.tsv"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"FP",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_FP_QC.tsv"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"FN",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_FN_QC.tsv"; } > "${TEMP_DIR}/${OUTPREFIX}_variant_calling_QC.tsv"
    if [[ -n "$USER_REGIONS" ]]; then
        run_py3 variant_calling_QC.py --vcf-orig "$GVCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "FN" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        run_py3 variant_calling_QC.py --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "FP" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        run_py3 variant_calling_QC.py --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "TP" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        { echo -e "#chrom\tpos\ttype\tvariant_type\tgenotype\tDP\tBAF"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"TP",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_regions_TP_QC.tsv"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"FP",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_regions_FP_QC.tsv"; awk -F'\t' -v OFS='\t' 'NR>1{print $1,$2,"FN",$3,$4,$5,$6}' "${TEMP_DIR}/${OUTPREFIX}_regions_FN_QC.tsv"; } > "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling_QC.tsv"
    fi
    run_py3 pathogenic.py --fn "${TEMP_DIR}/${OUTPREFIX}_FN_QC.tsv" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"
    if [[ -n "$USER_REGIONS" ]]; then run_py3 pathogenic.py --fn "${TEMP_DIR}/${OUTPREFIX}_regions_FN_QC.tsv" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}_region"; fi
fi

# -------------------------
# Variant Capture QC for clinical sample
# -------------------------
if [[ "$MODE" == "variant-capture" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "clinical" ) ]]; then
    echo "Executing clinical variant-site QC analysis..."
    if [[ "$SEQ_TYPE" == "ES" ]]; then run_py3 variant_capture_QC.py --vcf "$VCF_FILE" --regions "$TARGET_REGIONS_BED" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"; else run_py3 variant_capture_QC.py --vcf "$VCF_FILE" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"; fi
    if [[ -n "$USER_REGIONS" ]]; then run_py3 variant_capture_QC.py --vcf "$VCF_FILE" --regions "$USER_REGIONS" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}_region"; fi
fi

# -------------------------
# Copy main output files
# -------------------------
echo "Copying main output files to ${OUTPUT_DIR}..."
copy_if_exists() { [[ -f "$1" ]] && cp "$1" "$2"; }
if [[ "$MODE" == "exon" || "$MODE" == "both" ]]; then
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_base.tsv" "${OUT}/${OUTPREFIX}_base.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv" "${OUT}/${OUTPREFIX}_exon_capture_QC.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.summary.tsv" "${OUT}/${OUTPREFIX}_exon_capture_QC.summary.tsv"
fi
if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" "${OUT}/${OUTPREFIX}_trans.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" "${OUT}/${OUTPREFIX}_trans_cds.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.tsv" "${OUT}/${OUTPREFIX}_gene_capture_QC.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.summary.tsv" "${OUT}/${OUTPREFIX}_gene_capture_QC.summary.tsv"
fi
if [[ "$MODE" == "variant-calling" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "standard" ) ]]; then
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_variant_calling.summary.tsv" "${OUT}/${OUTPREFIX}_variant_calling.summary.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling.summary.tsv" "${OUT}/${OUTPREFIX}_region_variant_calling.summary.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_variant_calling_QC.tsv" "${OUT}/${OUTPREFIX}_variant_calling_QC.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling_QC.tsv" "${OUT}/${OUTPREFIX}_region_variant_calling_QC.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_pathogenic.tsv" "${OUT}/${OUTPREFIX}_pathogenic.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_region_pathogenic.tsv" "${OUT}/${OUTPREFIX}_region_pathogenic.tsv"
fi
if [[ "$MODE" == "variant-capture" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "clinical" ) ]]; then
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_variant_capture_QC.tsv" "${OUT}/${OUTPREFIX}_variant_capture_QC.tsv"
    copy_if_exists "${TEMP_DIR}/${OUTPREFIX}_region_variant_capture_QC.tsv" "${OUT}/${OUTPREFIX}_region_variant_capture_QC.tsv"
fi

# -------------------------
# HTML report
# -------------------------
if [[ "$REPORT" -eq 1 ]]; then
    echo "Generating TargetQC plots..."
    "$RSCRIPT_PATH" "$PLOT_R" --input-root "$OUT" --output-root "${OUT}/figures" --prefix "$OUTPREFIX" --sample-type "$SAMPLE_TYPE" --seq "$SEQ_TYPE" --dp-high "$UP_EXON" --dp-low "$DOWN_EXON" --baf-threshold 0.05 --gene-level-range "≥${UP_GENE}X" --selected-gene-type "protein_coding"
    echo "Rendering TargetQC HTML report..."
    "$RSCRIPT_PATH" "$RENDER_R" --rmd "$REPORT_RMD" --outdir "$OUT" --plot-dir "${OUT}/figures" --prefix "$OUTPREFIX" --sample-name "$OUTPREFIX" --sample-type "$SAMPLE_TYPE" --seq "$SEQ_TYPE" --platform "$PLATFORM" --dp-high "$UP_EXON" --dp-low "$DOWN_EXON" --baf-threshold 0.05 --gene-level-range "≥${UP_GENE}X" --selected-gene-type "protein_coding" --out-html "$REPORT_HTML"
fi

echo "TargetQC workflow completed!"
echo "Main output files are in: $OUT"
[[ "$REPORT" -eq 1 ]] && echo "HTML report: $REPORT_HTML"
