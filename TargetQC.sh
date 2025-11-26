#!/bin/bash

# TargetQC.sh - Main workflow control script
#!/bin/bash

# TargetQC.sh - Main workflow control script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CONDA_ENV_PATH="$CONDA_PREFIX"

# Version information
VERSION="1.0"

# Default parameters
PYTHON3_PATH="python3"
PYTHON2_PATH="python2"
BASH_PATH="bash"
HAP_PATH="$CONDA_ENV_PATH/bin/hap.py"
MOSDEPTH="$CONDA_ENV_PATH/bin/mosdepth"
BEDTOOLS="bedtools"
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
Hap="$CONDA_ENV_PATH/bin/hap.py"
Python2="python2"

show_version() {
    echo "TargetQC version $VERSION"
    exit 0
}

usage() {
    echo "Quick Start:"
    echo "Use default parameters to complete the analysis in one step"
    echo""
    echo "Usage: TargetQC.sh --sample-type standard|clinical --seq ES|WGS --bam <file.bam> --vcf <file.vcf> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --pathogenic <pathogenic.txt> --out <outprefix> [options]"
    echo ""
    echo "Required parameters:"
    echo ""
    echo "  --sample-type <standard|clinical>            Sample type, standard for sample with benchmark, clinical for sample without benchmark"
    echo "  --seq <ES|WGS>                 Sequencing type, ES for exome sequencing, WGS for whole genome sequencing"
    echo "  --bam <bam_file>               Input BAM files of sample"
    echo "  --vcf <vcf_file>               Input VCF files of sample"
    echo "  --capture <capture.bed>        Capture regions (.bed). Required for ES sequencing type"
    echo "  --gene <gene_file>             Target genes. Such as protein-coding genes with clear genetic patterns in OMIM database"
    echo "  --reference <reference_file>   Gene annotation file"
    echo "  --pathogenic <pathogenic_file> Target variant sites"
    echo "  --out <outprefix>              Output prefix"
    echo ""
    echo "Additional general options:"
    echo ""
    echo "  --thresholds <string>          Coverage calculation thresholds, must be integers (default: 10,20,30,40)"
    echo "  --up-exon <int>                Coverage threshold for well-covered exon regions (default: 30)"
    echo "  --down-exon <int>              Coverage threshold for poor-covered exon regions (default: 20)"
    echo "  --base-exon <int>              Number of poor-covered bases in exons (default: 5)"
    echo "  --upstream <int>               Upstream extension length for capture regions (default: 150)"
    echo "  --downstream <int>             Downstream extension length for capture regions (default: 150)"
    echo "  --up-gene <int>                Coverage threshold for well-covered gene regions (default: 30)"
    echo "  --down-gene <int>              Coverage threshold for poor-covered gene regions (default: 20)"
    echo "  --well <float>                 Percentage of well-covered regions (default: 95.0)"
    echo "  --poor <float>                 Percentage of poorly-covered regions (default: 5.0)"
    echo "  --chrom <string>               Only genotype regions on this chromosome"
    echo "  --regions <regions.bed>        Additional regions for variant detection (optional)"
    echo "  --ref-bed <ref_bed>            Benchmark BED file (required for WGS when --sample-type is standard)"
    echo "  --ref-vcf <ref_vcf>            Benchmark VCF file (required when --sample-type is standard)"
    echo "  --ref-fasta <ref_fasta>        Human reference genome fasta file (required when --sample-type is standard)"
    echo "  --gvcf <gvcf_file>             Input GVCF files of sample(required when --sample-type is standard)"
    echo "  --version                      Show version information"
    echo "  -h, --help                     Show this help message"
    echo ""
    echo "Modes:"
    echo "  1. TargetQC.sh --mode exon "
    echo "  2. TargetQC.sh --mode gene "
    echo "  3. TargetQC.sh --mode variant-calling "
    echo "  4. TargetQC.sh --mode variant-capture "
    echo ""
    echo "Use -h with specific mode for detailed help:"
    echo "  TargetQC.sh --mode exon -h"
    echo "  TargetQC.sh --mode gene -h"
    echo "  TargetQC.sh --mode variant-calling -h"
    echo "  TargetQC.sh --mode variant-capture -h"
    exit 1
}

usage_exon() {
    echo "1. Exon capture QC:"
    echo "The quality control of exon coverage"
    echo ""
    echo "Usage: TargetQC.sh --mode exon --seq ES|WGS --bam <file.bam> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --out <outprefix> [options]"
    echo ""
    echo "Required parameters:"
    echo ""
    echo "  --mode exon                    Run TargetQC in Exon capture QC mode. This mode should be used when performing quality statistics of exon coverage analysis"
    echo "  --seq <ES|WGS>                 Sequencing type, ES for exome sequencing, WGS for whole genome sequencing"
    echo "  --bam <bam_file>               Input BAM files of sample"
    echo "  --capture <capture.bed>        Capture regions (.bed). Required for ES sequencing type"
    echo "  --gene <gene_file>             Target genes. Such as protein-coding genes with clear genetic patterns in OMIM database"
    echo "  --reference <reference_file>   Gene annotation file"
    echo "  --out <outprefix>              Output prefix"
    echo ""
    echo "Additional general options:"
    echo ""
    echo "  --thresholds <string>          Coverage calculation thresholds, must be integers (default: 10,20,30,40)"
    echo "  --up-exon <int>                Coverage threshold for well-covered exon regions (default: 30)"
    echo "  --down-exon <int>              Coverage threshold for poor-covered exon regions (default: 20)"
    echo "  --base-exon <int>              Number of poor-covered bases in exons (default: 5)"
    echo "  --upstream <int>               Upstream extension length for capture regions (default: 150)"
    echo "  --downstream <int>             Downstream extension length for capture regions (default: 150)"
    echo "  --chrom <string>               Only genotype regions on this chromosome"
    echo "  -h, --help                     Show this help message"
    exit 1
}

usage_gene() {
    echo "2. Gene capture QC:"
    echo "The quality control of gene coverage"
    echo ""
    echo "Usage: TargetQC.sh --mode gene --seq ES|WGS --bam <file.bam> [--capture <capture.bed>] --gene <gene.txt> --reference <reference.txt> --out <outprefix> [options]"
    echo ""
    echo "Required parameters:"
    echo ""
    echo "  --mode gene                    Run TargetQC in Exon capture QC mode. This mode should be used when performing quality statistics of exon coverage analysis"
    echo "  --seq <ES|WGS>                 Sequencing type, ES for exome sequencing, WGS for whole genome sequencing"
    echo "  --bam <bam_file>               Input BAM files of sample"
    echo "  --capture <capture.bed>        Capture regions (.bed). Required for ES sequencing type"
    echo "  --gene <gene_file>             Target genes. Such as protein-coding genes with clear genetic patterns in OMIM database"
    echo "  --reference <reference_file>   Gene annotation file"
    echo "  --out <outprefix>              Output prefix"
    echo ""
    echo "Additional general options:"
    echo ""
    echo "  --up-gene <int>                Coverage threshold for well-covered gene regions (default: 30)"
    echo "  --down-gene <int>              Coverage threshold for poor-covered gene regions (default: 20)"
    echo "  --well <float>                 Percentage of well-covered regions (default: 95.0)"
    echo "  --poor <float>                 Percentage of poorly-covered regions (default: 5.0)"
    echo "  --chrom <string>               Only genotype regions on this chromosome"
    echo "  -h, --help                     Show this help message"
    exit 1
}

usage_variant_calling() {
    echo "3. Variant Calling QC:"
    echo "The quality control of variant calling"
    echo ""
    echo "Usage: TargetQC.sh --mode variant-calling --seq ES|WGS --vcf <file.vcf> --gvcf <file.gvcf> [--capture <capture.bed>|--ref-bed <ref_bed>] --ref-vcf <ref_vcf> --ref-fasta <ref_fasta> --pathogenic <pathogenic.txt> --out <outprefix> [options]"
    echo ""
    echo "Required parameters:"
    echo ""
    echo "  --mode variant-calling         Run TargetQC in Exon capture QC mode. This mode should be used when performing quality statistics of variant calling analysis"
    echo "  --seq <ES|WGS>                 Sequencing type, ES for exome sequencing, WGS for whole genome sequencing"
    echo "  --vcf <vcf_file>               Input VCF files of sample"
    echo "  --gvcf <gvcf_file>             Input GVCF files of sample"
    echo "  --capture <capture.bed>        Capture regions (.bed). Required for ES sequencing type"
    echo "  --ref-bed <ref_bed>            Benchmark BED file. Required for WGS sequencing type"
    echo "  --ref-vcf <ref_vcf>            Reference dataset VCF file (required when --sample-type is standard)"
    echo "  --ref-fasta <ref_fasta>        Human reference genome fasta file"
    echo "  --pathogenic <pathogenic_file> Target variant sites"
    echo "  --out <outprefix>              Output prefix"
    echo ""
    echo "Additional general options:"
    echo ""
    echo "  --regions <regions.bed>        Additional regions for variant calling benchmark (optional)"
    echo "  -h, --help                     Show this help message"
    exit 1
}

usage_variant_capture() {
    echo "4. Variant Capture QC:"
    echo "Usage: TargetQC.sh --mode variant-capture --seq ES|WGS --vcf <file.vcf> [--capture <capture.bed>] --pathogenic <pathogenic.txt> --out <outprefix> [options]"
    echo ""
    echo "Required parameters:"
    echo ""
    echo "  --mode variant-capture         Run TargetQC in Exon capture QC mode. This mode should be used when performing quality statistics of variant capture analysis"
    echo "  --seq <ES|WGS>                 Sequencing type, ES for exome sequencing, WGS for whole genome sequencing"
    echo "  --vcf <vcf_file>               Input VCF files of sample"
    echo "  --capture <capture.bed>        Capture regions (.bed). Required for ES sequencing type"
    echo "  --pathogenic <pathogenic_file> Target variant sites"
    echo "  --out <outprefix>              Output prefix"
    echo ""
    echo "Optional parameters:"
    echo ""
    echo "  --regions <regions.bed>        Additional regions for variant QC(optional)"
    echo "  -h, --help                     Show this help message"
    exit 1
}

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)
            BAM_FILE="$2"
            shift 2
            ;;
        --capture)
            TARGET_REGIONS_BED="$2"
            shift 2
            ;;
        --regions)
            USER_REGIONS="$2"
            shift 2
            ;;
        --seq)
            SEQ_TYPE="$2"
            shift 2
            ;;
        --out)
            OUT="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --thresholds)
            THRESHOLDS="$2"
            shift 2
            ;;
        --chrom)
            CHROM="$2"
            shift 2
            ;;
        --up-exon)
            UP_EXON="$2"
            shift 2
            ;;
        --down-exon)
            DOWN_EXON="$2"
            shift 2
            ;;
        --base-exon)
            BASE_EXON="$2"
            shift 2
            ;;
        --upstream)
            UPSTREAM="$2"
            shift 2
            ;;
        --downstream)
            DOWNSTREAM="$2"
            shift 2
            ;;
        --up-gene)
            UP_GENE="$2"
            shift 2
            ;;
        --down-gene)
            DOWN_GENE="$2"
            shift 2
            ;;
        --well)
            WELL="$2"
            shift 2
            ;;
        --poor)
            POOR="$2"
            shift 2
            ;;
        --gene)
            GENE="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --sample-type)
            SAMPLE_TYPE="$2"
            shift 2
            ;;
        --vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        --gvcf)
            GVCF_FILE="$2"
            shift 2
            ;;
        --ref-bed)
            REF_BED="$2"
            shift 2
            ;;
        --ref-vcf)
            REF_VCF="$2"
            shift 2
            ;;
        --ref-fasta)
            REF_FASTA="$2"
            shift 2
            ;;
        --pathogenic)
            PATHOGENIC="$2"
            shift 2
            ;;
        --version)
            show_version
            ;;
        -h|--help)
            case $MODE in
                "exon")
                    usage_exon
                    ;;
                "gene")
                    usage_gene
                    ;;
                "variant-calling")
                    usage_variant_calling
                    ;;
                "variant-capture")
                    usage_variant_capture
                    ;;
                *)
                    usage
                    ;;
            esac
            ;;
        *)
            echo "Unknow: $1"
            usage
            ;;
    esac
done

# Check required --seq parameter
if [[ -z "$SEQ_TYPE" ]]; then
    echo "Error: --seq parameter is required!"
    echo "Use TargetQC.sh -h for detailed help information"
    exit 1
fi

if [[ "$SEQ_TYPE" != "ES" && "$SEQ_TYPE" != "WGS" ]]; then
    echo "Error: --seq parameter must be either ES or WGS!"
    echo "Use TargetQC.sh -h for detailed help information"
    exit 1
fi

# Check capture region file for ES sequencing type
if [[ "$SEQ_TYPE" == "ES" && -z "$TARGET_REGIONS_BED" ]]; then
    echo "Error: ES sequencing type requires --capture parameter!"
    echo "Use TargetQC.sh -h for detailed help information"
    exit 1
fi

if [[ -z "$MODE" ]]; then
    MODE="both"
    echo "--mode parameter not specified, defaulting to full flow analysis"
fi

# 根据mode检查参数
case $MODE in
    "exon")
        if [[ -z "$BAM_FILE" || -z "$OUT" || -z "$GENE" || -z "$REFERENCE" ]]; then
            echo "Error: exon mode requires --bam, --gene, --reference and --out parameters!"
            echo "Use TargetQC.sh --mode exon -h for detailed help information"
            exit 1
        fi
        ;;
    "gene")
        if [[ -z "$BAM_FILE" || -z "$OUT" || -z "$GENE" || -z "$REFERENCE" ]]; then
            echo "Error: gene mode requires --bam, --gene, --reference and --out parameters!"
            echo "Use TargetQC.sh --mode gene -h for detailed help information"
            exit 1
        fi
        ;;
    "variant-calling")
        if [[ -z "$VCF_FILE" || -z "$GVCF_FILE" || -z "$OUT" || -z "$REF_VCF" || -z "$REF_FASTA" || -z "$PATHOGENIC" ]]; then
            echo "Error: variant-calling mode requires --vcf, --gvcf, --ref-vcf, --ref-fasta, --pathogenic and --out parameters!"
            echo "Use TargetQC.sh --mode variant-calling -h for detailed help information"
            exit 1
        fi
        if [[ "$SEQ_TYPE" == "WGS" && -z "$REF_BED" ]]; then
            echo "Error: WGS sequencing type requires --ref-bed parameter in variant-calling mode!"
            echo "Use TargetQC.sh --mode variant-calling -h for detailed help information"
            exit 1
        fi
        ;;
    "variant-capture")
        if [[ -z "$VCF_FILE" || -z "$OUT" || -z "$PATHOGENIC" ]]; then
            echo "Error: variant-capture mode requires --vcf, --pathogenic and --out parameters!"
            echo "Use TargetQC.sh --mode variant-capture -h for detailed help information"
            exit 1
        fi
        ;;
    "both")
        if [[ -z "$SAMPLE_TYPE" || -z "$BAM_FILE" || -z "$VCF_FILE" || -z "$OUT" || -z "$GENE" || -z "$REFERENCE" || -z "$PATHOGENIC" ]]; then
            echo "Error: full flow analysis mode requires --sample-type, --bam, --vcf, --gene, --reference, --pathogenic and --out parameters!"
            echo "Use TargetQC.sh -h for detailed help information"
            exit 1
        fi
        if [[ "$SAMPLE_TYPE" == "standard" && ( -z "$GVCF_FILE" || -z "$REF_VCF" || -z "$REF_FASTA" ) ]]; then
            echo "Error: standard sample type requires --gvcf, --ref-vcf and --ref-fasta parameters!"
            echo "Use TargetQC.sh -h for detailed help information"
            exit 1
        fi
        if [[ "$SAMPLE_TYPE" == "standard" && "$SEQ_TYPE" == "WGS" && -z "$REF_BED" ]]; then
            echo "Error: WGS sequencing type requires --ref-bed parameter for standard sample type!"
            echo "Use TargetQC.sh -h for detailed help information"
            exit 1
        fi
        ;;
    *)
        echo "Error: Unknown mode parameter: $MODE"
        echo "Use TargetQC.sh -h for detailed help information"
        exit 1
        ;;
esac

TEMP_DIR="${OUT}/temp"
OUTPUT_DIR="${OUT}"
OUTPREFIX=$(basename "$OUT")
mkdir -p "$TEMP_DIR"
mkdir -p "$OUTPUT_DIR"

echo "Starting TargetQC workflow..."
if [[ "$MODE" != "both" ]]; then
    echo "Analysis mode: $MODE"
fi
echo "Sequencing type: $SEQ_TYPE"
echo "Output prefix: $OUTPREFIX"

if [[ "$MODE" == "exon" || "$MODE" == "gene" || "$MODE" == "both" ]]; then
    if [[ ! -z "$TARGET_REGIONS_BED" ]]; then
        echo "Capture regions file: $TARGET_REGIONS_BED"
    fi
    if [[ ! -z "$USER_REGIONS" ]]; then
        echo "User regions file: $USER_REGIONS"
    fi
    if [[ ! -z "$GENE" ]]; then
      echo "Target genes file: $GENE"
    fi
    if [[ ! -z "$REFERENCE" ]]; then
      echo "Gene annotation file: $REFERENCE"
    fi
fi

if [[ "$MODE" == "variant-calling" || "$MODE" == "variant-capture" || "$MODE" == "both" ]]; then
    echo "Sample VCF file: $VCF_FILE"
    if [[ "$MODE" == "variant-calling" || ("$MODE" == "both" && "$SAMPLE_TYPE" == "standard") ]]; then
        if [[ ! -z "$REF_BED" ]]; then
            echo "Reference BED: $REF_BED"
        fi
        echo "Reference VCF: $REF_VCF"
        echo "Reference FASTA: $REF_FASTA"
    fi
fi

CHROM_PARAM=""
if [[ ! -z "$CHROM" ]]; then
    CHROM_PARAM="--chrom $CHROM"
fi

if [[ "$MODE" == "exon" || "$MODE" == "both" ]]; then
    echo "Upstream extension: $UPSTREAM bp"
    echo "Downstream extension: $DOWNSTREAM bp"
    echo "Coverage threshold for well-covered: $UP_EXON"
    echo "Coverage threshold for poor-covered: $DOWN_EXON"
    echo "Number of poor-covered bases in exons: $BASE_EXON"
    echo "Executing exon capture QC analysis..."

    # Step 1: Generate target exon region files
    echo "Step 1: Generating capture region files..."

    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next}
        $2=="exon" && $6 in transcript && $8==transcript[$6] {
            print $1 "\t" $3-1 "\t" $4
        }' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_exon.bed"

    # Step 2: Execute coverage analysis
    echo "Step 2: Executing exon coverage calculation analysis..."

    $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_EXON,$UP_EXON" --out "${TEMP_DIR}/${OUTPREFIX}_exon"
    $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_exon.sh"

    # Only execute capture region analysis for ES type
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "$TARGET_REGIONS_BED" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture"
        $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_capture.sh"
    fi

    # Step 3: Execute exon capture QC analysis workflow
    echo "Step 3: Executing exon coverage analysis workflow..."
    if [[ ! -z "$CHROM" ]]; then
        if [[ "$SEQ_TYPE" == "ES" ]]; then
            $PYTHON3_PATH "${SCRIPT_DIR}/src/base.py" --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture.regions.bed.gz" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
        else
            $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture_exon"
            $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_capture_exon.sh"
            $PYTHON3_PATH "${SCRIPT_DIR}/src/base.py" --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.regions.bed.gz" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
        fi
        $PYTHON3_PATH "${SCRIPT_DIR}/src/exon_coverage.py" --thresholds-exon "${TEMP_DIR}/${OUTPREFIX}_exon.thresholds.bed.gz" --regions-exon "${TEMP_DIR}/${OUTPREFIX}_exon.regions.bed.gz" --up-exon "$UP_EXON" --down-exon "$DOWN_EXON" --base-exon "$BASE_EXON" --chrom "$CHROM" --out "${TEMP_DIR}/${OUTPREFIX}"
    else
        if [[ "$SEQ_TYPE" == "ES" ]]; then
            $PYTHON3_PATH "${SCRIPT_DIR}/src/base.py" --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture.regions.bed.gz" --out "${TEMP_DIR}/${OUTPREFIX}"
        else
            $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$THRESHOLDS" --out "${TEMP_DIR}/${OUTPREFIX}_capture_exon"
            $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_capture_exon.sh"
            $PYTHON3_PATH "${SCRIPT_DIR}/src/base.py" --thresholds-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.thresholds.bed.gz" --regions-base "${TEMP_DIR}/${OUTPREFIX}_capture_exon.regions.bed.gz" --out "${TEMP_DIR}/${OUTPREFIX}"
        fi
        $PYTHON3_PATH "${SCRIPT_DIR}/src/exon_coverage.py" --thresholds-exon "${TEMP_DIR}/${OUTPREFIX}_exon.thresholds.bed.gz" --regions-exon "${TEMP_DIR}/${OUTPREFIX}_exon.regions.bed.gz" --up-exon "$UP_EXON" --down-exon "$DOWN_EXON" --base-exon "$BASE_EXON" --out "${TEMP_DIR}/${OUTPREFIX}"
    fi

    # Only execute exon capture QC for ES type
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        $BASH_PATH "${SCRIPT_DIR}/src/exon_capture_QC.sh" --coverage-exon "${TEMP_DIR}/${OUTPREFIX}_exon_coverage.tsv" --regions "$TARGET_REGIONS_BED" --out "${TEMP_DIR}/${OUTPREFIX}"
    else
        awk -F'\t' 'BEGIN{OFS="\t"} 
        NR==1 {print "#chrom", "start", "end", "coverage", "type", "capture"; next}
        {print $1, $2, $3, $4, $5, "full_capture"}' "${TEMP_DIR}/${OUTPREFIX}_exon_coverage.tsv" > "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv"
    fi

    $PYTHON3_PATH "${SCRIPT_DIR}/src/exon_summary.py" --qc-exon "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv" --out "${TEMP_DIR}/${OUTPREFIX}"
fi

if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    echo "Executing gene capture QC analysis..."

    # Step 1: Generate target exon and CDS region files
    echo "Step 1: Generating target exon and CDS region files..."

    # Generate exon region file
    echo "Generating target exon region file..."
    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next}
        $2=="exon" && $6 in transcript && $8==transcript[$6] {
            print $1 "\t" $3-1 "\t" $4
        }' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_exon.bed"

    # Generate CDS region file
    echo "Generating target CDS region file..."
    awk -F'\t' 'NR==FNR{transcript[$1]=$2; next}
        $2=="CDS" && $6 in transcript && $8==transcript[$6] {
            print $1 "\t" $3-1 "\t" $4
        }' "$GENE" "$REFERENCE" > "${TEMP_DIR}/${OUTPREFIX}_cds.bed"

    # Step 2: Execute coverage analysis
    echo "Step 2: Executing coverage analysis..."

    # Execute exon region analysis
    echo "Executing target exon region coverage analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_exon.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_GENE,$UP_GENE" --out "${TEMP_DIR}/${OUTPREFIX}_target_exon"
    $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_target_exon.sh"

    # Execute CDS region analysis
    echo "Executing target CDS region coverage analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/target_coverage.py" --bam "$BAM_FILE" --regions "${TEMP_DIR}/${OUTPREFIX}_cds.bed" --mosdepth-path "$MOSDEPTH" --thresholds "$DOWN_GENE,$UP_GENE" --out "${TEMP_DIR}/${OUTPREFIX}_target_cds"
    $BASH_PATH "${TEMP_DIR}/${OUTPREFIX}_target_cds.sh"

    # Step 3: Execute transcript and CDS analysis
    echo "Step 3: Executing target transcript and CDS analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/trans.py" --coverage "${TEMP_DIR}/${OUTPREFIX}_target_exon.thresholds.bed.gz" --gene "$GENE" --reference "$REFERENCE" --up-gene "$UP_GENE" --down-gene "$DOWN_GENE" --out "${TEMP_DIR}/${OUTPREFIX}"
    $PYTHON3_PATH "${SCRIPT_DIR}/src/cds.py" --coverage "${TEMP_DIR}/${OUTPREFIX}_target_cds.thresholds.bed.gz" --gene "$GENE" --reference "$REFERENCE" --up "$UP_GENE" --down "$DOWN_GENE" --out "${TEMP_DIR}/${OUTPREFIX}"

    awk -F'\t' '$3=="protein_coding"' "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" > "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon.tsv"

    # Step 4: Execute proportion analysis
    echo "Step 4: Executing proportion analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/trans_proportion.py" --gene "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans"
    $PYTHON3_PATH "${SCRIPT_DIR}/src/trans_proportion.py" --gene "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans_pc_exon"
    $PYTHON3_PATH "${SCRIPT_DIR}/src/trans_proportion.py" --gene "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" --out "${TEMP_DIR}/${OUTPREFIX}_trans_cds"

    # Step 5: Execute gene capture QC analysis
    echo "Step 5: Executing gene capture QC analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/gene_capture_QC.py" --trans-exon "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" --trans-cds "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" --well "$WELL" --poor "$POOR" --out "${TEMP_DIR}/${OUTPREFIX}"
    $PYTHON3_PATH "${SCRIPT_DIR}/src/gene_capture_QC_summary.py" --gene-qc "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.tsv" --out "${TEMP_DIR}/${OUTPREFIX}"
fi

if [[ "$MODE" == "variant-calling" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "standard" ) ]]; then
    echo "Executing variant Calling QC analysis..."

    # Benchmark analysis
    echo "Step 1: Executing Benchmark analysis..."
    
    # Use appropriate regions file based on sequencing type
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        REGIONS_FILE="$TARGET_REGIONS_BED"
    else
        # For WGS, use exon regions file
        REGIONS_FILE="$REF_BED"
    fi
    
    # Execute default benchmark analysis
    $PYTHON3_PATH "${SCRIPT_DIR}/src/benchmark.py" --vcf "$VCF_FILE" --regions "$REGIONS_FILE" --ref-vcf "$REF_VCF" --ref-fasta "$REF_FASTA" --tools "$HAP_PATH" --python2-path "$PYTHON2_PATH"  --out "${TEMP_DIR}/${OUTPREFIX}_capture"
    $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_summary.py" --input "${TEMP_DIR}/${OUTPREFIX}_capture.summary.csv" --out "${TEMP_DIR}/${OUTPREFIX}_variant_calling"
    
    # Execute additional benchmark analysis if user regions are provided
    if [[ ! -z "$USER_REGIONS" ]]; then
        echo "Executing additional benchmark analysis with user regions..."
        $PYTHON3_PATH "${SCRIPT_DIR}/src/benchmark.py" --vcf "$VCF_FILE" --regions "$USER_REGIONS" --ref-vcf "$REF_VCF" --ref-fasta "$REF_FASTA" --tools "$HAP_PATH" --python2-path "$PYTHON2_PATH" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_summary.py" --input "${TEMP_DIR}/${OUTPREFIX}_regions.summary.csv" --out "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling"
    fi
    
    # Step 2: Execute coverage analysis
    echo "Step 2: Executing various variant analysis..."
        # Use default capture VCF for analysis
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$GVCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "FN" --out "${TEMP_DIR}/${OUTPREFIX}"
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "FP" --out "${TEMP_DIR}/${OUTPREFIX}"
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_capture.vcf.gz" --type "TP" --out "${TEMP_DIR}/${OUTPREFIX}"
    {
#        echo -e "#chrom\tpos\tvar_type\tgenotype\tAD_ref\tAD_alt\tDP\tGQ\tAD_DP\tBAF"
#        echo -e "#chrom\tpos\ttype\tvariant_type\tgenotype\tDP\tGQ\tBAF"
        echo -e "#chrom\tpos\ttype\tvariant_type\tgenotype\tDP\tBAF"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "TP", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_TP_QC.tsv"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "FP", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_FP_QC.tsv"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "FN", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_FN_QC.tsv"
    } > "${TEMP_DIR}/${OUTPREFIX}_variant_calling_QC.tsv"

    if [[ ! -z "$USER_REGIONS" ]]; then
        # Use user regions VCF for analysis
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$GVCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "FN" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "FP" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_calling_QC.py" --vcf-orig "$VCF_FILE" --vcf-process "${TEMP_DIR}/${OUTPREFIX}_regions.vcf.gz" --type "TP" --out "${TEMP_DIR}/${OUTPREFIX}_regions"
    {
#        echo -e "#chrom\tpos\tvar_type\tgenotype\tAD_ref\tAD_alt\tDP\tGQ\tAD_DP\tBAF"
        echo -e "#chrom\tpos\ttype\tvariant_type\tgenotype\tDP\tBAF"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "TP", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_regions_TP_QC.tsv"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "FP", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_regions_FP_QC.tsv"
        awk -F'\t' -v OFS='\t' 'NR>1 {print $1, $2, "FN", $3, $4, $5, $6}' "${TEMP_DIR}/${OUTPREFIX}_regions_FN_QC.tsv"
    } > "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling_QC.tsv"

    fi

    # Step 3: Pathogenicity analysis
    echo "Step 3: Pathogenicity analysis..."
    $PYTHON3_PATH "${SCRIPT_DIR}/src/pathogenic.py" --fn "${TEMP_DIR}/${OUTPREFIX}_FN_QC.tsv" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"

    if [[ ! -z "$USER_REGIONS" ]]; then
        # Use user regions VCF for analysis
        $PYTHON3_PATH "${SCRIPT_DIR}/src/pathogenic.py" --fn "${TEMP_DIR}/${OUTPREFIX}_regions_FN_QC.tsv" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}_region"
    fi
fi

if [[ "$MODE" == "variant-capture" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "clinical")  ]]; then
    echo "Executing variant Calling QC analysis..."

    # Use appropriate regions file based on sequencing type
    if [[ "$SEQ_TYPE" == "ES" ]]; then
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_capture_QC.py" --vcf "$VCF_FILE" --regions "$TARGET_REGIONS_BED" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"
    else
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_capture_QC.py" --vcf "$VCF_FILE" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}"        
    fi
    if [[ ! -z "$USER_REGIONS" ]]; then
        # Use user regions VCF for analysis
        $PYTHON3_PATH "${SCRIPT_DIR}/src/variant_capture_QC.py" --vcf "$VCF_FILE" --regions "$USER_REGIONS" --pathogenic "$PATHOGENIC" --out "${TEMP_DIR}/${OUTPREFIX}_region"
    fi
fi
    
# Copy main output files to output directory
echo "Copying main output files to ${OUTPUT_DIR}..."

if [[ "$MODE" == "exon" || "$MODE" == "both" ]]; then
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_base.tsv" ]]; then
#        cp "${TEMP_DIR}/${OUTPREFIX}_base.tsv" "${OUT}_base.tsv"
        cp "${TEMP_DIR}/${OUTPREFIX}_base.tsv" "${OUT}/${OUTPREFIX}_base.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.tsv" "${OUT}/${OUTPREFIX}_exon_capture_QC.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.summary.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_exon_capture_QC.summary.tsv" "${OUT}/${OUTPREFIX}_exon_capture_QC.summary.tsv"
    fi
fi

if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_trans.tsv" "${OUT}/${OUTPREFIX}_trans.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_trans_cds.tsv" "${OUT}/${OUTPREFIX}_trans_cds.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.tsv" "${OUT}/${OUTPREFIX}_gene_capture_QC.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.summary.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_gene_capture_QC.summary.tsv" "${OUT}/${OUTPREFIX}_gene_capture_QC.summary.tsv"
    fi
fi

if [[ "$MODE" == "variant-calling" || "$MODE" == "both" ]]; then
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_variant_calling.summary.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_variant_calling.summary.tsv" "${OUT}/${OUTPREFIX}_variant_calling.summary.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling.summary.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling.summary.tsv" "${OUT}/${OUTPREFIX}_region_variant_calling.summary.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_variant_calling_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_variant_calling_QC.tsv" "${OUT}/${OUTPREFIX}_variant_calling_QC.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_region_variant_calling_QC.tsv" "${OUT}/${OUTPREFIX}_region_variant_calling_QC.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_pathogenic.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_pathogenic.tsv" "${OUT}/${OUTPREFIX}_pathogenic.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_region_pathogenic.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_region_pathogenic.tsv" "${OUT}/${OUTPREFIX}_region_pathogenic.tsv"
    fi
fi

if [[ "$MODE" == "variant-capture" || "$MODE" == "both" ]]; then
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_variant_capture_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_variant_capture_QC.tsv" "${OUT}/${OUTPREFIX}_variant_capture_QC.tsv"
    fi
    if [[ -f "${TEMP_DIR}/${OUTPREFIX}_region_variant_capture_QC.tsv" ]]; then
        cp "${TEMP_DIR}/${OUTPREFIX}_region_variant_capture_QC.tsv" "${OUT}/${OUTPREFIX}_region_variant_capture_QC.tsv"
    fi
fi

# Remove temporary directory
rm -rf "$TEMP_DIR"
#rm -rf "$OUT"

echo "TargetQC workflow completed!"
echo ""
echo "Main output files:"

if [[ "$MODE" == "exon" || "$MODE" == "both" ]]; then

    echo "${OUT}/${OUTPREFIX}_base.tsv"
    echo "${OUT}/${OUTPREFIX}_exon_capture_QC.tsv"
    echo "${OUT}/${OUTPREFIX}_exon_capture_QC.summary.tsv"
fi

if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then

    echo "${OUT}/${OUTPREFIX}_trans.tsv"
    echo "${OUT}/${OUTPREFIX}_trans_cds.tsv"
    echo "${OUT}/${OUTPREFIX}_gene_capture_QC.tsv"
    echo "${OUT}/${OUTPREFIX}_gene_capture_QC.summary.tsv"
fi

if [[ "$MODE" == "variant-calling" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "standard" ) ]]; then
    echo "${OUT}/${OUTPREFIX}_variant_calling.summary.tsv"
    echo "${OUT}/${OUTPREFIX}_variant_calling_QC.tsv"
    if [[ -f "${OUT}/${OUTPREFIX}_pathogenic.tsv" ]] && [[ $(wc -l < "${OUT}/${OUTPREFIX}_pathogenic.tsv") -gt 1 ]]; then
        echo "${OUT}/${OUTPREFIX}_pathogenic.tsv"
    if [[ ! -z "$USER_REGIONS" ]]; then
        echo "${OUT}/${OUTPREFIX}_region_variant_calling.summary.tsv"
        echo "${OUT}/${OUTPREFIX}_region_variant_calling_QC.tsv"
        if [[ -f "${OUT}/${OUTPREFIX}_region_pathogenic.tsv" ]] && [[ $(wc -l < "${OUT}/${OUTPREFIX}_region_pathogenic.tsv") -gt 1 ]]; then
            echo "${OUT}/${OUTPREFIX}_region_pathogenic.tsv"
    fi
fi

if [[ "$MODE" == "variant-capture" || ( "$MODE" == "both" && "$SAMPLE_TYPE" == "clinical" )  ]]; then
    echo "${OUT}/${OUTPREFIX}_variant_capture_QC.tsv"
    if [[ ! -z "$USER_REGIONS" ]]; then
        echo "${OUT}/${OUTPREFIX}_region_variant_capture_QC.tsv"
    fi
fi

