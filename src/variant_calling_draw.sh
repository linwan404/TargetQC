#!/bin/bash
# variant_calling_draw.sh - Extract variant calls for visualization

usage() {
    cat << EOF
usage:  variant_calling_draw.sh --vcf <file.vcf> --out <outprefix>

Extract variant calls for visualization

optional arguments: 
  -h, --help            show this help message and exit

Required options:
--vcf <file.vcf>  Input VCF files of the sample
--out        Prefix to name output files
EOF
}

# 参数解析
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)
            VCF_FILE="$2"
            shift 2
            ;;
        --out)
            OUT_PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$VCF_FILE" || -z "$OUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    usage
    exit 1
fi

# 检查输入文件是否存在
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file $VCF_FILE does not exist"
    exit 1
fi

# 创建输出目录
OUT_DIR=$(dirname "$OUT_PREFIX")
if [[ -n "$OUT_DIR" && ! -d "$OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
    echo "Created output directory: $OUT_DIR"
fi

# 定义输出文件
FN_FILE="${OUT_PREFIX}_FN.tsv"
FP_FILE="${OUT_PREFIX}_FP.tsv"
TP_FILE="${OUT_PREFIX}_TP.tsv"

echo "Processing VCF file: $VCF_FILE"
echo "Output files:"
echo "  FN: $FN_FILE"
echo "  FP: $FP_FILE"
echo "  TP: $TP_FILE"

# 使用Awk处理VCF内容
if [[ "$VCF_FILE" == *.gz ]]; then
    # 如果是gzip压缩文件，使用zcat
    DECOMPRESS_CMD="zcat"
else
    # 如果不是压缩文件，使用cat
    DECOMPRESS_CMD="cat"
fi

$DECOMPRESS_CMD "$VCF_FILE" | awk -v fn_file="$FN_FILE" -v fp_file="$FP_FILE" -v tp_file="$TP_FILE" '
BEGIN {
    header = "#chrom\tpos\tref\talt";
    
    # 初始化文件头
    printf("%s\n", header) > fn_file;
    printf("%s\n", header) > fp_file;
    printf("%s\n", header) > tp_file;
    
    # 初始化计数器
    fn_count = 0;
    fp_count = 0;
    tp_count = 0;
}
/^##/ { next; }  # 跳过注释行
/^#CHROM/ { next; }  # 跳过头行
{
    # 处理第11列检查FP类型
    split($11, fp_parts, ":");
    if (fp_parts[2] == "FP") {
        fp_count++;
        printf("%s\t%s\t%s\t%s\n", $1, $2, $4, $5) >> fp_file;
    }

    # 处理第10列检查FN和TP类型
    split($10, parts, ":");
    type_10 = parts[2];
    if (type_10 == "FN") {
        fn_count++;
        printf("%s\t%s\t%s\t%s\n", $1, $2, $4, $5) >> fn_file;
    } else if (type_10 == "TP") {
        tp_count++;
        printf("%s\t%s\t%s\t%s\n", $1, $2, $4, $5) >> tp_file;
    }
}
END {
    # 输出统计结果
    print "Processing completed:";
    print "  FP variants:", fp_count;
    print "  FN variants:", fn_count;
    print "  TP variants:", tp_count;
    print "  Total variants processed:", fp_count + fn_count + tp_count;
}'

echo "Variant extraction completed"
echo "Output files:"
echo "  $FN_FILE"
echo "  $FP_FILE"
echo "  $TP_FILE"
