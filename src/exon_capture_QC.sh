#!/bin/bash
# exon_capture_QC.sh - Process exon coverage with capture region analysis

usage() {
    cat << EOF
usage:  add_capture.sh --coverage-exon <coverage.tsv> --regions <regions.bed> --out <outprefix>

optional arguments:
  -h, --help            show this help message and exit

Required options:
  --coverage-exon    <coverage.tsv>             Exon coverage input file of the sample
  --regions    <regions.bed>             Capture regions files of the sample
  --out         <outprefix>       Prefix to name output files
EOF
}
UPSTREAM=150
DOWNSTREAM=150
# 参数解析
while [[ $# -gt 0 ]]; do
    case $1 in
        --coverage-exon)
            COVERAGE_EXON="$2"
            shift 2
            ;;
        --regions)
            REGIONS="$2"
            shift 2
            ;;
        --out)
            OUT_PREFIX="$2"
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
if [[ -z "$COVERAGE_EXON" || -z "$REGIONS" || -z "$OUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    usage
    exit 1
fi

# 检查输入文件是否存在
if [[ ! -f "$COVERAGE_EXON" ]]; then
    echo "Error: Coverage file $COVERAGE_EXON does not exist"
    exit 1
fi

if [[ ! -f "$REGIONS" ]]; then
    echo "Error: Regions file $REGIONS does not exist"
    exit 1
fi

# 设置文件路径
file_a="$COVERAGE_EXON"
file_b="$REGIONS"
output_file="${OUT_PREFIX}_exon_capture_QC.tsv"

# 从输出前缀提取目录路径和样本名
output_dir=$(dirname "$OUT_PREFIX")
sample_name=$(basename "$OUT_PREFIX")

# 创建临时目录
temp_dir="${output_dir}/temp"
mkdir -p "$temp_dir"

echo "Processing sample: $sample_name"
echo "Coverage file: $file_a"
echo "Regions file: $file_b"
echo "Output file: $output_file"
echo "Temp directory: $temp_dir"

# Step 1: 处理文件a并排序
tail -n +2 "$file_a" | sort -k1,1V -k2,2n > "${temp_dir}/${sample_name}_a_processed.tsv"

# Step 2: 处理原始和扩展后的bed文件
# 原始bed文件处理
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3}' "$file_b" | sort -k1,1V -k2,2n > "${temp_dir}/${sample_name}_b_original.tsv"

# 扩展bed文件（起始-150，终止+150）

awk -v upstream='"$UPSTREAM"' -v downstream='"$DOWNSTREAM"' 'BEGIN {FS=OFS="\t"} {
    $2 = ($2 - upstream < 0) ? 0 : $2 - upstream;
    $3 = $3 + downstream;
    print $1, $2, $3
}' "$file_b" | sort -k1,1V -k2,2n > "${temp_dir}/${sample_name}_b_extended.tsv"

# Step 3: 初始交集（原始bed）
bedtools intersect -a "${temp_dir}/${sample_name}_a_processed.tsv" -b "${temp_dir}/${sample_name}_b_original.tsv" -wo > "${temp_dir}/${sample_name}_overlap_original.tsv"

# Step 4: 生成初始结果（包含all/part/no_capture）
awk '
BEGIN {
    FS = OFS = "\t";
}
{
    key = $1 "\t" $2 "\t" $3;
    overlap_length[key] += $9;
    exon_info[key] = $1 OFS $2 OFS $3 OFS $4 OFS $5;
}
END {
    for (key in overlap_length) {
        split(key, arr, "\t");
        chrom = arr[1];
        start = arr[2];
        end = arr[3];
        region_length = end - start;
        raw_percent = (overlap_length[key] / region_length) * 100.00;

        split(exon_info[key], fields, "\t");
        capture = (raw_percent >= 100.00) ? "full_capture" : (raw_percent > 0) ? "partial_capture" : "no_capture";
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n", chrom, start, end, fields[4], fields[5], capture, raw_percent;
    }
}' "${temp_dir}/${sample_name}_overlap_original.tsv" > "${temp_dir}/${sample_name}_initial_result.tsv"

# Step 5: 提取no_capture区域
# 生成a_processed的键列表
cut -f1-3 "${temp_dir}/${sample_name}_a_processed.tsv" | sort > "${temp_dir}/${sample_name}_a_keys.tsv"
# 生成initial_result的键列表
cut -f1-3 "${temp_dir}/${sample_name}_initial_result.tsv" | sort > "${temp_dir}/${sample_name}_initial_keys.tsv"
# 使用comm找出差异
comm -23 "${temp_dir}/${sample_name}_a_keys.tsv" "${temp_dir}/${sample_name}_initial_keys.tsv" \
    > "${temp_dir}/${sample_name}_no_capture_keys.tsv"

# 从a_processed.tsv中提取完整的no_capture行
awk -v OFS="\t" '
    NR==FNR {
        # 将前三列拼接为字符串作为键
        key = $1 OFS $2 OFS $3
        keys[key] = 1
        next
    }
    {
        # 检查当前行前三列是否在keys中
        current_key = $1 OFS $2 OFS $3
        if (current_key in keys) {
            print $0
        }
    }
' "${temp_dir}/${sample_name}_no_capture_keys.tsv" "${temp_dir}/${sample_name}_a_processed.tsv" > "${temp_dir}/${sample_name}_no_capture.bed"

# Step 6: 与扩展bed文件重新交集
bedtools intersect -a "${temp_dir}/${sample_name}_no_capture.bed" -b "${temp_dir}/${sample_name}_b_extended.tsv" -wo > "${temp_dir}/${sample_name}_overlap_extended.tsv"

# Step 7: 计算邻近捕获（near_capture）
awk '
BEGIN {
    FS = OFS = "\t";
}
{
    chrom = $1;
    a_start = $2;
    a_end = $3;
    coverage = $4;
    type = $5;
    b_start = $7;
    b_end = $8;

    # 计算交集区域
    overlap_start = (b_start > a_start) ? b_start : a_start;
    overlap_end = (b_end < a_end) ? b_end : a_end;
    if (overlap_start >= overlap_end) next;

    key = chrom "\t" a_start "\t" a_end;
    # 保存交集区间
    starts[key] = (starts[key] ? starts[key] " " : "") overlap_start;
    ends[key] = (ends[key] ? ends[key] " " : "") overlap_end;

    # 记录exon信息
    if (!(key in exon_info)) {
        exon_info[key] = chrom OFS a_start OFS a_end OFS coverage OFS type;
    }
}
END {
    for (key in exon_info) {
        split(key, parts, "\t");
        chrom = parts[1];
        a_start = parts[2] + 0;
        a_end = parts[3] + 0;
        exon_length = a_end - a_start;

        total_length = 0;
        split(starts[key], s, " ");
        split(ends[key], e, " ");
        n = length(s);

        if (n > 0) {
            # 对区间进行排序
            for (i = 1; i <= n; i++) {
                indices[i] = i;
            }
            for (i = 1; i <= n; i++) {
                for (j = i+1; j <= n; j++) {
                    if (s[indices[j]] < s[indices[i]]) {
                        temp = indices[i];
                        indices[i] = indices[j];
                        indices[j] = temp;
                    }
                }
            }

            # 合并区间
            merged_starts = s[indices[1]];
            merged_ends = e[indices[1]];
            for (i = 2; i <= n; i++) {
                idx = indices[i];
                current_start = s[idx];
                current_end = e[idx];
                last_end = merged_ends + 0;

                if (current_start <= last_end) {
                    if (current_end > last_end) {
                        merged_ends = current_end;
                    }
                } else {
                    merged_starts = merged_starts " " current_start;
                    merged_ends = merged_ends " " current_end;
                }
            }

            # 计算总长度
            split(merged_starts, ms, " ");
            split(merged_ends, me, " ");
            for (i = 1; i <= length(ms); i++) {
                total_length += me[i] - ms[i];
            }
        }

        # 计算覆盖百分比
        coverage_percent = (exon_length == 0) ? 0 : (total_length / exon_length) * 100.00;
        coverage_percent = sprintf("%.2f", coverage_percent);
        capture = "near_capture";
        print exon_info[key], capture, coverage_percent;
    }
}' "${temp_dir}/${sample_name}_overlap_extended.tsv" > "${temp_dir}/${sample_name}_near_capture.tsv"

# Step 8: 合并结果
awk -v OFS="\t" '
BEGIN {
    # 加载initial_result数据 (all_capture/part_capture)
    while (getline < "'"${temp_dir}/${sample_name}_initial_result.tsv"'" > 0) {
        key = $1 "\t" $2 "\t" $3
        initial[key] = $0  # 存储完整行
    }
    # 加载near_capture数据 (near_capture)
    while (getline < "'"${temp_dir}/${sample_name}_near_capture.tsv"'" > 0) {
        key = $1 "\t" $2 "\t" $3
        capture[key] = $6  # 捕获类型
        percent[key] = $7  # 覆盖百分比
    }
}
# 处理a_processed.tsv中的所有区域
{
    key = $1 "\t" $2 "\t" $3
    if (key in initial) {
        # 输出initial_result中的条目
        print initial[key]
    } else if (key in capture) {
        # 输出near_capture条目
        print $1, $2, $3, $4, $5, capture[key], percent[key]
    } else {
        # 标记为no_capture
        print $1, $2, $3, $4, $5, "no_capture", "0.00"
    }
}
' "${temp_dir}/${sample_name}_a_processed.tsv" \
| sort -k1,1V -k2,2n -k3,3n \
| awk 'BEGIN {print "#chrom\tstart\tend\tcoverage\ttype\tcapture\tcapture_percentage"} 1' > "$output_file"

# 清理临时文件（可选，调试时可注释掉）
 rm -rf "$temp_dir"

echo "Processing completed: $output_file"
