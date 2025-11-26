#!/usr/bin/env python3
"""
variant_capture_QC.py - Generate variant capture quality control metrics
"""

import os
import gzip
import bisect
import argparse
from collections import defaultdict
from pathlib import Path

def load_bed_intervals(bed_path):
    """加载BED文件并转换为1-based的区间字典，键为染色体，值为排序后的区间列表"""
    intervals = defaultdict(list)
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            start_1based = start + 1
            end_1based = end
            intervals[chrom].append((start_1based, end_1based))
    return intervals

def merge_intervals(intervals):
    """合并重叠或相邻的区间"""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def calculate_genotype(gt):
    """根据GT值判断基因型"""
    if gt in ["1|1", "1/1", "1/2", "1|2", "2|2", "2/2"]:
        return "hom"
    elif gt in ["0|1", "0/1", "1|0", "0/2", "2|0"]:
        return "het"
    else:
        return "unknown"

def calculate_baf(genotype, ad_alt, dp):
    """根据基因型和等位基因深度计算BAF"""
    if dp == 0:
        return 0.0

    allele_freq = ad_alt / dp

    if genotype == "hom":
        return abs(1 - allele_freq)
    elif genotype == "het":
        return abs(0.5 - allele_freq)
    else:
        return 0.0

def determine_var_type(ref, alt):
    """根据REF和ALT确定变异类型"""
    if len(ref) == len(alt):
        return "SNV"
    else:
        return "INDEL"

def process_variant_capture_qc(vcf_file, regions_file, pathogenic_file, output_file):
    """处理VCF文件，提取捕获QC指标"""

    # 检查必需的文件是否存在
    if not all(os.path.exists(f) for f in [vcf_file, pathogenic_file]):
        print(f"Error: Required input files do not exist")
        return False

    # 检查可选的regions文件是否存在（如果提供了的话）
    if regions_file and not os.path.exists(regions_file):
        print(f"Error: Regions file {regions_file} does not exist")
        return False

    # 加载致病位点区间
    pathogenic_intervals = load_bed_intervals(pathogenic_file)
    for chrom in pathogenic_intervals:
        pathogenic_intervals[chrom] = merge_intervals(pathogenic_intervals[chrom])

    # 加载捕获区域区间（如果提供了regions文件）
    capture_intervals = {}
    if regions_file:
        capture_intervals = load_bed_intervals(regions_file)
        for chrom in capture_intervals:
            capture_intervals[chrom] = merge_intervals(capture_intervals[chrom])
        print(f"Using capture regions from: {regions_file}")
    else:
        print("No capture regions provided, processing all variants in pathogenic regions")

    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    # 确定VCF文件打开方式
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    results = []

    try:
        with open_func(vcf_file, mode) as vcf:
            for line in vcf:
                if line.startswith('##'):
                    continue

                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    try:
                        format_idx = headers.index('FORMAT')
                        sample_idx = format_idx + 1
                    except ValueError:
                        continue
                    continue

                parts = line.strip().split('\t')
                if len(parts) <= sample_idx:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                format_str = parts[format_idx]
                sample_data = parts[sample_idx]

                # 检查捕获区间（如果提供了regions文件）
                in_capture = True  # 默认值为True，如果没有提供regions文件
                if regions_file:
                    in_capture = False
                    chrom_capture = capture_intervals.get(chrom, [])
                    if chrom_capture:
                        starts = [s for s, _ in chrom_capture]
                        i = bisect.bisect_right(starts, pos) - 1
                        if i >= 0 and chrom_capture[i][0] <= pos <= chrom_capture[i][1]:
                            in_capture = True

                # 检查致病位点区间
                in_pathogenic = False
                chrom_pathogenic = pathogenic_intervals.get(chrom, [])
                if chrom_pathogenic:
                    starts = [s for s, _ in chrom_pathogenic]
                    i = bisect.bisect_right(starts, pos) - 1
                    if i >= 0 and chrom_pathogenic[i][0] <= pos <= chrom_pathogenic[i][1]:
                        in_pathogenic = True

                # 根据是否提供了regions文件调整处理逻辑
                if regions_file:
                    # 如果提供了regions文件，只处理同时在捕获区域和致病区域的位点
                    if in_capture and in_pathogenic:
                        process_variant = True
                    else:
                        process_variant = False
                else:
                    # 如果没有提供regions文件，只处理在致病区域的位点
                    if in_pathogenic:
                        process_variant = True
                    else:
                        process_variant = False

                if process_variant:
                    format_fields = format_str.split(':')
                    sample_fields = sample_data.split(':')

                    # 提取GT
                    gt_idx = format_fields.index('GT') if 'GT' in format_fields else -1
                    gt = sample_fields[gt_idx] if gt_idx != -1 and gt_idx < len(sample_fields) else './.'

                    # 提取DP
                    dp_idx = format_fields.index('DP') if 'DP' in format_fields else -1
                    dp = int(sample_fields[dp_idx]) if dp_idx != -1 and dp_idx < len(sample_fields) and sample_fields[dp_idx] else 0

                    # 提取AD
                    ad_idx = format_fields.index('AD') if 'AD' in format_fields else -1
                    ad_str = sample_fields[ad_idx] if ad_idx != -1 and ad_idx < len(sample_fields) else '0,0'

                    # 解析AD_ref和AD_alt
                    ad_parts = ad_str.split(',')
                    ad_ref = int(ad_parts[0]) if len(ad_parts) > 0 and ad_parts[0] else 0
                    ad_alt = int(ad_parts[1]) if len(ad_parts) > 1 and ad_parts[1] else 0

                    # 计算其他指标
                    var_type = determine_var_type(ref, alt)
                    genotype = calculate_genotype(gt)
                    ad_alt_dp = ad_alt / dp if dp > 0 else 0.0
                    baf = calculate_baf(genotype, ad_alt, dp)

                    results.append((
                        chrom, pos, var_type, genotype,
#                        dp, ad_ref, ad_alt, ad_alt_dp, baf
                        dp,baf
                    ))

    except Exception as e:
        print(f"Error processing VCF file: {e}")
        return False

    # 写入结果文件
    if results:
        with open(output_file, 'w') as fout:
#            fout.write("#chrom\tpos\tvar_type\tgenotype\tDP\tAD_ref\tAD_alt\tAD_alt/DP\tBAF\n")
            fout.write("#chrom\tpos\tvar_type\tgenotype\tDP\tBAF\n")
            for res in results:
                fout.write("\t".join(map(str, res)) + "\n")

        print(f"Processing completed: {output_file}")
        print(f"Extracted {len(results)} variants")
        return True
    else:
        print("No variants found matching the criteria")
        return False

def main():
    # 创建自定义格式化类，隐藏大写参数名
    class CustomHelpFormatter(argparse.HelpFormatter):
        def _format_action_invocation(self, action):
            if not action.option_strings:
                metavar, = self._metavar_formatter(action, action.dest)(1)
                return metavar
            else:
                parts = []
                # 显示选项，但不显示大写的metavar
                if action.nargs == 0:
                    parts.extend(action.option_strings)
                else:
                    # 对于有参数的选项，只显示选项名
                    default = self._get_default_metavar_for_optional(action)
                    args_string = self._format_args(action, default)
                    for option_string in action.option_strings:
                        parts.append(option_string)
                    # 不添加args_string，这样就不会显示大写的metavar
                return ', '.join(parts)

    parser = argparse.ArgumentParser(
        description='Generate variant capture quality control metrics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )

    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')

    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--vcf', required=True,
                         metavar='<file.vcf>',
                         help='Input VCF files of the sample')
    required.add_argument('--out', required=True,
                         metavar='<outprefix>',
                         help='Prefix to name output files')

    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--regions',
                         metavar='<regions.bed>',
                         help='BED file for the target region (optional)')
    optional.add_argument('--pathogenic',
                         metavar='<pathogenic.tsv>',
                         default='/ocean/home/lanff/reference/variant/Clinvar_1_22.hg38.txt',
                         help='Pathogenic site files (default: Clinvar_1_22.hg38.txt)')

    args = parser.parse_args()

    # 构建输出文件路径
    output_file = f"{args.out}_variant_capture_QC.tsv"

    # 处理Variant Capture QC
    success = process_variant_capture_qc(
        args.vcf,
        args.regions,
        args.pathogenic,
        output_file
    )

    if success:
        print(f"Variant capture QC completed: {output_file}")
    else:
        print(f"Variant capture QC failed")

if __name__ == "__main__":
    main()
