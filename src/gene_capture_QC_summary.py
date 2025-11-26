#!/usr/bin/env python3
"""
gene_capture_QC_summary.py - Generate gene QC summary statistics
"""

import os
import csv
import argparse
from collections import defaultdict

def load_genes_from_file(gene_file):
    genes = set()
    with open(gene_file, 'r') as f:
        for line in f:
            genes.add(line.strip().split()[0])  # 假设每行第一个字段是基因名
    return genes

def safe_divide(numerator, denominator):
    return numerator / denominator if denominator != 0 else 0.0

def process_gene_qc_summary(gene_qc_file, output_file):
    coverage_stats = {
        "well_covered": defaultdict(int),
        "middle_covered": defaultdict(int),
        "poor_covered": defaultdict(int)
    }

    totals = {
        "gene": 0,
        "protein_coding": 0,
        "pseudogene": 0,
        "others": 0
    }

    sample_genes = set()
    try:
        with open(gene_qc_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row['gene_name']
                coverage = row['coverage']
                gene_type = row['gene_type']
                
                sample_genes.add(gene_name)
                
                if coverage not in coverage_stats:
                    continue

                coverage_stats[coverage]['gene'] += 1
                if gene_type == "protein_coding":
                    coverage_stats[coverage]['protein_coding'] += 1
                elif "pseudogene" in gene_type:
                    coverage_stats[coverage]['pseudogene'] += 1
                else:
                    coverage_stats[coverage]['others'] += 1
    except Exception as e:
        print(f"Error processing gene QC file: {e}")
        return

    # 计算总数
    for stat_type in ['gene', 'protein_coding', 'pseudogene', 'others']:
        totals[stat_type] = sum(
            coverage_stats[coverage_type][stat_type]
            for coverage_type in coverage_stats
        )

    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    # 生成每个覆盖类型的输出行
    with open(output_file, 'w') as f:
        # 写入表头
        header = [
            "coverage_quality",
            "all", "all_percent", 
            "protein_coding", "protein_coding_percent", 
            "pseudogene", "pseudogene_percent", 
            "others", "others_percent"
        ]
        f.write("\t".join(header) + "\n")
        
        for coverage_type in coverage_stats:
            stats = coverage_stats[coverage_type]

            # 计算百分比
            all_percent = safe_divide(stats['gene'], totals['gene']) * 100
            protein_coding_percent = safe_divide(stats['protein_coding'], totals['protein_coding']) * 100
            pseudogene_percent = safe_divide(stats['pseudogene'], totals['pseudogene']) * 100
            others_percent = safe_divide(stats['others'], totals['others']) * 100

            # 构建输出行
            output_row = [
                coverage_type,
                stats['gene'], f"{all_percent:.2f}",
                stats['protein_coding'], f"{protein_coding_percent:.2f}",
                stats['pseudogene'], f"{pseudogene_percent:.2f}",
                stats['others'], f"{others_percent:.2f}"
            ]

            # 写入结果
            f.write("\t".join(map(str, output_row)) + "\n")
        
        # 添加total行
        total_row = [
            "total",
            totals['gene'], "100.00",
            totals['protein_coding'], "100.00",
            totals['pseudogene'], "100.00",
            totals['others'], "100.00"
        ]
        f.write("\t".join(map(str, total_row)) + "\n")

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
        description='Generate gene QC summary statistics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--gene-qc', required=True, 
                         metavar='<gene_qc.tsv>',
                         help='Gene coverage')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.gene_qc):
        print(f"Error: Gene QC file {args.gene_qc} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}_gene_capture_QC.summary.tsv"
    
    process_gene_qc_summary(args.gene_qc, output_file)

if __name__ == '__main__':
    main()

