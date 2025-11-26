#!/usr/bin/env python3
"""
gene_capture_QC.py - Generate gene quality control statistics
"""

import os
import csv
import argparse
from collections import defaultdict

def process_gene_qc(trans_exon_file, trans_cds_file, output_file, well_threshold, poor_threshold):
    # 处理CDS文件，构建gene_name到max_full和max_poor的映射
    cds_data = defaultdict(lambda: {'max_full': 0.0, 'max_poor': 0.0})
    
    try:
        with open(trans_cds_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # 跳过表头
            for row in reader:
                if len(row) < 7:
                    continue
                gene_name = row[1]
                try:
                    full = float(row[4])
                    poor = float(row[6])
                except ValueError:
                    continue
                if full > cds_data[gene_name]['max_full']:
                    cds_data[gene_name]['max_full'] = full
                if poor > cds_data[gene_name]['max_poor']:
                    cds_data[gene_name]['max_poor'] = poor
    except Exception as e:
        print(f"Error processing CDS file: {e}")
        return

    # 处理转录本文件并写入输出
    well, poor_cov, other = 0, 0, 0

    try:
        with open(trans_exon_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
            reader = csv.reader(fin, delimiter='\t')
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow(['#gene_id', 'gene_name', 'gene_type', 'transcript_id', 'coverage'])
            next(reader)  # 跳过表头

            for row in reader:
                if len(row) < 7:
                    continue
                gene_id = row[0]
                gene_name = row[1]
                gene_type = row[2]
                transcript_id = row[3]

                # 判断覆盖类型
                if gene_type == 'protein_coding':
                    data = cds_data.get(gene_name, {'max_full': 0.0, 'max_poor': 0.0})
                    max_full = data['max_full']
                    max_poor = data['max_poor']
                else:
                    try:
                        max_full = float(row[4])
                        max_poor = float(row[6])
                    except ValueError:
                        max_full, max_poor = 0.0, 0.0

                # 确定覆盖类型
                if max_full >= well_threshold:
                    coverage = 'well_covered'
                    well += 1
                elif max_poor >= poor_threshold:
                    coverage = 'poor_covered'
                    poor_cov += 1
                else:
                    coverage = 'middle_covered'
                    other += 1

                writer.writerow([gene_id, gene_name, gene_type, transcript_id, coverage])
        
        print(f"Gene QC results: well_covered={well}, poor_covered={poor_cov}, middle_covered={other}")
        
    except Exception as e:
        print(f"Error processing exon file: {e}")
        return

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
        description='Generate gene quality control statistics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--trans-exon', required=True, 
                         metavar='<trans_exon.tsv>',
                         help='Exon region coverage in each gene')
    required.add_argument('--trans-cds', required=True, 
                         metavar='<trans_cds.tsv>',
                         help='CDS region coverage in each gene')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--well', type=float, default=95.0,
                         metavar='<float>',
                         help='Proportion of bases with well coverage in this region (default: 95.0)')
    optional.add_argument('--poor', type=float, default=5.0,
                         metavar='<float>',
                         help='Proportion of bases with poor coverage in this region (default: 5.0)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.trans_exon):
        print(f"Error: Trans-exon file {args.trans_exon} does not exist.")
        return
    
    if not os.path.exists(args.trans_cds):
        print(f"Error: Trans-cds file {args.trans_cds} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}_gene_capture_QC.tsv"
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    print(f"Using well threshold: {args.well}")
    print(f"Using poor threshold: {args.poor}")
    
    # 处理基因QC
    process_gene_qc(
        args.trans_exon,
        args.trans_cds,
        output_file,
        args.well,
        args.poor
    )
    
    print(f"Processing completed: {output_file}")

if __name__ == '__main__':
    main()
