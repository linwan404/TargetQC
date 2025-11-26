#!/usr/bin/env python3
"""
trans.py - Process transcript coverage data
"""

import gzip
import os
import csv
import argparse
from collections import defaultdict

def process_files(file_a_path, file_b_path, file_c_path, output_path, up_threshold, down_threshold):
    # 读取文件a，建立trans_id到gene_name的映射
    trans_info = {}
    with open(file_a_path, 'r') as f:
        for line in f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 1:
                continue
            trans_id = parts[0]
            gene_name = parts[1]
            trans_info[trans_id] = gene_name

    # 读取文件b，建立transcript信息字典
    transcripts = defaultdict(dict)
    with open(file_b_path, 'r') as f:
        current_trans_id = None
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            feature_type = parts[1]
            if feature_type == 'transcript':
                current_trans_id = parts[5]
                if current_trans_id not in trans_info:
                    continue
                transcripts[current_trans_id] = {
                    'gene_id': parts[4],
                    'gene_type': parts[6],
                    'exons': [],
                    'gene_name': trans_info[current_trans_id]
                }
            elif feature_type == 'exon' and current_trans_id and current_trans_id in transcripts:
                chrom = parts[0]
                start = int(parts[2]) - 1  # 转换为0-based
                end = int(parts[3])
                transcripts[current_trans_id]['exons'].append((chrom, start, end))

    # 读取文件c，建立位置索引
    c_index = defaultdict(list)
    
    # 动态检测阈值列名
    with gzip.open(file_c_path, 'rt') as f:
        # 查找包含列名的行（以#开头的行）
        for line in f:
            if line.startswith('#'):
                header_line = line.strip()
                break
        
        if header_line is None:
            raise ValueError("No header line found in coverage file")
        print(header_line)
        # 去掉开头的#并分割列名
        columns = header_line[1:].strip().split('\t')
        
        # 前4列固定为 chrom, start, end, region
        fixed_columns = columns[:4]
        # 剩余的列为阈值列
        threshold_columns = columns[4:]
        
        # 构建阈值列名
        up_column = f"{up_threshold}X"
        down_column = f"{down_threshold}X"
        
        # 检查阈值列是否存在
        if up_column not in threshold_columns:
            raise ValueError(f"Threshold column {up_column} not found in coverage file")
        if down_column not in threshold_columns:
            raise ValueError(f"Threshold column {down_column} not found in coverage file")
        
        # 获取列索引
        up_idx = threshold_columns.index(up_column) + 4  # +4 因为前4列是固定的
        down_idx = threshold_columns.index(down_column) + 4
        # 重置文件指针到数据开始处
        f.seek(0)

        # 处理数据行
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < max(up_idx, down_idx) + 1:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            x_up = int(parts[up_idx])
            x_down = int(parts[down_idx])
            key = (chrom, start, end)
            # 直接覆盖，每个位置只保留最后出现的覆盖度数据
            c_index[key] = (x_down, x_up)
#            print(c_index)
    # 计算结果
    results = []
    for trans_id in trans_info:
        if trans_id not in transcripts:
            continue

        t_data = transcripts[trans_id]
        total_poor = total_other = total_full = total_bases = 0.0
        for (chrom, start, end) in t_data['exons']:
            key = (chrom, start, end)
            exon_length = end - start
            total_bases += exon_length
            
            if key not in c_index:
#                print("not key/",key)
                total_poor += exon_length
                continue
            
            # 只处理一次，无论重复多少次
            x_up, x_down = c_index[key]
            poor = exon_length - x_down
            other = x_up - x_down
            full = x_up
            total_poor += poor
            total_other += other
            total_full += full

        # 计算百分比
        if total_bases == 0:
            full_pct = other_pct = poor_pct = 0.0
        else:
            full_pct = (total_full / total_bases) * 100
            other_pct = (total_other / total_bases) * 100
            poor_pct = (total_poor / total_bases) * 100
        
        results.append([
            t_data['gene_id'],
            t_data['gene_name'],
            t_data['gene_type'],
            trans_id,
            f"{full_pct:.2f}",
            f"{other_pct:.2f}",
            f"{poor_pct:.2f}"
        ])

    # 写入输出文件
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['#gene_id', 'gene_name', 'gene_type', 'transcript_id', f'≥{up_threshold}X', f'{up_threshold}_{down_threshold}X', f'≤{down_threshold}X'])
        writer.writerows(results)

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
        description='Process transcript coverage data',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--coverage', required=True, 
                         metavar='<coverage.tsv>',
                         help='Exon coverage input file of the sample')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--gene', 
                         metavar='<transcript.tsv>',
                         default='./target_transcript.tsv',
                         help='Specifies the gene transcript (default: target_transcript.tsv)')
    optional.add_argument('--reference', 
                         metavar='<gene_reference.tsv>',
                         default='./genes.gtf',
                         help='Gene annotation (default: genes.gtf)')
    optional.add_argument('--up-gene', type=int, default=30,
                         metavar='<int>',
                         help='Exon coverage reached the standard threshold (default: 30)')
    optional.add_argument('--down-gene', type=int, default=20,
                         metavar='<int>',
                         help='Exon coverage falls below the standard threshold (default: 20)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.coverage):
        print(f"Error: Coverage file {args.coverage} does not exist.")
        return
    
    if not os.path.exists(args.gene):
        print(f"Error: Transcript file {args.gene} does not exist.")
        return
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference file {args.reference} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}_trans.tsv"
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    print(f"Processing with up threshold: {args.up_gene}X")
    print(f"Processing with down threshold: {args.down_gene}X")
    
    # 处理样本
    process_files(
        args.gene,
        args.reference,
        args.coverage,
        output_file,
        args.up_gene,
        args.down_gene
    )
    
    print(f"Processing completed: {output_file}")

if __name__ == '__main__':
    main()
