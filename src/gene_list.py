#!/usr/bin/env python3
"""
gene_list.py - Generate gene lists based on coverage and gene type
"""

import os
import csv
import argparse
from collections import defaultdict

def read_gene_sets(gene_file):
    """读取基因文件，返回基因集合"""
    genes = set()
    with open(gene_file, 'r') as f:
        for line in f:
            genes.add(line.strip())
    return genes

def process_gene_list(gene_qc_file, output_prefix, omim_genes, gene_type):
    """处理基因列表"""
    # 初始化数据结构：每个覆盖类型维护基因集合
    coverage_data = {
        'well_covered': set(),
        'poor_covered': set()
    }

    try:
        with open(gene_qc_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row['gene_name']
                gene_type_row = row['gene_type']
                coverage = row['coverage']
                
                # 只处理目标覆盖类型
                if coverage not in coverage_data:
                    continue
                
                # 根据基因类型筛选
                if gene_type == 'all_gene':
                    # 所有基因都包含
                    coverage_data[coverage].add(gene_name)
                elif gene_type == 'protein_coding_gene' and gene_type_row == 'protein_coding':
                    coverage_data[coverage].add(gene_name)
                elif gene_type == 'pseudogene' and 'pseudogene' in gene_type_row:
                    coverage_data[coverage].add(gene_name)
                elif gene_type == 'others' and gene_type_row not in ['protein_coding'] and 'pseudogene' not in gene_type_row:
                    coverage_data[coverage].add(gene_name)
                elif gene_type == 'omim_gene' and gene_name in omim_genes and gene_type_row == 'protein_coding':
                    coverage_data[coverage].add(gene_name)
    
    except Exception as e:
        print(f"Error processing gene QC file: {e}")
        return
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # 为每个覆盖类型创建文件
    for coverage_type in ['well_covered', 'poor_covered']:
        output_file = f"{output_prefix}_gene_{coverage_type}.tsv"
        
        with open(output_file, 'w') as f:
            # 按字母顺序写入基因名称
            for gene in sorted(coverage_data[coverage_type]):
                f.write(f"{gene}\n")
        
        print(f"Created {output_file} with {len(coverage_data[coverage_type])} genes")
    
    return coverage_data

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
        description='Generate gene lists based on coverage and gene type',
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
    required.add_argument('--gene-type', required=True,
                         choices=['all_gene', 'protein_coding_gene', 'pseudogene', 'others', 'omim_gene'],
                         metavar='<gene-type>',
                         help='Gene type (all_gene, protein_coding_gene, pseudogene, others, omim_gene)')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--omim', 
                         metavar='<omim_gene>',
                         default='./omim.txt',
                         help='Protein-coding genes in OMIM database have clear inheritance patterns (default: omim.txt)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.gene_qc):
        print(f"Error: Gene QC file {args.gene_qc} does not exist.")
        return
    
    # 如果基因类型是omim_gene，检查OMIM文件是否存在
    if args.gene_type == 'omim_gene' and not os.path.exists(args.omim):
        print(f"Error: OMIM gene file {args.omim} does not exist.")
        return
    
    # 加载OMIM基因列表（仅在需要时）
    omim_genes = set()
    if args.gene_type == 'omim_gene':
        omim_genes = read_gene_sets(args.omim)
        print(f"Loaded {len(omim_genes)} OMIM genes from {args.omim}")
    
    # 处理基因列表
    coverage_data = process_gene_list(
        args.gene_qc,
        args.out,
        omim_genes,
        args.gene_type
    )
    
    # 输出统计信息
    if coverage_data:
        well_count = len(coverage_data['well_covered'])
        poor_count = len(coverage_data['poor_covered'])
        print(f"Processing completed:")
        print(f"  Well covered genes: {well_count}")
        print(f"  Poor covered genes: {poor_count}")
        print(f"  Output files: {args.out}_gene_well_covered.tsv, {args.out}_gene_poor_covered.tsv")

if __name__ == '__main__':
    main()
