#!/usr/bin/env python3
"""
exon_coverage.py - Analyze exon coverage based on mosdepth output files
"""

import os
import pandas as pd
import argparse

def process_exon_coverage(thresholds_file, regions_file, output_file, valid_chromosomes, up_threshold, down_threshold, base_threshold):
    try:
        # 读取阈值文件
        df_thresholds = pd.read_csv(thresholds_file, sep='\t', header=0, compression='gzip')
        
        # 重命名列，将#chrom改为chrom
        df_thresholds = df_thresholds.rename(columns={'#chrom': 'chrom'})
        
        if valid_chromosomes is not None:
            df_thresholds = df_thresholds[df_thresholds['chrom'].isin(valid_chromosomes)]
        
        # 读取区域文件
        df_exons = pd.read_csv(regions_file, sep='\t', header=None, 
                              names=['chrom', 'start', 'end', 'coverage'], 
                              compression='gzip')
        
        if valid_chromosomes is not None:
            df_exons = df_exons[df_exons['chrom'].isin(valid_chromosomes)]
        
        print(f"Merging dataframes...")
        # 合并数据框
        merged_df = pd.merge(df_thresholds, df_exons, on=['chrom', 'start', 'end'])
        print(f"Merged DataFrame has {len(merged_df)} rows")
        
        # 构建阈值列名
        up_column = f"{up_threshold}X"
        down_column = f"{down_threshold}X"
        base_column = base_threshold
        # 检查阈值列是否存在
        if up_column not in merged_df.columns:
            raise ValueError(f"Threshold column {up_column} not found in thresholds file")
        if down_column not in merged_df.columns:
            raise ValueError(f"Threshold column {down_column} not found in thresholds file")
        
        # 筛选不同类型的区域
        full_exon = merged_df[(merged_df['end'] - merged_df['start']) == merged_df[up_column]]
        poor_exon = merged_df[(merged_df['end'] - merged_df['start'] - merged_df[down_column]) >= base_column]
        other_exon = merged_df[~(((merged_df['end'] - merged_df['start']) == merged_df[up_column]) | 
                                ((merged_df['end'] - merged_df['start'] - merged_df[down_column]) >= base_column))]
        
        print(f"Full exon regions: {len(full_exon)}")
        print(f"Poor exon regions: {len(poor_exon)}")
        print(f"Other exon regions: {len(other_exon)}")
        
        # 准备输出数据
        full_exon_data = full_exon[['chrom', 'start', 'end', 'coverage']].copy()
        full_exon_data['Type'] = f'well_covered'
        
        poor_exon_data = poor_exon[['chrom', 'start', 'end', 'coverage']].copy()
        poor_exon_data['Type'] = f'poor_covered'
        
        other_exon_data = other_exon[['chrom', 'start', 'end', 'coverage']].copy()
        other_exon_data['Type'] = f'middle_covered'
        
        # 合并结果
        results = pd.concat([full_exon_data, poor_exon_data, other_exon_data], ignore_index=True)
        
        # 确保输出目录存在
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        
        # 保存结果
        results.to_csv(output_file, sep='\t', index=False)
        print(f"Results saved to {output_file}")
        
    except Exception as e:
        print(f"Error processing files: {e}")
        raise

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
        description='Analyze exon coverage based on mosdepth output files',
        usage='exon_coverage.py --thresholds-exon <thresholds.bed.gz> --regions-exon <regions.bed.gz> --out <outprefix>',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--thresholds-exon', required=True, 
                         metavar='<thresholds.bed.gz>',
                         help='Thresholds files of the sample')
    required.add_argument('--regions-exon', required=True, 
                         metavar='<regions.bed.gz>',
                         help='Target regions files of the sample')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--up-exon', type=int, default=30,
                         metavar='<int>',
                         help='Exon coverage reached the standard threshold (default: 30)')
    optional.add_argument('--down-exon', type=int, default=20,
                         metavar='<int>',
                         help='Exon coverage falls below the standard threshold (default: 20)')
    optional.add_argument('--base-exon', type=int, default=5,
                         metavar='<int>',
                         help='The number of bases in an exon below the threshold for not meeting the standard (default: 5)')
    optional.add_argument('--chrom', nargs='+', default=None,
                         metavar='[chr1,chr2,...]',
                         help='List of target chromosomes (eg.chr1 chr2 ...)')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.thresholds_exon):
        print(f"Error: Thresholds file {args.thresholds_exon} does not exist.")
        return
    
    if not os.path.exists(args.regions_exon):
        print(f"Error: Regions file {args.regions_exon} does not exist.")
        return
    
    valid_chromosomes = args.chrom

    if valid_chromosomes is not None:
        print(f"Processing with chromosomes: {valid_chromosomes}")
    else:
        print("Processing all chromosomes (no chromosome filter applied)")
    
    print(f"Processing with chromosomes: {valid_chromosomes}")
    print(f"Using up threshold: {args.up_exon}X")
    print(f"Using down threshold: {args.down_exon}X")
    print(f"Using base threshold: {args.base_exon}bp")
    output_file = f"{args.out}_exon_coverage.tsv"    
    # 处理样本
    process_exon_coverage(
        args.thresholds_exon, 
        args.regions_exon, 
        output_file, 
        valid_chromosomes, 
        args.up_exon, 
        args.down_exon,
        args.base_exon
    )

if __name__ == '__main__':
    main()
