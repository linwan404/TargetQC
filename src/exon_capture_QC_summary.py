#!/usr/bin/env python3
"""
exon_capture_QC_summary.py - Generate exon quality control statistics
"""

import os
import pandas as pd
import argparse

def process_exon_qc(qc_exon_file, output_file):
    capture_tags = ['full', 'partial', 'near', 'no']  # capture_tag类型顺序
    valid_cover_tags = ['well', 'middle', 'poor']     
    
    results = []
    
    try:
        df = pd.read_csv(
            qc_exon_file,
            sep=r'\s+',
            comment='#',
            header=0,
            keep_default_na=False  # 关键修改点1：保留原始NA字符串
        )
        
        # 关键修改点1：只保留region=exon的行
        df = df[df['region'] == 'exon']
        # 关键修改点2：过滤无效cover_tag
        df = df[df['cover_tag'].isin(valid_cover_tags)]
        
        print(f"Processed {len(df)} rows from {qc_exon_file}")
        
    except Exception as e:
        print(f"Error processing {qc_exon_file}: {e}")
        return

    # 从输出文件路径提取样本名
#    sample_name = os.path.splitext(os.path.basename(output_file))[0].replace('_exon_capture_QC', '')
    
    # 统计各capture和type的数量
    grouped = df.groupby(['capture_tag', 'cover_tag'])['count'].sum().unstack(fill_value=0)
    print(f"Grouped data: {grouped}")
    print(f"Capture tags: {capture_tags}")
    
    has_na = 'NA' in grouped.index  # 判断是否存在NA标签

    if has_na:
        # 只处理NA标签的统计
        counts_na = grouped.loc['NA']
        full = counts_na.get('well', 0)
        other = counts_na.get('middle', 0)
        poor = counts_na.get('poor', 0)
        total = full + other + poor
        # 生成4个capture_tag行
        for capture in capture_tags:
            if capture == 'well':
                # all类型使用实际统计值
                results.append([
#                    sample_name,
                    capture,
                    full,
                    other,
                    poor,
                    total
                ])
            else:
                # all类型使用实际统计值
                results.append([
#                    sample_name,
                    capture,
                    0,
                    0,
                    0,
                    0
                ])

    else:
        # 无NA时，处理原四类capture标签
        for capture in capture_tags:
            if capture in grouped.index:
                counts = grouped.loc[capture]
            else:
                counts = pd.Series({t: 0 for t in valid_cover_tags})

            full = counts.get('well', 0)
            other = counts.get('middle', 0)
            poor = counts.get('poor', 0)
            total = full + other + poor

            results.append([
#                sample_name,
                capture,
                full,
                other,
                poor,
                total
            ])

    # 统计所有capture类型的总和
    total_counts = df.groupby('cover_tag')['count'].sum()
    full_total = total_counts.get('well', 0)
    other_total = total_counts.get('middle', 0)
    poor_total = total_counts.get('poor', 0)
    total_all = full_total + other_total + poor_total
    results.append([
#        sample_name,
        'coverage_total',  # 新增的capture类型标识
        full_total,
        other_total,
        poor_total,
        total_all
    ])
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # 写入结果文件（显式指定换行符和分隔符）
    with open(output_file, 'w', newline='\n') as f:  # 强制使用Linux换行符
        # 写入表头（用\t分隔）
        f.write("#capture\twell\tmiddle\tpoor\tcapture_total\n")
        # 写入数据行（逐行处理，确保使用\t分隔）
        for row in results:
            print(f"Row: {row}")
            line = "\t".join(map(str, row)) + "\n"  # \n换行，Linux格式
            f.write(line)
    
    print(f"Results written to {output_file}")

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
        description='Generate exon quality control statistics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--qc-exon', required=True, 
                         metavar='<exon.tsv>',
                         help='Exon coverage and capture input file of the sample')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.qc_exon):
        print(f"Error: QC exon file {args.qc_exon} does not exist.")
        return
    
    # 构建输出文件路径，添加_exon_capture_QC.tsv后缀
    output_file = f"{args.out}_exon_capture_QC.summary.tsv"
    
    # 处理样本
    process_exon_qc(args.qc_exon, output_file)

if __name__ == '__main__':
    main()
