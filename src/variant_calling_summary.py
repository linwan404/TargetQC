#!/usr/bin/env python2
"""
variant_calling_summary.py - Generate variant calling summary report
"""

import os
import csv
import argparse

def process_variant_summary(input_csv, output_file):
    """处理单个样本的变异检测总结文件"""
    results = []
    
    try:
        # 读取CSV文件内容
        with open(input_csv, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row.get('Filter') == 'PASS':
                    variant_type = row['Type']
                    if variant_type == 'SNP':
                        variant_type = 'SNV'
                    print(variant_type)
                    record = {
                        '#type': variant_type,
                        'TP': row['TRUTH.TP'],
                        'FN': row['TRUTH.FN'],
                        'FP': row['QUERY.FP'],
                        'Recall': row['METRIC.Recall'],
                        'Precision': row['METRIC.Precision'],
                        'F1_Score': row['METRIC.F1_Score']
                    }
                    results.append(record)
    
    except Exception as e:
        print(f"Error processing CSV file {input_csv}: {e}")
        return False
    
    # 写入输出文件
    if results:
        # 定义输出文件的列顺序
        fieldnames = [
            '#type',
            'TP', 'FN', 'FP',
            'Recall', 'Precision', 'F1_Score'
        ]

        # 确保输出目录存在
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")

        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
        print(f"Summary report generated: {output_file}")
        print(f"Processed {len(results)} records with Filter=PASS")
        return True
    else:
        print("No valid data with Filter=PASS to output.")
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
        description='Generate variant calling summary report',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--input', required=True, 
                         metavar='<file.csv>',
                         help='Detection accuracy file')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}.summary.tsv"
    
    # 处理变异检测总结
    success = process_variant_summary(args.input, output_file)
    
    if success:
        print(f"Processing completed: {output_file}")
    else:
        print(f"Processing failed: {args.input}")

if __name__ == '__main__':
    main()
