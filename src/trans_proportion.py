#!/usr/bin/env python3
"""
trans_proportion.py - Calculate gene coverage proportion statistics
"""

import os
import argparse

def process_gene_coverage(gene_file, output_file):
    # 定义阈值范围
    thresholds = list(range(0, 101))  # 生成0到100的整数列表
    
    # 读取输入文件的第一行获取列名
    with open(gene_file, 'r') as f:
        header = f.readline().strip().split('\t')
    
    # 获取阈值列名（从第5列开始）
    threshold_columns = header[4:]
    
    # 从输出文件路径提取样本名
#    sample_name = os.path.splitext(os.path.basename(output_file))[0].replace('_trans_proportion', '')
    
    # 初始化计数器
    count_dict = {}
    for col in threshold_columns:
        count_dict[col] = [0] * len(thresholds)

    # 处理基因文件
    with open(gene_file, 'r') as f:
        next(f)  # 跳过标题行
        for line_num, line in enumerate(f, 2):
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < len(header):
                print(f"{gene_file} 第{line_num}行列数不足")
                continue

            try:
                # 处理每个阈值列
                for i, col_name in enumerate(threshold_columns):
                    value = float(cols[4+i])
                    # 更新计数器
                    for j, thres in enumerate(thresholds):
                        if value >= thres:
                            count_dict[col_name][j] += 1
            except ValueError:
                print(f"{gene_file} 第{line_num}行数值格式错误")
                continue

    # 写入结果
    with open(output_file, 'w', newline='\n') as f_out:
        # 写入标题行
        header_out = ['tag'] + [f"{thres}%" for thres in thresholds]
        f_out.write('\t'.join(header_out) + '\n')
        
        # 写入每个标签的结果
        for tag, counts in count_dict.items():
            row = [
                tag,
                *[str(c) for c in counts]
            ]
            f_out.write('\t'.join(row) + '\n')
    
    print(f"Processing completed: {output_file}")

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
        description='Calculate gene coverage proportion statistics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--gene', required=True, 
                         metavar='<gene.tsv>',
                         help='Gene coverage input file of the sample')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.gene):
        print(f"Error: Gene file {args.gene} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}_proportion.tsv"
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # 处理基因覆盖度文件
    process_gene_coverage(args.gene, output_file)

if __name__ == '__main__':
    main()
