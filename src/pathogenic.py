#!/usr/bin/env python3
"""
pathogenic.py - Annotate false negative sites with pathogenic information
"""

import os
import argparse
from collections import defaultdict

def load_reference(file_path):
    """加载参考文件到集合，键为(chrom, start+1)"""
    ref_set = set()
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue  # 忽略不完整的行
            chrom = parts[0]
            start = int(parts[1])
            pos = start + 1
            ref_set.add((chrom, pos))
    return ref_set

def process_single_sample(fn_file, output_prefix, pathogenic_file):
    """处理单个样本的假阴性位点"""
    
    # 检查输入文件是否存在
    if not os.path.exists(fn_file):
        print(f"Error: FN file {fn_file} does not exist.")
        return False
    
    if not os.path.exists(pathogenic_file):
        print(f"Error: Pathogenic file {pathogenic_file} does not exist.")
        return False
    
    # 加载致病位点文件
    pathogenic_set = load_reference(pathogenic_file)
    
    # 构建输出文件路径
    output_file = f"{output_prefix}_pathogenic.tsv"
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    counts = defaultdict(int)
    
    try:
        with open(fn_file, 'r') as fin, \
             open(output_file, 'w') as fout:
            # 写入表头
            header = "#chrom\tpos\tclinical\n"
            fout.write(header)
            
            for line in fin:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                chrom = parts[0]
                pos = int(parts[1])
                
                # 检查是否在致病位点集合中
                in_pathogenic = (chrom, pos) in pathogenic_set
                
                clinical = 'Pathogenic' if in_pathogenic else 'None'
                counts[clinical] += 1
                
                # 仅输出致病位点
                if clinical == 'Pathogenic':
                    output_line = f"{chrom}\t{pos}\t{clinical}\n"
                    fout.write(output_line)
        
        # 输出统计信息
        print(f"Processing completed: {output_file}")
        print(f"Pathogenic sites: {counts['Pathogenic']}")
        print(f"Non-pathogenic sites: {counts['None']}")
        print(f"Total sites processed: {sum(counts.values())}")
        
        return True
        
    except Exception as e:
        print(f"Error processing files: {e}")
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
        description='Annotate false negative sites with pathogenic information',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--fn', required=True, 
                         metavar='<FN.tsv>',
                         help='File of false negative sites')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--pathogenic', 
                         metavar='<pathogenic.tsv>',
                         default='./clinvar.txt',
                         help='Pathogenic site files (default: clinvar.txt)')
    
    args = parser.parse_args()
    
    # 处理单个样本
    success = process_single_sample(
        args.fn,
        args.out,
        args.pathogenic
    )
    
    if success:
        print(f"Pathogenic annotation completed: {args.out}_pathogenic.tsv")
    else:
        print(f"Pathogenic annotation failed")

if __name__ == "__main__":
    main()
