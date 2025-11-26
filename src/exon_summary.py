#!/usr/bin/env python3
import argparse
import pandas as pd
from collections import defaultdict
import sys
import os

def process_exon_qc(input_file, output_file):
    """
    处理exon QC文件，统计不同type和capture的组合，生成矩阵文件
    """
    try:

        df = pd.read_csv(input_file, sep='\t')
        required_columns = ['type', 'capture']
        for col in required_columns:
            if col not in df.columns:
                sys.exit(1)
        capture_types = ['full', 'partial', 'near', 'no']
        type_categories = ['well', 'middle', 'poor']
        matrix = defaultdict(lambda: defaultdict(int))
        for _, row in df.iterrows():
            capture = row['capture']
            type_val = row['type']
            
            capture_base = capture.replace('_capture', '') if '_capture' in capture else capture
            
            type_base = type_val.split('_')[0] if '_' in type_val else type_val
            
            if capture_base in capture_types and type_base in type_categories:
                matrix[capture_base][type_base] += 1
        
        row_totals = {}
        col_totals = defaultdict(int)
        
        for capture in capture_types:
            row_total = 0
            for type_cat in type_categories:
                row_total += matrix[capture][type_cat]
                col_totals[type_cat] += matrix[capture][type_cat]
            row_totals[capture] = row_total
        
        grand_total = sum(row_totals.values())
        
        output_lines = []
        
        header = "#capture\twell\tmiddle\tpoor\tcapture_total"
        output_lines.append(header)
        
        for capture in capture_types:
            line = f"{capture}"
            for type_cat in type_categories:
                line += f"\t{matrix[capture][type_cat]}"
            line += f"\t{row_totals[capture]}"
            output_lines.append(line)
        
        total_line = "coverage_total"
        for type_cat in type_categories:
            total_line += f"\t{col_totals[type_cat]}"
        total_line += f"\t{grand_total}"
        output_lines.append(total_line)
        
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_lines))
        
        print('\n'.join(output_lines))
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='统计exon QC文件中不同type和capture的组合')
    parser.add_argument('--qc-exon', required=True, help='输入exon QC文件路径')
    parser.add_argument('--out', required=True, help='输出文件路径，包含后缀')
    
    args = parser.parse_args()
    
    output_file = f"{args.out}_exon_capture_QC.summary.tsv"
    
    process_exon_qc(args.qc_exon, output_file)

if __name__ == '__main__':
    main()
