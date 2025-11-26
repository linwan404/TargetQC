#!/usr/bin/env python3
"""
base.py - Process mosdepth output files for coverage analysis
"""

import os
import pandas as pd
import argparse
import numpy as np
import gzip

def process_sample(thresholds_file, regions_file, output_file, allowed_chroms):
    try:
        with gzip.open(thresholds_file, 'rt') as f:
            header_line = f.readline().strip()
            columns = header_line[1:].strip().split('\t')
            fixed_columns = columns[:4]
            threshold_columns = columns[4:]
            
            print(f"Detected threshold columns: {threshold_columns}")        
        df_a = pd.read_csv(
            thresholds_file,
            sep='\t',
            comment='#',
            compression='gzip',
            names=fixed_columns + threshold_columns
        )
        
        if allowed_chroms is not None:
            df_a = df_a[df_a['chrom'].isin(allowed_chroms)]
        
        total_bases = (df_a['end'] - df_a['start']).sum()
        
        threshold_results = {}
        for col in threshold_columns:
            total_count = df_a[col].sum()
            percent = (total_count / total_bases * 100) if total_bases != 0 else 0.0
            threshold_value = col.replace('X', 'X')
            threshold_results[threshold_value] = {
                'count': total_count,
                'percent': percent
            }
            
    except Exception as e:
        print(f"Error processing thresholds file {thresholds_file}: {e}")
        return

    try:
        df_b = pd.read_csv(
            regions_file,
            sep='\t',
            header=None,
            compression='gzip',
            names=['chrom', 'start', 'end', 'avg_cov']
        )
        if allowed_chroms is not None:
            df_b = df_b[df_b['chrom'].isin(allowed_chroms)]

        df_b['length'] = df_b['end'] - df_b['start']
        
        weighted_sum = (df_b['avg_cov'] * df_b['length']).sum()
        avg_coverage = weighted_sum / total_bases if total_bases != 0 else 0.0
        
    except Exception as e:
        print(f"Error processing regions file {regions_file}: {e}")
        return
    
    with open(output_file, 'w') as out_f:
        out_f.write("#thresholds\tsize\tavg\tcount\tpercent\n")
        for threshold, results in threshold_results.items():
            line = f"{threshold}\t{total_bases}\t{avg_coverage:.2f}\t{results['count']}\t{results['percent']:.2f}\n"
            out_f.write(line)
    
    print(f"Results written to {output_file}")

def main():
    class CustomHelpFormatter(argparse.HelpFormatter):
        def _format_action_invocation(self, action):
            if not action.option_strings:
                metavar, = self._metavar_formatter(action, action.dest)(1)
                return metavar
            else:
                parts = []
                if action.nargs == 0:
                    parts.extend(action.option_strings)
                else:
                    default = self._get_default_metavar_for_optional(action)
                    args_string = self._format_args(action, default)
                    for option_string in action.option_strings:
                        parts.append(option_string)
                return ', '.join(parts)

    parser = argparse.ArgumentParser(
        description='Process mosdepth output files for coverage analysis',
        usage='base.py --thresholds-base <thresholds.bed.gz> --regions-base <regions.bed.gz> --out <outprefix>',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )

    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    required = parser.add_argument_group('Required options')
    required.add_argument('--thresholds-base', required=True, 
                         metavar='<thresholds.bed.gz>',
                         help='Thresholds files of the sample')
    required.add_argument('--regions-base', required=True,
                         metavar='<regions.bed.gz>',
                         help='Target regions files of the sample')
    required.add_argument('--out', required=True,
                         metavar='<outprefix>', 
                         help='Prefix to name output files')
    
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--chrom', nargs='+', default=None,
                         metavar='[chr1,chr2,...]',
                         help='List of target chromosomes (eg.chr1 chr2 ...)(default: chr1-22)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.thresholds_base):
        print(f"Error: Thresholds file {args.thresholds_base} does not exist.")
        return
    
    if not os.path.exists(args.regions_base):
        print(f"Error: Regions file {args.regions_base} does not exist.")
        return
    
    allowed_chroms = args.chrom

    if allowed_chroms is not None:
        print(f"Processing with chromosomes: {allowed_chroms}")
    else:
        print("Processing all chromosomes (no chromosome filter applied)")    

    output_file = f"{args.out}_base.tsv"    
    process_sample(args.thresholds_base, args.regions_base, output_file, allowed_chroms)

if __name__ == '__main__':
    main()
