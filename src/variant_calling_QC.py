#!/usr/bin/env python3
"""
variant_calling_QC.py - Generate variant calling quality control metrics
"""

import os
import pysam
import argparse

def calculate_baf(genotype, ad_alt, dp):
    """根据基因型和等位基因深度计算BAF"""
    if dp == 0:
        return "NA"
    
    allele_freq = ad_alt / dp
    
    if genotype == "hom":
        return abs(1 - allele_freq)
    elif genotype == "het":
        return abs(0.5 - allele_freq)
    else:
        return "NA"

def process_variant_qc(vcf_process, vcf_orig, variant_type, output_file):
    """处理VCF文件，提取QC指标"""
    
    # 检查文件是否存在
    if not all(os.path.exists(f) for f in [vcf_process, vcf_orig]):
        print(f"Error: One or more input files do not exist")
        return False
    
    # 处理文件1：提取符合条件的位点
    positions = {}
    
    try:
        with pysam.TabixFile(vcf_process) as vcf1:
            for row in vcf1.fetch():
                if row.startswith('#'):
                    continue
                fields = row.split('\t')
                # 检查过滤条件
                if (fields[6] == '.' and  # FILTER列
                    len(fields) > 10):
                    
                    # 根据variant_type选择不同的列
                    if variant_type == "FP":
                        fmt_data = fields[10].split(':')  # 第11列
                    else:  # FN或TP
                        fmt_data = fields[9].split(':')   # 第10列
                    
                    if len(fmt_data) > 5:
                        trio_type = fmt_data[1]
                        var_type = fmt_data[4]  # 变异类型
                        if var_type == "SNP":
                            var_type = "SNV"
                        genotype = fmt_data[5]  # 基因型
                        
                        # 如果genotype是homalt，则改为hom
                        if genotype == "homalt":
                            genotype = "hom"
                        if genotype == "hetalt":
                            genotype = "het"                        
                        if trio_type == variant_type:
                            chrom, pos = fields[0], fields[1]
                            positions[(chrom, pos)] = {
                                'var_type': var_type,
                                'genotype': genotype
                            }
    except Exception as e:
        print(f"Error processing {vcf_process}: {str(e)}")
        return False

    if not positions:
        print(f"No positions found for variant type: {variant_type}")
        return False

    # 处理文件2：提取深度信息
    results = []
    try:
        with pysam.TabixFile(vcf_orig) as vcf2:
            for (chrom, pos), info in positions.items():
                try:
                    pos_int = int(pos)
                    # 使用tabix索引快速定位
                    records = list(vcf2.fetch(chrom, pos_int-1, pos_int))
                except ValueError:
                    continue

                if not records:
                    continue

                # 取第一条记录（通常一个位置只有一条）
                record = records[0].split('\t')
                if len(record) < 10:
                    continue

                # 解析FORMAT和样本列
                fmt_keys = record[8].split(':')
                sample_info = record[9]
                fmt_vals = sample_info.split(':')
                fmt_dict = dict(zip(fmt_keys, fmt_vals))

                # 提取所需字段
                ad_str = fmt_dict.get('AD', 'NA')
                dp_str = fmt_dict.get('DP', 'NA')
                gq_str = fmt_dict.get('GQ', 'NA')
                var_type = info['var_type']
                genotype = info['genotype']

                if ad_str != 'NA' and ',' in ad_str:
                    ad_parts = ad_str.split(',')
                    ad_ref = int(ad_parts[0]) if ad_parts[0] else 0
                    ad_alt = int(ad_parts[1]) if len(ad_parts) > 1 and ad_parts[1] else 0
                    
                    try:
                        dp = int(dp_str) if dp_str and dp_str != 'NA' else 0
                        gq = int(gq_str) if gq_str and gq_str != 'NA' else 0
                        
                        # 计算AD_DP比值
                        ad_dp = ad_alt / dp if dp > 0 else 0.0
                        
                        # 计算BAF
                        baf = calculate_baf(genotype, ad_alt, dp)
                        
                        results.append((
                            chrom, pos, var_type, genotype,
                        #    ad_ref, ad_alt, dp, gq, ad_dp, baf
                           dp, baf
                        ))
                    except (ValueError, TypeError, ZeroDivisionError) as e:
                        continue
    except Exception as e:
        print(f"Error processing {vcf_orig}: {str(e)}")
        return False

    # 写入结果文件
    if results:
        # 确保输出目录存在
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"Created output directory: {output_dir}")
            
        with open(output_file, 'w') as fout:
#            fout.write("#chrom\tpos\tvar_type\tgenotype\tAD_ref\tAD_alt\tDP\tGQ\tAD_DP\tBAF\n")
            fout.write("#chrom\tpos\tvar_type\tgenotype\tDP\tBAF\n")
            for res in results:
                fout.write("\t".join(map(str, res)) + "\n")
        print(f"Processing completed: {output_file}")
        print(f"Extracted {len(results)} variants of type {variant_type}")
        return True
    else:
        print("No valid data to output.")
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
        description='Generate variant calling quality control metrics',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--vcf-orig', required=True, 
                         metavar='<file.vcf>',
                         help='Original VCF file')
    required.add_argument('--vcf-process', required=True, 
                         metavar='<file.vcf>',
                         help='Processed VCF file')
    required.add_argument('--type', required=True,
                         choices=['FN', 'FP', 'TP'],
                         metavar='<type>',
                         help='Variant site type (FN, FP, TP)')
    required.add_argument('--out', required=True, 
                         metavar='<outprefix>',
                         help='Prefix to name output files')
    
    args = parser.parse_args()
    
    # 构建输出文件路径
    output_file = f"{args.out}_{args.type}_QC.tsv"
    
    # 处理Variant QC
    success = process_variant_qc(
        args.vcf_process,
        args.vcf_orig,
        args.type,
        output_file
    )
    
    if success:
        print(f"Variant QC completed: {output_file}")
    else:
        print(f"Variant QC failed")

if __name__ == "__main__":
    main()
