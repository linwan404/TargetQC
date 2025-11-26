#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exon_gene.py - Process gene and exon coverage data
"""

import os
import csv
import argparse
from collections import defaultdict

def read_targets(file_a_path):
    """返回：
    - targets: [(trans_id, gene_name), ...]
    - gene_order: {gene_name: 出现顺序索引}
    """
    gene_order = {}
    targets = []
    order_idx = 0
    with open(file_a_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            trans_id = parts[0]
            gene_name = parts[1]
            targets.append((parts[0], parts[1]))
            # 记录首次出现的基因顺序
            if gene_name not in gene_order:
                gene_order[gene_name] = order_idx
                order_idx += 1
    return targets, gene_order

def parse_file_d(file_d_path, targets):
    # 初始化数据结构
    precise_index = defaultdict(list)  # {(chrom,start,end): [trans_info]}
    transcript_regions = defaultdict(lambda: {'exon': []})  # {tid: regions}
    target_tids = {tid for tid, _ in targets}

    with open(file_d_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 9: continue

            chrom = parts[0]
            feature_type = parts[1]
            start = int(parts[2]) - 1  # 转换为0-based
            end = int(parts[3])
            tid = parts[5]
            if tid not in target_tids:
                continue

            # 提取基因信息（关键修复）
            if feature_type == "transcript":
                transcript_regions[tid]['gene_id'] = parts[4]          # 第5列是gene_id
                transcript_regions[tid]['gene_type'] = parts[6]        # 第7列是gene_type
                transcript_regions[tid]['gene_name'] = parts[7]        # 第8列是gene_name

            # 处理exon
            elif feature_type == 'exon':
                exon_num = int(parts[8])
                region = {
                    'chrom': chrom, 'start': start, 'end': end,
                    'exon_number': exon_num, 'type': 'exon'
                }
                # 更新三级索引
                precise_key = (chrom, start, end)
                precise_index[precise_key].append(('exon', tid, region))
                transcript_regions[tid]['exon'].append(region)
    return precise_index, transcript_regions

def find_matching_region(sample_row, precise_idx):
    chrom = sample_row['#chrom']
    s = int(sample_row['start'])
    e = int(sample_row['end'])
    # 存储未匹配信息用于最终输出
    mismatch_info = {
        'chrom': chrom,
        'start': s,
        'end': e,
        'match_stage': None
    }
    # 第一级：精准匹配
    precise_matches = precise_idx.get((chrom, s, e), [])
    # 优先匹配当前区域类型
    type_matches = [m for m in precise_matches]
    if type_matches:
        return type_matches
    # 构建详细未匹配信息
    mismatch_details = (
        f"NO_MATCH: | "
        f"location={chrom}:{s}-{e} | "
        f"attempted_stages={['precise', 'fuzzy', 'position_fuzzy']}"
    )
    # 打印未匹配信息（实际使用时建议写入日志文件）
#    print(mismatch_details)
    return []

def process_sample(sample_file, region_type, precise_idx, trans_regions):
    results = defaultdict(lambda: {
        'total_len': 0, 'total_cov': 0, 'total_cap': 0,
        'tags': defaultdict(lambda: {'count': 0, 'exon': set()})
    })

    with open(sample_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            matches = find_matching_region(row, precise_idx)
            if not matches:
                continue
            # 使用文件D中的实际坐标
            for match in matches:
                rtype, tid, region = match
                region_len = region['end'] - region['start']
                exon_number = region['exon_number']
                coverage = float(row['coverage'])
                cover_tag = row['type'].split('_')[0]

                # 处理capture数据
                capture_tag = row['capture'].split('_')[0]
                capture_pct = float(row['capture_percentage'])

                # 更新统计
                key = (tid, rtype)
                results[key]['total_len'] += region_len
                results[key]['total_cov'] += coverage * region_len
                results[key]['total_cap'] += capture_pct * region_len

                tag_key = (cover_tag, capture_tag)
                # 核心修复：使用集合去重
                exon_set = results[key]['tags'][tag_key]['exon']
                before_size = len(exon_set)
                exon_set.add(region['exon_number'])
                after_size = len(exon_set)
                # 仅当exon编号是新增时更新count
                if after_size > before_size:
                    results[key]['tags'][tag_key]['count'] = after_size
                results[key]['tags'][tag_key]['exon'].add(region['exon_number'])
    return results

def generate_output(data, transcript_info, gene_order):
    output = []
    sorted_output = []
    for (tid, rtype), stats in data.items():
        gene_info = transcript_info.get(tid, {})
        total_regions = len(transcript_info[tid][rtype])

        # 计算加权平均值
        avg_cov = stats['total_cov'] / stats['total_len'] if stats['total_len'] else 0
        avg_cap = stats['total_cap'] / stats['total_len'] if stats['total_len'] else 0
        cap_str = f"{avg_cap:.2f}"

        # 生成各tag组合
        for (ctag, cptag), tag_stats in stats['tags'].items():
            prop = (tag_stats['count'] / total_regions) * 100
            output.append({
                'gene_id': gene_info.get('gene_id', 'NA'),
                'gene_name': gene_info.get('gene_name', 'NA'),
                'gene_type': gene_info.get('gene_type', 'NA'),
                'transcript_id': tid,
                'region': rtype,
                'cover_tag': ctag,
                'capture_tag': cptag,
                'num': ','.join(map(str, sorted(tag_stats['exon']))),
                'count': tag_stats['count'],
                'total': total_regions,
                'proportion': f"{prop:.2f}%",
                'capture': cap_str,
                'coverage': f"{avg_cov:.2f}"
            })
    # Step 2: 自定义排序逻辑
    def sort_key(item):
        #排序优先级：
       # 1. 文件A中的基因出现顺序
       # 2. 同一基因内exon在前
       # 3. 转录本ID字母顺序
        gene_name=item['gene_name']
        return (
        gene_order.get(gene_name, float('inf')),
            0 if item['region'] == 'exon' else 1,
            item['transcript_id']
        )
        # 执行排序
    sorted_output = sorted(output, key=sort_key)
    return sorted_output

def merge_data(target_dict, source_dict):
    """安全合并两个数据集"""
    for key, src_value in source_dict.items():
        tgt_value = target_dict.setdefault(key, {
            'total_len': 0,
            'total_cov': 0,
            'total_cap': 0,
            'tags': defaultdict(lambda: {'exon': set(), 'count': 0})
        })

        # 合并基础统计量
        tgt_value['total_len'] += src_value['total_len']
        tgt_value['total_cov'] += src_value['total_cov']
        tgt_value['total_cap'] += src_value['total_cap']

        # 合并标签数据
        for tag_key, src_tag in src_value['tags'].items():
            tgt_tag = tgt_value['tags'][tag_key]

            # 合并exon集合（自动去重）
            tgt_tag['exon'].update(src_tag['exon'])

            # 重新计算count（基于合并后的集合大小）
            tgt_tag['count'] = len(tgt_tag['exon'])

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
        description='Process gene and exon coverage data',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )
    
    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')
    
    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--exon', required=True, 
                         metavar='<exon.tsv>',
                         help='Exon coverage and capture input file of the sample')
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
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.exon):
        print(f"Error: Exon file {args.exon} does not exist.")
        return
    
    if not os.path.exists(args.gene):
        print(f"Error: Transcript file {args.gene} does not exist.")
        return
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference file {args.reference} does not exist.")
        return
    
    # 构建输出文件路径
    output_file = f"{args.out}_exon_gene.tsv"
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # 读取全局数据
    print("Reading transcript data...")
    targets, gene_order = read_targets(args.gene)
    print("Reading gene reference data...")
    precise_idx, trans_regions = parse_file_d(args.reference, targets)
    
    # 处理样本数据
    print(f"Processing exon file: {args.exon}")
    exon_data = process_sample(args.exon, 'exon', precise_idx, trans_regions)
    
    # 合并结果
    all_results = defaultdict(lambda: {
        'total_len': 0, 'total_cov': 0, 'total_cap': 0,
        'tags': defaultdict(lambda: {'count': 0, 'exon': set()})
    })
    merge_data(all_results, exon_data)
    
    # 生成输出
    output = generate_output(all_results, trans_regions, gene_order)
    
    # 写入文件
    with open(output_file, 'w', newline='\n', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            '#gene_id', 'gene_name', 'gene_type', 'transcript_id', 'region',
            'cover_tag', 'capture_tag', 'num', 'count', 'total', 'proportion', 'capture', 'coverage'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(output)
    
    print(f"Processing completed: {output_file}")

if __name__ == '__main__':
    main()
