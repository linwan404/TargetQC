#!/usr/bin/env python3
"""
benchmark.py - Run hap.py for variant benchmarking
"""

import os
import subprocess
import errno
import argparse

def safe_makedirs(path):
    """创建目录，如果已存在则忽略"""
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def run_benchmark(query_vcf, regions_bed, ref_vcf, ref_fasta, tools_path, python2_path, output_prefix):
    """运行hap.py进行变异检测基准测试"""

    # 检查输入文件存在性
    missing_files = []
    for file_path, file_desc in [
        (query_vcf, "query VCF"),
        (regions_bed, "regions BED"),
        (ref_vcf, "reference VCF"),
        (ref_fasta, "reference FASTA"),
        (tools_path, "hap.py tool"),
        (python2_path, "Python 2 interpreter")
    ]:
        if not os.path.isfile(file_path) and file_desc != "Python 2 interpreter":
            # 对于 Python 解释器，检查是否可执行
            if file_desc == "Python 2 interpreter" and not (os.path.isfile(python2_path) and os.access(python2_path, os.X_OK)):
                missing_files.append((python2_path, file_desc))
            else:
                missing_files.append((file_path, file_desc))

    if missing_files:
        print("Error: Missing files or invalid executables:")
        for file_path, file_desc in missing_files:
            print(f"  {file_desc}: {file_path}")
        return False

    # 创建输出目录
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        safe_makedirs(output_dir)

    # 构建hap.py命令
    cmd = [
        python2_path,
        tools_path,
        ref_vcf,
        query_vcf,
        "-f", regions_bed,
        "-r", ref_fasta,
        "--threads", "5",
        "-o", output_prefix
    ]

    # 执行命令
    try:
        print(f"Processing {query_vcf}...")
        print(f"Command: {' '.join(cmd)}")
        subprocess.check_call(cmd)
        print(f"Completed: {output_prefix}\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to process {query_vcf}: {str(e)}\n")
        return False
    except Exception as e:
        print(f"Unexpected error processing {query_vcf}: {str(e)}\n")
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
        description='Run hap.py for variant benchmarking',
        formatter_class=CustomHelpFormatter,
        add_help=False
    )

    # 手动添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help='show this help message and exit')

    # Required arguments
    required = parser.add_argument_group('Required options')
    required.add_argument('--vcf', required=True,
                         metavar='<file.vcf>',
                         help='Input VCF files of the sample')
    required.add_argument('--regions', required=True,
                         metavar='<regions.bed>',
                         help='BED file for the target region')
    required.add_argument('--tools', required=True,
                         metavar='<hap_py_path>',
                         default='/ocean/home/lanff/miniconda3/envs/happy-env/bin/hap.py',
                         help='The path of hap.py (default: /ocean/home/lanff/miniconda3/envs/happy-env/bin/hap.py)')
    required.add_argument('--out', required=True,
                         metavar='<outprefix>',
                         help='Prefix to name output files')

    # Optional arguments
    optional = parser.add_argument_group('Additional general options')
    optional.add_argument('--ref-vcf',
                         metavar='<ref.vcf>',
                         default='/ocean/home/lanff/reference/variant/NISTv4.2.1/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz',
                         help='Benchmark dataset for the GRCh38/hg38 human reference genome of NA12878 (default: HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz)')
    optional.add_argument('--ref-fasta',
                         metavar='<ref.fasta>',
                         default='/ocean/home/lanff/reference/variant/NISTv4.2.1/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta',
                         help='FASTA file of the GRCh38/hg38 human reference genome (default: GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta)')
    optional.add_argument('--python2-path',
                         metavar='<python2_path>',
                         default='/home/lanff/miniconda3/envs/happy-env/bin/python',
                         help='Path to Python 2 interpreter for running hap.py (default: /usr/bin/python2)')

    args = parser.parse_args()

    # 运行基准测试
    success = run_benchmark(
        args.vcf,
        args.regions,
        args.ref_vcf,
        args.ref_fasta,
        args.tools,
        args.python2_path,
        args.out
    )

    if success:
        print(f"Benchmark completed successfully: {args.out}")
    else:
        print(f"Benchmark failed: {args.out}")

if __name__ == "__main__":
    main()
