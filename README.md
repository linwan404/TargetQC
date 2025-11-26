# TargetQC: A Targeted Quality Control Framework for Clinical Genomic Testing

## Introduction
TargetQC is clinically oriented quality control framework that quantifies exon- and gene-level coverage, evaluates variant calling performance (including accuracy at known pathogenic sites), and supports multiple ES and WGS platforms.

## Installation
```bash
conda env create -f environment.yml
conda activate targetqc
chmod +x TargetQC.sh
```
## Usage
### Target Exon QC
To run the target exon QC module of TargetQC:
```bash
./TargetQC.sh \
  --mode exon \
  --seq ES|WGS \
  --bam <file.bam> \
  [--capture <capture.bed>] \
  --gene <gene.txt> \
  --reference <reference.txt> \
  --out <outprefix>
```
Required Parameters  

- `--mode exon`: Run TargetQC in **Exon capture QC mode**. This mode should be used when performing quality statistics of exon coverage analysis  
- `--seq <ES|WGS>`: Sequencing type, `ES` for exome sequencing, `WGS` for whole genome sequencing  
- `--bam <bam_file>`: Input BAM file(s) of the sample  
- `--capture <capture.bed>`: Capture regions (`.bed`); required for `ES` sequencing type  
- `--gene <gene_file>`: Target genes (e.g., protein-coding genes with clear genetic patterns in OMIM)  
- `--reference <reference_file>`: Gene annotation file 
- `--out <outprefix>`: Output file prefix  
Additional General Options  

- `--thresholds <string>`: Coverage calculation thresholds (must be integers; default: `10,20,30,40`)  
- `--up-exon <int>`: Coverage threshold for well-covered exon regions (default: `30`)  
- `--down-exon <int>`: Coverage threshold for poorly-covered exon regions (default: `20`)  
- `--base-exon <int>`: Number of poorly-covered bases in exons (default: `5`)  
- `--upstream <int>`: Upstream extension length for capture regions (default: `150`)  
- `--downstream <int>`: Downstream extension length for capture regions (default: `150`)  
- `--chrom <string>`: Restrict analysis to regions on this chromosome  
- `-h, --help`: Show this help message
### Target Gene QC
To run the target gene QC module of TargetQC:
```bash
./TargetQC.sh \
  --mode variant-capture \
  --seq ES|WGS \
  --vcf <file.vcf> \
  [--capture <capture.bed>] \
  --pathogenic <pathogenic.txt> \
  --out <outprefix>
```
Required Parameters  
- `--mode gene`: Run TargetQC in **Gene Capture QC mode**. This mode should be used when performing quality statistics of gene coverage analysis  
- `--seq <ES|WGS>`: Sequencing type (`ES` for exome sequencing; `WGS` for whole-genome sequencing)  
- `--bam <bam_file>`: Input BAM file of the sample  
- `--capture <capture.bed>`: Capture regions (`.bed`); **required for `ES` sequencing type**  
- `--gene <gene_file>`: Target genes (e.g., protein-coding genes with clear genetic patterns in OMIM)  
- `--reference <reference_file>`: Gene annotation file
- `--out <outprefix>`: Output file prefix  

Additional General Options  
- `--up-gene <int>`: Coverage threshold for well-covered gene regions (default: `30`)  
- `--down-gene <int>`: Coverage threshold for poorly-covered gene regions (default: `20`)  
- `--well <float>`: Percentage of well-covered regions (default: `95.0`)  
- `--poor <float>`: Percentage of poorly-covered regions (default: `5.0`)  
- `--chrom <string>`: Restrict analysis to regions on this chromosome  
- `-h, --help`: Show this help message
### Variant Calling QC
To run the variant calling QC module of TargetQC:
```bash
./TargetQC.sh \
  --mode variant-calling \
  --seq ES|WGS \
  --vcf <file.vcf> \
  --gvcf <file.gvcf> \
  [--capture <capture.bed>|--ref-bed <ref_bed>] \
  --ref-vcf <ref_vcf> \
  --ref-fasta <ref_fasta> \
  --pathogenic <pathogenic.txt> \
  --out <outprefix>
```
Required Parameters  
- `--mode variant-calling`: Run TargetQC in **Variant Calling QC mode**. This mode should be used when performing quality statistics of variant calling analysis  
- `--seq <ES|WGS>`: Sequencing type (`ES` for exome sequencing; `WGS` for whole-genome sequencing)  
- `--vcf <vcf_file>`: Input VCF file of the sample  
- `--gvcf <gvcf_file>`: Input GVCF file of the sample  
- `--capture <capture.bed>`: Capture regions (`.bed`); required for `ES` sequencing type  
- `--ref-bed <ref_bed>`: Benchmark BED file; required for `WGS` sequencing type
- `--ref-vcf <ref_vcf>`: Reference dataset VCF file (required when `--sample-type` is `standard`)  
- `--ref-fasta <ref_fasta>`: Human reference genome FASTA file  
- `--pathogenic <pathogenic_file>`: Target variant sites
- `--out <outprefix>`: Output file prefix  

Additional General Options  
- `--regions <regions.bed>`: Additional regions for variant calling benchmark (optional)  
- `-h, --help`: Show this help message
### Variant Capture QC
To run the variant capture QC module of TargetQC:
```bash
./TargetQC.sh \
  --mode variant-capture \
  --seq ES|WGS \
  --vcf <file.vcf> \
  [--capture <capture.bed>] \
  --pathogenic <pathogenic.txt> \
  --out <outprefix>
```
Required Parameters  
- `--mode variant-capture`: Run TargetQC in **Variant Capture QC mode**. This mode should be used when performing quality statistics of variant capture analysis  
- `--seq <ES|WGS>`: Sequencing type (`ES` for exome sequencing; `WGS` for whole-genome sequencing)  
- `--vcf <vcf_file>`: Input VCF file of the sample  
- `--capture <capture.bed>`: Capture regions (`.bed`); **required for `ES` sequencing type**  
- `--pathogenic <pathogenic_file>`: Target variant sites
- `--out <outprefix>`: Output file prefix
  
Optional Parameters  
- `--regions <regions.bed>`: Additional regions for variant QC (optional)  
- `-h, --help`: Show this help message
### Quick start
The following commands were used to complete the analysis in one step:
```bash
./TargetQC.sh \
  --sample-type standard|clinical \
  --seq ES|WGS \
  --bam <file.bam> \
  --vcf <file.vcf> \
  [--capture <capture.bed>] \
  --gene <gene.txt> \
  --reference <reference.txt> \
  --pathogenic <pathogenic.txt> \
  --out <outprefix>
```
Required Parameters  
- `--sample-type <standard|clinical>`: **Sample type** (`standard` for benchmarked samples; `clinical` for non-benchmarked clinical samples)  
- `--seq <ES|WGS>`: Sequencing type (`ES` for exome sequencing; `WGS` for whole-genome sequencing)  
- `--bam <bam_file>`: Input BAM file of the sample  
- `--vcf <vcf_file>`: Input VCF file of the sample  
- `--capture <capture.bed>`: Capture regions (`.bed`); **required for `ES` sequencing type**  
- `--gene <gene_file>`: Target genes (e.g., protein-coding genes with clear genetic patterns in OMIM)  
- `--reference <reference_file>`: Gene annotation file
- `--pathogenic <pathogenic_file>`: Target variant sites
- `--out <outprefix>`: Output file prefix  

Additional General Options  
- `--thresholds <string>`: Coverage calculation thresholds (must be integers; default: `10,20,30,40`)  
- `--up-exon <int>`: Coverage threshold for well-covered exon regions (default: `30`)  
- `--down-exon <int>`: Coverage threshold for poorly-covered exon regions (default: `20`)  
- `--base-exon <int>`: Number of poorly-covered bases in exons (default: `5`)  
- `--upstream <int>`: Upstream extension length for capture regions (default: `150`)  
- `--downstream <int>`: Downstream extension length for capture regions (default: `150`)  
- `--up-gene <int>`: Coverage threshold for well-covered gene regions (default: `30`)  
- `--down-gene <int>`: Coverage threshold for poorly-covered gene regions (default: `20`)  
- `--well <float>`: Percentage of well-covered regions (default: `95.0`)  
- `--poor <float>`: Percentage of poorly-covered regions (default: `5.0`)  
- `--chrom <string>`: Restrict analysis to regions on this chromosome  
- `--regions <regions.bed>`: Additional regions for variant detection (optional)  
- `--ref-bed <ref_bed>`: Benchmark BED file (required for `WGS` when `--sample-type` is `standard`)  
- `--ref-vcf <ref_vcf>`: Benchmark VCF file (required when `--sample-type` is `standard`)  
- `--ref-fasta <ref_fasta>`: Human reference genome FASTA file (required when `--sample-type` is `standard`)  
- `--gvcf <gvcf_file>`: Input GVCF file of the sample (**required when `--sample-type` is `standard`**)  
- `-h, --help`: Show this help message

## Input file formats
TargetQC takes as input a BAM, GVCF and VCF file of short read alignments, a gene annotation file and a reference genome, and outputs genotypes in a VCF file. All input files must use the same reference genome version.
### BAM
TargetQC requires a BAM file produced by an indel-sensitive aligner. The BAM file must be sorted and indexed. TargetQC currently only processes a single sample at a time.
### VCF/GVCF
Generated by variant calling software from the BAM file. Refer to the VCF specification for format details.
### Target Genes
TargetQC requires a **target genes file** to specify transcripts and their associated genes for analysis. This is a tab-delimited file with the following columns:

1. Target transcript ID  
2. Gene name  

An example file containing 3 target genes is shown below.  
**NOTE: The table header is for descriptive purposes. The actual file should not have a header.**  

| **Target transcript ID** | **Gene name** |  
|--------------------------|---------------|  
| ENST00000641515.2        | OR4F5         |  
| ENST00000616016.5        | SAMD11        |  
| ENST00000379410.8        | PLEKHN1       |  


#### Key Notes  
- Use **tab separation** between columns (not spaces).  
- Ensure transcript IDs (e.g., `ENST00000641515.2`) and gene names (e.g., `OR4F5`) exactly match definitions in the gene annotation file (`--reference`).  
- This example contains 3 target genes, but the file can include any number of entries.  
- Default file: `target_transcript.tsv` (if not specified via `--gene`).
### Gene annotation
TargetQC requires a **gene annotation file** to define gene/transcript structure details (chromosome, start/end positions, gene/transcript IDs, exon numbers, etc.). You can use the built-in GENCODE v47 annotation or provide a custom file, but the format must strictly follow the specifications below.  

This is a **tab-delimited file** with the following columns:  

1. Chromosome  
2. Gene element (e.g., exon, CDS)  
3. Start position (1-based)  
4. End position (1-based)  
5. Gene ID  
6. Transcript ID  
7. Gene type  
8. Gene name  
9. Exon number (use `NA` for non-exon elements)  

Below is an example file containing 2 genes (truncated for brevity; actual files may include more entries).  

**NOTE: The table header is for descriptive purposes. The actual file should not have a header.**  

| **Chromosome** | **Gene element** | **Start position (1-based)** | **End position (1-based)** | **Gene ID**          | **Transcript ID**        | **Gene type**   | **Gene name** | **Exon number** |  
|----------------|------------------|------------------------------|----------------------------|----------------------|--------------------------|-----------------|---------------|-----------------|  
| chr1           | exon             | 65419                        | 65433                      | ENSG00000186092.7    | ENST00000641515.2        | protein_coding  | OR4F5         | 1               |  
| chr1           | CDS              | 65565                        | 65573                      | ENSG00000186092.7    | ENST00000641515.2        | protein_coding  | OR4F5         | 2               |  
| chr1           | transcript       | 923923                       | 944574                     | ENSG00000187634.13   | ENST00000616016.5        | protein_coding  | SAMD11        | NA               |  
| chr1           | exon             | 923923                       | 924948                     | ENSG00000187634.13   | ENST00000616016.5        | protein_coding  | SAMD11        | 1               |  


#### Key Notes  
- **Tab separation**: Columns must be separated by tabs (not spaces).  
- **1-based coordinates**: Start/end positions use 1-based indexing.  
- **Column order**: Strictly follow the 9-column sequence listed above.  
- **Exon number**: Use `NA` for non-exon elements (e.g., gene, transcript).  
- **Built-in option**: Default file is `genes.gtf` (derived from [GENCODE Release v47](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz)).  
- **Consistency**: Gene/transcript IDs (e.g., `ENSG000001`, `ENST000001`) and gene names (e.g., `GeneA`) must match entries in the **target genes file** (`--gene`).  
- **Scalability**: The file can include any number of genes/transcripts (not limited to the 2-gene example).
### Pathogenic variant sites
TargetQC requires a **pathogenic variant sites file** (0-based coordinates) to identify target variant loci for analysis. You can use the built-in file (derived from the [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar)) or provide a custom file, but the format must strictly follow the BED-like specifications below.  

This is a **BED-like file** with the following columns:  

1. Chromosome  
2. Start position (0-based)  
3. End position (0-based)  
 
Below is an example file containing 2 pathogenic variant loci.  

**NOTE: The table header is for descriptive purposes. The actual BED file should not have a header.**  

| **Chromosome** | **Start position** | **End position** |  
|----------------|--------------------|------------------|  
| chr1           | 943994             | 943995           |  
| chr1           | 964511             | 964512           |  


#### Key Notes  
- **0-based coordinates**: Start/end positions use 0-based indexing (consistent with BED standards).  
- **Tab separation**: Columns must be separated by tabs (not spaces).  
- **Built-in option**: Default file is `clinvar.txt` (derived from the [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar)).  
- **Custom files**: Ensure chromosome names (e.g., `chr1`), start/end positions, and coordinate system (1-based) match your reference genome.  
- **Scalability**: The file can include any number of pathogenic variant loci (not limited to the 2-locus example).

### Capture Regions
BED file defining sequencing capture regions (required for ES data).

## Output File Formats
All outputs are in tab-separated `.tsv` format.
### Exon Capture QC Output Files  
#### `base.tsv`  
Coverage metrics by threshold.  

| **Column number**  | **Description**                                                                 |  
|--------------------|---------------------------------------------------------------------------------|  
| 1                  | Coverage threshold                                                              |  
| 2                  | Target region size                                                             |  
| 3                  | Mean depth in target region                                                     |  
| 4                  | Bases above threshold                                                          |  
| 5                  | Proportion of bases above threshold                                            |  


#### `exon_capture_QC.tsv`  
Per-exon coverage quality.  

| **Column number**  | **Description**                          |  
|--------------------|------------------------------------------|  
| 1                  | Chromosome                               |  
| 2                  | Start position                           |  
| 3                  | End position                             |  
| 4                  | Mean coverage                            |  
| 5                  | Coverage quality label                   |  
| 6                  | Capture label                            |  
| 7                  | Capture region ratio                     |  


#### `exon_capture_QC.summary.tsv`  
Summary matrix of exon QC results.  

### Gene Capture QC Output Files  

#### `trans.tsv`, `trans_cds.tsv`  
Coverage distribution by exon/CDS region in each transcript, describing the proportion of each coverage level within gene exon/CDS regions.  

| **Column number**  | **Description**                          |  
|--------------------|------------------------------------------|  
| 1                  | Gene ID                                  |  
| 2                  | Gene name                                |  
| 3                  | Gene type                                |  
| 4                  | Target transcript ID                     |  
| 5                  | Well-covered region ratio                |  
| 6                  | Middle-covered region ratio              |  
| 7                  | Poorly-covered region ratio              |  

#### `gene_capture_QC.tsv`  
Target genes' coverage summary.  

| **Column number**  | **Description**                          |  
|--------------------|------------------------------------------|  
| 1                  | Gene ID                                  |  
| 2                  | Gene name                                |  
| 3                  | Gene type                                |  
| 4                  | Target transcript ID                     |  
| 5                  | Coverage quality                         |  

#### `gene_capture_QC.summary.tsv`  
Summary matrix of target gene QC results.
### Variant Calling QC Output Files  

#### `variant_calling.summary.tsv`  
Overall variant detection accuracy.  

#### `region_variant_calling.summary.tsv(--regions)`  
Accuracy within target regions.   

| **Column number**  | **Description**       |  
|--------------------|-----------------------|  
| 1                  | Variant type          |  
| 2                  | TP (True Positives)   |  
| 3                  | FN (False Negatives)  |  
| 4                  | FP (False Positives)  |  
| 5                  | Recall                |  
| 6                  | Precision             |  
| 7                  | F1_Score              |  

#### `variant_calling_QC.tsv`  
TP/FP/FN calls with sequencing quality metrics.  
#### `region_variant_calling_QC.tsv(--regions)`  
Accuracy within target regions.   

| **Column number**  | **Description**    |  
|--------------------|--------------------|  
| 1                  | Chromosome         |  
| 2                  | Position           |  
| 3                  | TP/FP/FN type      |  
| 4                  | Variant type       |  
| 5                  | Genotype           |  
| 6                  | DP (Read depth)    |  
| 7                  | BAF (B-allele freq)|  

#### `FN_pathogenic.tsv`  
False negatives in pathogenic variants.  
#### `region_FN_pathogenic.tsv(--regions)`  
Accuracy within target regions.   

- **Generated only when false-negative pathogenic variants exist**.  
| **Column number**  | **Description**    |  
|--------------------|--------------------|  
| 1                  | Chromosome         |  
| 2                  | Position           |
| 3                  | Pathogenic         |

### Variant Capture QC Output Files  

#### `variant_capture_QC.tsv`  
Sequencing quality for pathogenic variants.  
#### `region_variant_capture_QC.tsv(--regions)`  
Sequencing quality for pathogenic variants in target regions. 

| **Column number**  | **Description**    |  
|--------------------|--------------------|  
| 1                  | Chromosome         |  
| 2                  | Position           |  
| 3                  | Variant type       |  
| 4                  | Genotype           |  
| 5                  | DP (Read depth)    |  
| 6                  | BAF (B-allele freq)|  

## Reference Files
- **Target Gene**: Derived from [OMIM](https://www.omim.org/).
- **Gene Annotation**: Derived from [GENCODE Release v47](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz).  
- **Known Pathogenic Variants**: Derived from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar).  
- **NA12878 Benchmark Dataset**: Derived from [Genome in a Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle).This dataset provides the necessary files for the `--ref-bed`, `--ref-vcf`, and `--ref-fasta` parameters. 
## Support
For technical support or questions, please contact 24211240001@m.fudan.edu.cn
