#!/usr/bin/env Rscript

# Render TargetQC_report.Rmd with command-line parameters.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = "") { idx <- which(args == key); if (length(idx) == 0 || idx[[1]] == length(args)) return(default); args[[idx[[1]] + 1]] }
flag_exists <- function(key) key %in% args

rmd <- get_arg("--rmd")
outdir <- path.expand(get_arg("--outdir"))
plot_dir <- path.expand(get_arg("--plot-dir", file.path(outdir, "figures")))
prefix <- get_arg("--prefix")
sample_name <- get_arg("--sample-name", prefix)
sample_type <- get_arg("--sample-type", "clinical")
seq_type <- get_arg("--seq", "ES")
platform <- get_arg("--platform", "")
out_html <- path.expand(get_arg("--out-html", file.path(outdir, paste0(prefix, "_TargetQC_report.html"))))
dp_high <- as.numeric(get_arg("--dp-high", "30"))
dp_low <- as.numeric(get_arg("--dp-low", "20"))
baf_threshold <- as.numeric(get_arg("--baf-threshold", "0.05"))
gene_level_range <- get_arg("--gene-level-range", "≥30X")
selected_gene_type <- get_arg("--selected-gene-type", "protein_coding")

if (!nzchar(rmd) || !file.exists(rmd)) stop("Rmd template not found: ", rmd)
if (!nzchar(outdir) || !dir.exists(outdir)) stop("Output directory not found: ", outdir)
if (!nzchar(prefix)) stop("--prefix is required")
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("R package 'rmarkdown' is required")

dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(input = rmd, output_file = basename(out_html), output_dir = dirname(out_html), params = list(sample_name = sample_name, sample_type = sample_type, assay_type = seq_type, platform = platform, outdir = outdir, plot_dir = plot_dir, prefix = prefix, dp_high = dp_high, dp_low = dp_low, baf_threshold = baf_threshold, gene_level_range = gene_level_range, selected_gene_type = selected_gene_type), envir = new.env(parent = globalenv()), quiet = !flag_exists("--verbose"))
message("[Saved] ", out_html)
