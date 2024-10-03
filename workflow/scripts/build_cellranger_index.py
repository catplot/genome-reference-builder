"""Snakemake wrapper for CellRanger index."""

__author__ = "Songqi Duan | 段松岐"
__copyright__ = "Copyright 2024, Songqi Duan | 段松岐"
__email__ = "duan@songqi.org"
__license__ = "MIT"

import os
from snakemake.shell import shell

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get paths to input fasta and gtf files
fasta = snakemake.input.get("fasta")
gtf = snakemake.input.get("gtf")

# Ensure that the input files are provided
assert fasta is not None, "input-> a FASTA file is required"
assert gtf is not None, "input-> a GTF file is required"

# Get output directory
output_dir = snakemake.output[0]

# Get the path to CellRanger from the params
cellranger_path = snakemake.params.get("cellranger", "cellranger")
assert cellranger_path is not None, "params-> cellranger path must be provided"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Run the CellRanger mkref command
shell(
    "{cellranger_path} mkref "
    "--genome=genome "
    "--fasta={fasta} "
    "--genes={gtf} "
    "--localmem 30 "
    "--nthreads {snakemake.threads} "
    "> {log} 2>&1"
)
