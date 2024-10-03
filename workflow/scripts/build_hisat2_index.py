"""Snakemake wrapper for HISAT2 index."""

__author__ = "Joël Simoneau"
__copyright__ = "Copyright 2019, Joël Simoneau"
__email__ = "simoneaujoel@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get optional extra parameters
extra = snakemake.params.get("extra", "")

# Handle multiple FASTA files
fasta = snakemake.input.get("fasta")
assert fasta is not None, "input-> a FASTA file is required"
input_seq = ",".join(fasta) if isinstance(fasta, list) else fasta

# Create the output directory if it doesn't exist
hisat_dir = snakemake.output[0]
os.makedirs(hisat_dir, exist_ok=True)

# Run HISAT2 index command
shell(
    "hisat2-build {extra} "
    "-p {snakemake.threads} "
    "{input_seq} "
    "{hisat_dir} "
    "{log}"
)
