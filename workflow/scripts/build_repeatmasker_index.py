"""Snakemake wrapper for RepeatMasker index."""

__author__ = "Songqi Duan | 段松岐"
__copyright__ = "Copyright 2024, Songqi Duan | 段松岐"
__email__ = "duan@songqi.org"
__license__ = "MIT"

import os
from snakemake.shell import shell

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get optional extra parameters and species
extra = snakemake.params.get("extra", "")
species = snakemake.params.get("species", None)
assert species is not None, "params-> a species is required for RepeatMasker"

# Handle the input FASTA file for the genome
fasta = snakemake.input.get("fasta")
assert fasta is not None, "input-> a FASTA file is required"
input_seq = ",".join(fasta) if isinstance(fasta, list) else fasta

# Create the output directory if it doesn't exist
repeatmasker_dir = snakemake.output[0]
os.makedirs(repeatmasker_dir, exist_ok=True)

# Run RepeatMasker command
shell(
    "RepeatMasker {extra} "
    "-species {species} "
    "-pa {snakemake.threads} "
    "-dir {repeatmasker_dir} "
    "{input_seq} "
    "{log}"
)