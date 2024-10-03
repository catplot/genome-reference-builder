"""Snakemake wrapper for Kallisto index."""

__author__ = "Joël Simoneau"
__copyright__ = "Copyright 2019, Joël Simoneau"
__email__ = "simoneaujoel@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get optional extra parameters
extra = snakemake.params.get("extra", "")

# Handle multiple FASTA files
fasta = snakemake.input.get("fasta")
assert fasta is not None, "input-> a FASTA file is required"
fasta = " ".join(fasta) if isinstance(fasta, list) else fasta

# Run Kallisto index command
shell(
    "kallisto index --threads {snakemake.threads} "
    "{extra} "
    "--index {snakemake.output.index} "
    "{fasta} "
    "{log}"
)
