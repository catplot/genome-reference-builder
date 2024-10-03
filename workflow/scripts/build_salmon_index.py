"""Snakemake wrapper for Salmon index."""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2018, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os.path import dirname
from snakemake.shell import shell
from tempfile import TemporaryDirectory

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

# Handle decoys if provided
decoys = snakemake.input.get("decoys", "")
if decoys:
    decoys = f"--decoys {decoys}"

# Handle output directory
output = snakemake.output
if len(output) > 1:
    output = dirname(snakemake.output[0])

# Run Salmon index command inside a temporary directory
with TemporaryDirectory() as tempdir:
    shell(
        "salmon index "
        "--transcripts {snakemake.input.sequences} "
        "--index {output} "
        "--threads {snakemake.threads} "
        "--tmpdir {tempdir} "
        "{decoys} "
        "{extra} "
        "{log}"
    )
