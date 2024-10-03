"""Snakemake wrapper for STAR index."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import tempfile
from snakemake.shell import shell

# Create log for stdout and stderr
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

# Handle sjdbOverhang parameter if provided
sjdb_overhang = snakemake.params.get("sjdbOverhang", "")
if sjdb_overhang:
    sjdb_overhang = f"--sjdbOverhang {sjdb_overhang}"

# Handle GTF file if provided
gtf = snakemake.input.get("gtf", "")
if gtf:
    gtf = f"--sjdbGTFfile {gtf}"

# Run STAR index command inside a temporary directory
with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "STAR --runThreadN {snakemake.threads} "
        "--runMode genomeGenerate "
        "--genomeFastaFiles {snakemake.input.fasta} "
        "{sjdb_overhang} "
        "{gtf} "
        "{extra} "
        "--outTmpDir {tmpdir}/STARtmp "
        "--genomeDir {snakemake.output} "
        "{log}"
    )
