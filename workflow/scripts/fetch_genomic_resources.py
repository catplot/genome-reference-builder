__author__ = "Songqi Duan"
__copyright__ = "Copyright (C) 2024 by Songqi Duan | 段松岐"
__email__ = "songqi.duan@outlook.com"
__license__ = "MIT"

import subprocess as sp
import sys
from snakemake.shell import shell

# Fetch wildcards from Snakemake
species = snakemake.wildcards.species.lower()  # Get species in lowercase
release = snakemake.wildcards.release          # Release version
genome = snakemake.wildcards.genome            # Genome identifier
database = snakemake.wildcards.database        # Database source

# Setup log file for capturing error messages
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Base URL for downloading data from GenCode FTP server
url_prefix = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{release}"

# Determine file suffixes based on the data type (genome, transcript, annotation)
datatype = snakemake.params.get("datatype", "")
if datatype == "genome":
    suffixes = [f"{genome}.primary_assembly.genome.fa.gz"]
elif datatype == "transcript":
    suffixes = [f"gencode.v{release}.transcripts.fa.gz"]
elif datatype == "annotation.gtf":
    suffixes = [f"gencode.v{release}.primary_assembly.annotation.gtf.gz"]
elif datatype == "annotation.gff3":
    suffixes = [f"gencode.v{release}.primary_assembly.annotation.gff3.gz"]
else:
    raise ValueError("Invalid datatype, must be one of genome, transcript, annotation.gtf, or annotation.gff3")

# Flag to track if download is successful
success = False

# Try to download files by iterating over possible suffixes
for suffix in suffixes:
    url = f"{url_prefix}/{suffix}"
    try:
        # Check if URL is accessible
        shell("curl -sSfI {url} > /dev/null 2> /dev/null")
    except sp.CalledProcessError:
        # If URL check fails, continue to the next suffix
        continue

    # Download and decompress the file, logging output and errors
    shell("(curl -L {url} | zcat >> {snakemake.output[0]}) {log}")
    success = True  # Mark download as successful
    break  # Break loop as soon as a file is successfully downloaded

# Handle the case where no download was successful
if not success:
    url = f"{url_prefix}/{suffixes[0]}"
    print(
        f"Unable to download requested sequence data from GenCode ({url}). "
        "Please check if the above URL is accessible (might be a temporary server issue). "
        "Also, ensure that the specified combination of species, genome, and release is provided by GenCode.",
        file=sys.stderr,
    )
    exit(1)  # Exit with error if no valid data was downloaded
