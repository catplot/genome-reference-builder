# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-10-03

### Added

- **Initial Release** of the **Genome Reference Builder** Snakemake pipeline.
- Automated genome indexing for:
  - **Human**: GENCODE Release 46 (GRCh38.p14)
  - **Mouse**: GENCODE Release M35 (GRCm39)
- Integration of essential bioinformatics tools:
  - `STAR` (v2.7.11b)
  - `RSEM` (v1.3.3)
  - `Cell Ranger` (v8.0.1)
  - `HISAT2` (v2.2.1)
  - `Kallisto` (v0.51.0)
  - `Salmon` (v1.10.1)
- GTF filtering step based on the 10x Genomics filtering scheme to ensure consistency between Bulk RNA-seq and Single-cell RNA-seq analyses.
- Configuration file (`config.yaml`) pre-configured with May 2024 GENCODE versions for human and mouse genomes.
- Support for both direct execution and SLURM job scheduler submission.
- Comprehensive `README.md` with detailed instructions on installation, configuration, usage, and troubleshooting.
- Project structure setup, including directories for configuration, logs, and output.
