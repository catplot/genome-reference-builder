index_keys = config["index"].keys()

# Salmon index rule
if "salmon" in index_keys:
    rule build_salmon_index:
        input:
            sequences="{outdir}/{database}/{species}/{release}/{genome}.v{release}.transcripts.fa"
        output:
            directory("{outdir}/{database}/{species}/{release}/index/salmon/{salmon_kmer_size}/{genome}")
        log:
            "{outdir}/{database}/{species}/{release}/logs/salmon_index_{genome}_{salmon_kmer_size}.log"
        threads: config["index_threads"]
        params:
            extra=lambda wildcards: f"--kmerLen {wildcards.salmon_kmer_size}"
        conda: "../envs/salmon.yaml"
        message:
            "Building Salmon index for {wildcards.species} ({wildcards.release}) with kmer size {wildcards.salmon_kmer_size}."
        script:
            "../scripts/build_salmon_index.py"

# Kallisto index rule
if "kallisto" in index_keys:
    rule build_kallisto_index:
        input:
            fasta="{outdir}/{database}/{species}/{release}/{genome}.v{release}.transcripts.fa"
        output:
            index="{outdir}/{database}/{species}/{release}/index/kallisto/{kallisto_kmer_size}/{genome}.idx"
        params:
            extra=lambda wildcards: f"--kmer-size {wildcards.kallisto_kmer_size}"
        log:
            "{outdir}/{database}/{species}/{release}/logs/kallisto_index_{genome}_{kallisto_kmer_size}.log"
        threads: config["index_threads"]
        conda: "../envs/kallisto.yaml"
        message:
            "Building Kallisto index for {wildcards.species} ({wildcards.release}) with kmer size {wildcards.kallisto_kmer_size}."
        script:
            "../scripts/build_kallisto_index.py"

# RSEM index rule
if "rsem" in index_keys:
    rule build_rsem_index:
        input:
            reference_genome="{outdir}/{database}/{species}/{release}/{genome}.primary_assembly.genome.fa",
            gtf="{outdir}/{database}/{species}/{release}/{genome}.v{release}.primary_assembly.annotation.filtered.gtf"
        output:
            seq="{outdir}/{database}/{species}/{release}/index/rsem/{genome}.seq",
            grp="{outdir}/{database}/{species}/{release}/index/rsem/{genome}.grp",
            ti="{outdir}/{database}/{species}/{release}/index/rsem/{genome}.ti"
        params:
            extra=""
        log:
            "{outdir}/{database}/{species}/{release}/logs/rsem_index_{genome}.log"
        threads: config["index_threads"]
        conda: "../envs/rsem.yaml"
        message:
            "Building RSEM index for {wildcards.species} ({wildcards.release})."
        script:
            "../scripts/build_rsem_index.py"

# HISAT2 index rule
if "hisat2" in index_keys:
    rule build_hisat2_index:
        input:
            fasta="{outdir}/{database}/{species}/{release}/{genome}.primary_assembly.genome.fa"
        output:
            directory("{outdir}/{database}/{species}/{release}/index/hisat2/{genome}")
        params:
            prefix = "{outdir}/{database}/{species}/{release}/index/hisat2/{genome}/genome"
        log:
            "{outdir}/{database}/{species}/{release}/logs/hisat2_index_{genome}.log"
        threads: config["index_threads"]
        conda: "../envs/hisat2.yaml"
        message:
            "Building HISAT2 index for {wildcards.species} ({wildcards.release})."
        script:
            "../scripts/build_hisat2_index.py"

# STAR index rule
if "star" in index_keys:
    rule build_star_index:
        input:
            fasta="{outdir}/{database}/{species}/{release}/{genome}.primary_assembly.genome.fa",
            gtf="{outdir}/{database}/{species}/{release}/{genome}.v{release}.primary_assembly.annotation.filtered.gtf"
        output:
            directory("{outdir}/{database}/{species}/{release}/index/star/{sjdbOverhang}/{genome}")
        threads: config["index_threads"]
        conda: "../envs/star.yaml"
        params:
            sjdbOverhang=lambda wildcards: wildcards.sjdbOverhang
        log:
            "{outdir}/{database}/{species}/{release}/logs/star_index_{genome}_{sjdbOverhang}.log"
        message:
            "Building STAR index for {wildcards.species} ({wildcards.release}) with sjdbOverhang {wildcards.sjdbOverhang}."
        script:
            "../scripts/build_star_index.py"

# CellRanger index rule
if "cellranger" in index_keys:
    rule build_cellranger_index:
        input:
            fasta="{outdir}/{database}/{species}/{release}/{genome}.primary_assembly.genome.fa",
            gtf="{outdir}/{database}/{species}/{release}/{genome}.v{release}.primary_assembly.annotation.filtered.gtf"
        output:
            directory("{outdir}/{database}/{species}/{release}/index/cellranger/{genome}")
        threads: config["index_threads"]
        params:
            cellranger=config["index"]["cellranger"]["path"]  # Path to CellRanger executable
        log:
            "{outdir}/{database}/{species}/{release}/logs/cellranger_index_{genome}.log"
        message:
            "Building CellRanger index for {wildcards.species} ({wildcards.release}) using CellRanger."
        script:
            "../scripts/build_cellranger_index.py"
