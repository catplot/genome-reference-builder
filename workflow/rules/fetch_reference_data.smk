rule get_genome:
    output:
        "{outdir}/{database}/{species}/{release}/{genome}.primary_assembly.genome.fa"
    params:
        datatype="genome"
    log:
        "{outdir}/{database}/{species}/{release}/logs/{genome}.get_genome.log"
    message:
        "Downloading genome for {wildcards.species} from {wildcards.database} (release {wildcards.release})..."
    script:
        "../scripts/fetch_genomic_resources.py"

rule get_transcript:
    output:
        "{outdir}/{database}/{species}/{release}/{genome}.v{release}.transcripts.fa"
    params:
        datatype="transcript"
    log:
        "{outdir}/{database}/{species}/{release}/logs/{genome}.get_transcript.log"
    message:
        "Downloading transcriptome for {wildcards.species} from {wildcards.database} (release {wildcards.release})..."
    script:
        "../scripts/fetch_genomic_resources.py"

rule get_annotation_gtf:
    output:
        "{outdir}/{database}/{species}/{release}/{genome}.v{release}.primary_assembly.annotation.gtf"
    params:
        datatype="annotation.gtf"
    log:
        "{outdir}/{database}/{species}/{release}/logs/{genome}.get_annotation_gtf.log"
    message:
        "Downloading GTF annotation for {wildcards.species} from {wildcards.database} (release {wildcards.release})..."
    script:
        "../scripts/fetch_genomic_resources.py"

rule get_annotation_gff3:
    output:
        "{outdir}/{database}/{species}/{release}/{genome}.v{release}.primary_assembly.annotation.gff3"
    params:
        datatype="annotation.gff3"
    log:
        "{outdir}/{database}/{species}/{release}/logs/{genome}.get_annotation_gff3.log"
    message:
        "Downloading GFF3 annotation for {wildcards.species} from {wildcards.database} (release {wildcards.release})..."
    script:
        "../scripts/fetch_genomic_resources.py"
