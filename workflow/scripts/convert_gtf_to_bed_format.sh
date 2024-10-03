#!/usr/bin/env bash

# Input and output paths
gtf="${snakemake_input[0]}"
id_mapping="${snakemake_output[0]}"
temp_id_mapping="$id_mapping.temp"

bed12="${snakemake_params[bed12]}"
temp_bed12="$bed12.temp"
chr_list="${snakemake_params[chr_list]}"
bed_gene_ids="${snakemake_params[bed_gene_ids]}"
gene_pred="${snakemake_params[gene_pred]}"

# Step 1: Convert GTF to BED12 format
gtfToGenePred "$gtf" "$gene_pred" && genePredToBed "$gene_pred" "$temp_bed12"

# Step 2: Remove non-standard chromosomes
awk '$1 ~ /^chr[0-9XYM]+$/' "$temp_bed12" > "$bed12"
cut -f1 "$bed12" | sort | uniq > "$chr_list"

# Step 3: Generate gene ID mapping with gene name
awk '$3 == "transcript" {for(i=9;i<=NF;i++){if($i ~ /transcript_id/){transcript_id=$(i+1)}; if($i ~ /gene_name/){gene_name=$(i+1); print transcript_id"\t"gene_name; break;}}}' \
    "$gtf" | tr -d '";' | sort -u > "$temp_id_mapping"
cut -f4 "$bed12" > "$bed_gene_ids"
awk 'NR==FNR{order[$1]=NR; next} ($1 in order){print $0"\t"order[$1]}' \
    "$bed_gene_ids" "$temp_id_mapping" | sort -k3,3n | awk '{print $1"\t"$2}' > "$id_mapping"

# Clean up temporary files
rm -rf "$bed_gene_ids" "$gene_pred" "$temp_id_mapping" "$temp_bed12"
