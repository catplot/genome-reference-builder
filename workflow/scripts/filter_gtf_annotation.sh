#!/usr/bin/env bash

# Input: modified GTF file
gtf_modified="${snakemake_params[modified_temp]}"

# Pattern to match Ensembl gene, transcript, and exon IDs for human or mouse
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"

# Step 1: Remove version suffix from transcript, gene, and exon IDs
cat "${snakemake_input[0]}" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for GTF tags
BIOTYPE_PATTERN="(protein_coding|protein_coding_LoF|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Step 2: Construct the gene ID allowlist
gene_allow_list="${snakemake_params[gene_allow_list]}"
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "$gene_allow_list"

# Step 3: Filter the GTF file based on the gene allowlist
gtf_filtered="${snakemake_output[0]}"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
grep -Ff "$gene_allow_list" "$gtf_modified" >> "$gtf_filtered"

# Clean up
rm -rf "$gtf_modified" "$gene_allow_list"
