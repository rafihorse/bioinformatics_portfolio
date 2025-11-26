#!/usr/bin/env bash

# Load environments
module load python/3.10.12-fasrc01
mamba activate reliability

# Check if a GFF file is provided as a command-line argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_gff_file>"
    exit 1
fi

# Set the variable for the GFF file
gff_file="$1"

# Sort the GFF file
bedtools sort -i "$gff_file" > "${gff_file%.gff}_sorted.gff"

# Extract relevant terms and create RNA_terms.txt
awk '{print $3}' "$gff_file" | sort -u | grep 'RNA\|transcript' > RNA_terms.txt

# Filter and create a new GFF file
awk 'NR==FNR{a[$0]; next} $3 in a && $1 ~ /^chr/' RNA_terms.txt "${gff_file%.gff}_sorted.gff" \
> "${gff_file%.gff}_sorted_RNA.gff"

# Merge overlapping features and filter for non-overlapping ones
bedtools merge -i "${gff_file%.gff}_sorted_RNA.gff" -c 1 -o count | awk ' { if($4==1) print $0} ' \
> "${gff_file%.gff}_sorted_RNA_no-overlaps.bed"

# Remove overlaps from the original GFF file
bedtools intersect -a "${gff_file%.gff}_sorted_RNA.gff" -b \
"${gff_file%.gff}_sorted_RNA_no-overlaps.bed" -wa \
> "${gff_file%.gff}_sorted_RNA_no-overlaps.gff"

# Get only well supported transcripts
cat "${gff_file%.gff}_sorted_RNA_no-overlaps.gff" | grep transcript_support_level=1
> "${gff_file%.gff}_sorted_RNA_no-overlaps_transcript-support-1.gff"

# Get only protein coding transcripts
cat "${gff_file%.gff}_sorted_RNA_no-overlaps_transcript-support-1.gff" | grep gene_type=protein_coding
> "${gff_file%.gff}_sorted_RNA_no-overlaps_transcript-support-1_protein-coding.gff"

# Unload environment
mamba deactivate