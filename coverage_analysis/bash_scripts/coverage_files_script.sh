#!/bin/bash

# Define the path to the BED file
BED_FILE="./msClockSites_mm10.bed"

# Parent directory containing subfolders
PARENT_DIR="$WORK/data/khiltunen/demuxed_thesis_samples/S183"

# Output directory
OUTPUT_DIR="."

# Iterate over subfolders sample_qc_1 to sample_qc_4
for SUBFOLDER in sample_qc_{1..4}; do
    # Directory containing BAM files
    BAM_DIR="$PARENT_DIR/$SUBFOLDER/aligned/bam/deduplicated"

    # Loop through each BAM file in the subfolder
    for BAM_FILE in "$BAM_DIR"/*.bam; do
        # Extract the base name of the file for naming the output
        BASE_NAME=$(basename "$BAM_FILE")

        # Run bedtools coverage
        bedtools coverage -a "$BED_FILE" -b "$BAM_FILE" -d > "$OUTPUT_DIR/${BASE_NAME}_coverage.txt"
    done
done

echo "Coverage analysis complete."

