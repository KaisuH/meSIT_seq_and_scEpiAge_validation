#!/bin/bash

# Define the path to the BED file
BED_FILE="./msClockSites_mm10.bed"

# Directory containing BAM files
BAM_DIR="$WORK/data/khiltunen/demuxed_thesis_samples/S183/sample_qc_4/aligned/bam/deduplicated"

# Output directory
OUTPUT_DIR="."

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam;
do
    # Extract the base name of the file for naming the output
    BASE_NAME=$(basename "$BAM_FILE")

    # Run bedtools coverage
    bedtools coverage -a "$BED_FILE" -b "$BAM_FILE" -d > "$OUTPUT_DIR/${BASE_NAME}_coverage.txt"
done

echo "Coverage analysis complete."

