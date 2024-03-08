#!/bin/bash

# Define the path to the BED file
BED_FILE="./msClockSites_mm10.bed"

# Parent directory containing BAM files
PARENT_DIR="./S183_bam_coverage"

# Output directory
OUTPUT_DIR="./output_coverage_analysis"

# Loop through each BAM file in the directory
for BAM_FILE in "$PARENT_DIR"/*.bam; do
    # Extract the base name of the file for naming the output
    BASE_NAME=$(basename "$BAM_FILE" .bam)

    # Run bedtools coverage
    bedtools coverage -a "$BED_FILE" -b "$BAM_FILE" -d > "$OUTPUT_DIR/${BASE_NAME}_coverage.txt"
    
    # Plot coverage using gnuplot
    gnuplot -persist <<-EOF
        set terminal pngcairo enhanced font "arial,10" size 800,600
        set output "$OUTPUT_DIR/${BASE_NAME}_coverage_plot.png"
        set title "${BASE_NAME} Coverage Plot"
        set xlabel "Genomic Position"
        set ylabel "Coverage"
        plot "$OUTPUT_DIR/${BASE_NAME}_coverage.txt" using 2:3 with lines title "Sample QC 1"
EOF
done

echo "Coverage analysis and plotting complete."

