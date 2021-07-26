#!/usr/bin/env bash
#
# 03 - Extract Reads
#
# This file extracts the location of interest for each disease, converts to bam and sam, and indexes. The script is based on 03_extract_location.ipynb
#
# INPUTS:
#   run_number
#   disease
#   nanopore reads for location: e.g. bioliquid_run2.bam
# 
# OUTPUTS:
#   BAM file: e.g. run1_chr11_sca.bam
#   SAM file: e.g. run1_chr11_sca.sam
echo ""
echo "------------------------------------------------------------------------------------------"
echo "03 - EXTRACT READS (${0})"
echo "Run: ${1}, disease: ${2}"
echo ""
echo ""

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2
source src/setup_variables.sh $dis

# Setup
run_reads="${datadir}/bioliquid_${run_number}.bam"
output="${datadir}/${run_number}_${chrom_dis}.bam"
begin=$(expr $location - $location_padding)
end=$(expr $location + $location_padding)

# Extract location of interest
# -----------------------------------------------------------------------------
echo "Extracting location $chrom:$begin-$end from $run_reads"
/home/fer/miniconda3/envs/genomics/bin/samtools view -b $run_reads "chr11:$begin-$end" > $output
# Index
/home/fer/miniconda3/envs/genomics/bin/samtools index $output

# Convert to Sam
# -----------------------------------------------------------------------------
output="${datadir}/${run_number}_${chrom_dis}.bam"
output_sam="${datadir}/${run_number}_${chrom_dis}.sam"
/home/fer/miniconda3/envs/genomics/bin/samtools view $output > $output_sam
echo "Created $output"
echo "Created $output_sam"
echo "Done extracting location."