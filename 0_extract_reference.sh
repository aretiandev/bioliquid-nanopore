#!/usr/bin/env bash
# 0 - Extract Reference
#
# This script extracts a chromosome from the reference genome.
#
# INPUTS:
#   run_number
#   disease
#   Reference Genome
# 
# OUTPUTS:
#   Selected fasta file: chr11_selected.fa
echo ''
echo "0 - EXTRACT REFERENCE"
echo ''

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2
echo "Run: ${1}, disease: $dis"
source src/setup_variables.sh $dis

# Extract selected chromosome from reference genome
#-------------------------------------------------------------------------------------------
echo "Extracting $chrom from /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
/home/fer/miniconda3/envs/genomics/bin/samtools faidx "/mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" $chrom > "${rootdir}/${chrom}_selected.fa"
echo "Saved: ${rootdir}/${chrom}_selected.fa"
