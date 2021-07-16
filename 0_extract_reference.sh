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
window_width=2000000

case $dis in
    cf)
        echo "Run: ${1}, disease: $dis"
        chrom="chr7"
        location=
        ;;
    sca)
        echo "Run: ${1}, disease: $dis"
        chrom="chr11"
        location=5227002
        ;;
    sma)
        echo "Run: ${1}, disease: $dis"
        chrom="chr5"
        location=
        ;;
    thal)
        echo "Run: ${1}, disease: $dis"
        chrom="chr16"
        location=
        ;;
    pompe)
        echo "Run: ${1}, disease: $dis"
        chrom="chr17"
        location=25000000
        ;;
    *)
        echo "Disease should be in disease list: cf, sca, sma, thal, pompe."
        exit 1
esac

chrom_dis="${chrom}_${dis}"

# Setup
datadir="/mnt/aretian/genomics/nanopore/${run_number}"

# Extract selected chromosome from reference genome
#-------------------------------------------------------------------------------------------
echo "Extracting $chrom from /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
/home/fer/miniconda3/envs/genomics/bin/samtools faidx "/mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" $chrom > "${datadir}/${chrom}_selected.fa"
echo "Saved: ${datadir}/${chrom}_selected.fa"
