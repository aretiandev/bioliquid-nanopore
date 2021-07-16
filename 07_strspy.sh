#!/usr/bin/env bash
# 07 - STRspy
#
# This script runs STRspy on the Bioliquid Nanopore data. The script is based on 07_STRspy.ipynb
#
# INPUTS:
#   run_number
#   disease
#   BAM, BED and fasta files needed to run STRspy
# 
# OUTPUTS:
#   Identified STRs: e.g. {datadir}/strspy/${dis}/output/Countings
echo ''
echo "07 - STRSPY"
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
cd $datadir

# Run STRspy
# -----------------------------------------------------------------------------
echo "Running STRspy."
echo "Config file: /home/fer/genomics/strspy/config/${dis}_InputConfig.txt"

rm -rf "${datadir}/strspy/${dis}/output"
mkdir "${datadir}/strspy/${dis}/output"
cd /home/fer/genomics/strspy
bash STRspy_run_v1.0.sh "config/${dis}_InputConfig.txt" config/UserToolsConfig.txt

echo "Saved STRspy output in: ${datadir}/strspy/${dis}/output"