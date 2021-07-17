#!/usr/bin/env bash
# 09 - STRspy
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
echo "09 - STRSPY"

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2
echo "Run: ${1}, disease: $dis"
echo ''
source src/setup_variables.sh $dis

# Run STRspy
# -----------------------------------------------------------------------------
echo "Running STRspy."
echo "Config file: /home/fer/genomics/strspy/config/${run_number}_${dis}_inputconfig.txt"

rm -rf "${datadir}/strspy/${dis}/output"
mkdir "${datadir}/strspy/${dis}/output"

cd /home/fer/genomics/strspy
bash STRspy_run_v1.0.sh "config/${run_number}_${dis}_inputconfig.txt" config/UserToolsConfig.txt

echo "Saved STRspy output in: ${datadir}/strspy/${dis}/output"

# Concatenate all output
# -----------------------------------------------------------------------------
echo "Concatenating all STRs in ${datadir}/strspy/${dis}/output/Countings"
cd ${datadir}/strspy/${dis}/output/Countings
cat *person0*Allele_freqs.txt > ${run_number}_${chrom}_person0_strs_raw.txt
cat *person1*Allele_freqs.txt > ${run_number}_${chrom}_person1_strs_raw.txt

grep Human_STR ${run_number}_${chrom}_person0_strs_raw.txt" > ${run_number}_${chrom}_person0_strs.txt
grep Human_STR ${run_number}_${chrom}_person1_strs_raw.txt" > ${run_number}_${chrom}_person1_strs.txt
echo "Saved: ${datadir}/strspy/${dis}/output/Countings/${run_number}_${chrom}_person0_strs.txt"
echo "Saved: ${datadir}/strspy/${dis}/output/Countings/${run_number}_${chrom}_person1_strs.txt"
