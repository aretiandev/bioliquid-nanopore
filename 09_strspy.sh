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
echo ''

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2
echo "Run: ${1}, disease: $dis"
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