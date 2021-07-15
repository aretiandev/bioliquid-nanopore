#!/usr/bin/env bash
# 07 - STRspy
#
# This script runs STRspy on the Bioliquid Nanopore data. The script is based on 07_STRspy.ipynb
#
# INPUTS:
#   run_number
#   disease
#   nanopore reads for location: e.g. run1_chr11_sca.sam
#   reference genome fasta for location: e.g. chr11_selected.fa
# 
# OUTPUTS:
#   reads without gaps: e.g. run1_chr11_sca_clean.csv
#   reference genome without gaps: e.g. run1_chr11_reference_genome.json
echo ''
echo "07 - STRSPY"
echo ''

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2

case $dis in
    cf)
        echo "Run: ${1}, disease: $dis"
        chrom="chr7"
        location=
        window_width=
        ;;
    sca)
        echo "Run: ${1}, disease: $dis"
        chrom="chr11"
        location=5227002
        window_width=2000000
        ;;
    sma)
        echo "Run: ${1}, disease: $dis"
        chrom="chr5"
        location=
        window_width=
        ;;
    thal)
        echo "Run: ${1}, disease: $dis"
        chrom="chr16"
        location=
        window_width=
        ;;
    pompe)
        echo "Run: ${1}, disease: $dis"
        chrom="chr17"
        location=
        window_width=
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
echo "Content of Config file /home/fer/genomics/strspy/config/${chrom_dis}_InputConfig.txt"
cat /home/fer/genomics/strspy/config/${chrom_dis}_InputConfig.txt

rm -rf "${datadir}/strspy/output"
mkdir "${datadir}/strspy/output"
cd /home/fer/genomics/strspy
bash STRspy_run_v1.0.sh "config/${chrom_dis}_InputConfig.txt" config/UserToolsConfig.txt
echo "Done."
echo "Saved STRspy output in: ${datadir}/strspy/output"