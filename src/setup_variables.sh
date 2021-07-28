#!/usr/bin/env bash
# SETUP VARIABLES
#
# This script sets up the common variables in all bash scripts.
#
# INPUTS:
#   disease
# 
# OUTPUTS:
#   chrom
#   location
#   directories: rootdir, datadir


# Input variables
# -----------------------------------------------------------------------------
location_padding=2000000

dis=$1

case $dis in
    cf)
        chrom="chr7"
        location=117559590
        ;;
    sca)
        chrom="chr11"
        location=5227002
        ;;
    sma1)
        chrom="chr5"
        location=$(( (70924941 + 70953015)/2 ))
        ;;
    sma2)
        chrom="chr5"
        location=$(( (70049523 + 70077595)/2 ))
        ;;
    thal1)
        chrom="chr16"
        location=$(( (176680   + 177522)/2  )) 
        ;;
    thal2)
        chrom="chr16"
        location=$(( (172876   + 173710)/2  ))
        ;;
    thal3)
        chrom="chr11"
        location=$(( (5225464  + 5227071)/2 ))
        ;;
    pompe)
        chrom="chr17"
        location=25000000
        ;;
    *)
        echo "Disease should be in disease list: cf, sca, sma1/2, thal1/2/3, pompe."
        exit 1
esac

chrom_dis="${chrom}_${dis}"

rootdir="/mnt/aretian/genomics/nanopore"
datadir="/mnt/aretian/genomics/nanopore/${run_number}"
