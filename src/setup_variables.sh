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
    sma)
        chrom="chr5"
        location=
        ;;
    thal)
        chrom="chr16"
        location=
        ;;
    pompe)
        chrom="chr17"
        location=25000000
        ;;
    *)
        echo "Disease should be in disease list: cf, sca, sma, thal, pompe."
        exit 1
esac

chrom_dis="${chrom}_${dis}"

rootdir="/mnt/aretian/genomics/nanopore"
datadir="/mnt/aretian/genomics/nanopore/${run_number}"
