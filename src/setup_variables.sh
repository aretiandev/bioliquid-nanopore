#!/usr/bin/env bash
# SETUP VARIABLES
#
# This script sets up the common variables in all bash scripts.
#
# INPUTS:
#   disease
# 
# OUTPUTS:
#   chr
#   location
#   directories: rootdir, datadir


# Input variables
# -----------------------------------------------------------------------------
window_width=2000000

dis=$1
case $dis in
    cf)
        chr="chr7"
        location=
        ;;
    sca)
        chr="chr11"
        location=5227002
        ;;
    sma)
        chr="chr5"
        location=
        ;;
    thal)
        chr="chr16"
        location=
        ;;
    pompe)
        chr="chr17"
        location=25000000
        ;;
    *)
        echo "Disease should be in disease list: cf, sca, sma, thal, pompe."
        exit 1
esac

chr_dis="${chr}_${dis}"

rootdir="/mnt/aretian/genomics/nanopore"
datadir="/mnt/aretian/genomics/nanopore/${run_number}"
