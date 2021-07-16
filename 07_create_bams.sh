#!/usr/bin/env bash
# 06 - Create Bams
#
# This script extracts the reads for each cluster and creates two BAM files. The script is based on 06_read_assignment.ipynb
#
# INPUTS:
#   run_number
#   disease
#   List of read IDs for each person: e.g. run1_chr11_person0_uniqueids.txt
#   nanopore reads for location: e.g. run1_chr11_sca.bam
# 
# OUTPUTS:
#   BAM files: e.g. run1_chr11_person0.bam" 
echo ''
echo "06 - CREATE BAMS"
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

# Create bam files
# -----------------------------------------------------------------------------
echo "Extracting SAM files for person0 and person1 out of ${run_number}_${chrom_dis}.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools view "${run_number}_${chrom_dis}.bam" | grep -f "${run_number}_${chrom}_person0_uniqueids.txt" > "${run_number}_${chrom}_person0_reads.sam"
/home/fer/miniconda3/envs/genomics/bin/samtools view "${run_number}_${chrom_dis}.bam" | grep -f "${run_number}_${chrom}_person1_uniqueids.txt" > "${run_number}_${chrom}_person1_reads.sam"
# Get header
/home/fer/miniconda3/envs/genomics/bin/samtools view -H "${run_number}_${chrom_dis}.bam" > "${run_number}_${chrom}_person_header.txt"
# Concatenate
cat "${run_number}_${chrom}_person_header.txt" "${run_number}_${chrom}_person0_reads.sam" > "${run_number}_${chrom}_person0.sam"
cat "${run_number}_${chrom}_person_header.txt" "${run_number}_${chrom}_person1_reads.sam" > "${run_number}_${chrom}_person1.sam"
echo "Saved: ${datadir}/${run_number}_${chrom}_person0.sam"
echo "Saved: ${datadir}/${run_number}_${chrom}_person1.sam"
# Convert to bam file
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chrom}_person0.sam" > "strspy/${dis}/input/${run_number}_${chrom}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chrom}_person1.sam" > "strspy/${dis}/input/${run_number}_${chrom}_person1.bam"
# Index bam file
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/${dis}/input/${run_number}_${chrom}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/${dis}/input/${run_number}_${chrom}_person1.bam"
echo "Saved: ${datadir}/strspy/${dis}/input/${run_number}_${chrom}_person0.bam"
echo "Saved: ${datadir}/strspy/${dis}/input/${run_number}_${chrom}_person1.bam"
echo "Indexed BAM files."