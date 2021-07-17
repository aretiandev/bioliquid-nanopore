#!/usr/bin/env bash
# 07 - Create Bams
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
echo "07 - CREATE BAMS"
echo ''

# Input variables
# -----------------------------------------------------------------------------
run_number="run${1}"
dis=$2
echo "Run: ${1}, disease: $dis"
source src/setup_variables.sh $dis

# Create bam files
# -----------------------------------------------------------------------------
cd $datadir
echo "Extracting SAM files for person0 and person1 out of ${run_number}_${chr_dis}.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools view "${run_number}_${chr_dis}.bam" | grep -f "${run_number}_${chr}_person0_uniqueids.txt" > "${run_number}_${chr}_person0_reads.sam"
/home/fer/miniconda3/envs/genomics/bin/samtools view "${run_number}_${chr_dis}.bam" | grep -f "${run_number}_${chr}_person1_uniqueids.txt" > "${run_number}_${chm}_person1_reads.sam"
# Get header
/home/fer/miniconda3/envs/genomics/bin/samtools view -H "${run_number}_${chr_dis}.bam" > "${run_number}_${chr}_person_header.txt"
# Concatenate
cat "${run_number}_${chr}_person_header.txt" "${run_number}_${chr}_person0_reads.sam" > "${run_number}_${chr}_person0.sam"
cat "${run_number}_${chr}_person_header.txt" "${run_number}_${chr}_person1_reads.sam" > "${run_number}_${chr}_person1.sam"
echo "Saved: ${datadir}/${run_number}_${chr}_person0.sam"
echo "Saved: ${datadir}/${run_number}_${chr}_person1.sam"
# Convert to bam file
mkdir -p strspy/${dis}/input
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chr}_person0.sam" > "strspy/${dis}/input/${run_number}_${chr}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chr}_person1.sam" > "strspy/${dis}/input/${run_number}_${chr}_person1.bam"
# Index bam file
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/${dis}/input/${run_number}_${chr}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/${dis}/input/${run_number}_${chr}_person1.bam"
echo "Saved: ${datadir}/strspy/${dis}/input/${run_number}_${chr}_person0.bam"
echo "Saved: ${datadir}/strspy/${dis}/input/${run_number}_${chr}_person1.bam"
echo "Indexed BAM files."