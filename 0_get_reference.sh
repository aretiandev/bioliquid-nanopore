#!/usr/bin/env bash
# 0 - Get Reference
#
# Get Reference Genome
#
# INPUTS:
#   run_number
#   disease
# 
# OUTPUTS:
#   Reference Genome
echo ''
echo "0 - GET REFERENCE"
echo ''

# Download reference genome
#-------------------------------------------------------------------------------------------
# Pull reference genome from S3
echo "Downloading reference genome from S3"
s3cmd get --force s3://aretian-genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz     /mnt/aretian/genomics/nanopore/
s3cmd get --force s3://aretian-genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.fai /mnt/aretian/genomics/nanopore/
s3cmd get --force s3://aretian-genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.gzi /mnt/aretian/genomics/nanopore/
touch /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
touch /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.fai
touch /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.gzi
echo "Saved: /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
echo "Saved: /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.fai"
echo "Saved: /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz.gzi"

# Not needed
#gunzip /mnt/aretian/genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
