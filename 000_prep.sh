# PREP script
# This script should be run before everything else
# It creates the base files that are inputs to the pipeline
# - Get reference genome and index
# - Extract locations of interest and index

# FROM 04_remove_gaps

# Pull reference genome from S3
# !s3cmd get s3://aretian-genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# !gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Index reference genome
# !samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

%%bash
echo "Extracting $chrom from ${datadir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
echo "Writing to ${datadir}/${chrom}_selected.fa"

/home/fer/miniconda3/envs/genomics/bin/samtools faidx "${datadir}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" $chrom > "${datadir}/${chrom}_selected.fa"