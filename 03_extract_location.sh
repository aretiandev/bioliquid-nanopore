#!/usr/bin/env bash
#
# 03 - Extract location
#
# This file extracts the location of interest for each disease, converts to bam and sam, and indexes. The script is based on 03_extract_location.ipynb
#
# INPUTS:
#   run_number
#   disease
#   nanopore reads for location: e.g. run1_chr11_sca.sam
# 
# OUTPUTS:
#   BAM file: e.g. run1_chr11_sca.bam
#   SAM file: e.g. run1_chr11_sca.sam
echo ""
echo "03 - EXTRACT LOCATION"
echo ""

# Input variables
run_number="run${1}"
dis=$2

case $dis in
    cf)
        echo "Disease: $dis"
        chrom="chr7"
        location=
        window_width=
        ;;
    sca)
        echo "Disease: $dis"
        chrom="chr11"
        location=5227002
        window_width=2000000
        ;;
    sma)
        echo "Disease: $dis"
        chrom="chr5"
        location=
        window_width=
        ;;
    thal)
        echo "Disease: $dis"
        chrom="chr16"
        location=
        window_width=
        ;;
    pompe)
        echo "Disease: $dis"
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
run_reads="${datadir}/bioliquid_${run_number}.bam"
output="${datadir}/${run_number}_${chrom_dis}.bam"
begin=$(expr $location - $window_width)
end=$(expr $location + $window_width)

# Extract location of interest
echo "Extracting location $chrom:$begin-$end..."
/home/fer/miniconda3/envs/genomics/bin/samtools view -b $run_reads "chr11:$begin-$end" > $output
# Index
/home/fer/miniconda3/envs/genomics/bin/samtools index $output

# Convert to Sam
output="${datadir}/${run_number}_${chrom_dis}.bam"
output_sam="${datadir}/${run_number}_${chrom_dis}.sam"
/home/fer/miniconda3/envs/genomics/bin/samtools view $output > $output_sam
echo "Created $output"
echo "Created $output_sam"
echo "Done extracting location."