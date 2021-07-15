#!/usr/bin/env bash
#
# Extract location
#
# This file extracts the location of interest for each disease.

# Input variables
run_number=$1
dis=$2

case $dis in
    cf)
        echo "Disease: $dis."
        chrom=""
        location=
        window_width=
        ;;
    sca)
        echo "Disease: $dis."
        chrom="chr11"
        location=5227002
        window_width=2000000
        ;;
    sma)
        echo "Disease: $dis."
        chrom=""
        location=
        window_width=
        ;;
    thal)
        echo "Disease: $dis."
        chrom=""
        location=
        window_width=
        ;;
    pompe)
        echo "Disease: $dis."
        chrom=""
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

# Run Samtools
echo "Extracting location $chrom:$begin-$end into $output."
/home/fer/miniconda3/envs/genomics/bin/samtools view -b $run_reads "chr11:$begin-$end" > $output

# Index the file
/home/fer/miniconda3/envs/genomics/bin/samtools index $output
echo "Done extracting location."