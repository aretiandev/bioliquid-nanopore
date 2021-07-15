# 06 - Create Bams
#
# This script makes sure the clusters do not flip as the sliding window moves to the right. The script is based on 06_read_assignment.ipynb
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
echo ""
echo "06 - CREATE BAMS"
echo ""

# Input variables
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

# Get reads
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
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chrom}_person0.sam" > "strspy/input/${run_number}_${chrom}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools view -b "${run_number}_${chrom}_person1.sam" > "strspy/input/${run_number}_${chrom}_person1.bam"
# Index bam file
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/input/${run_number}_${chrom}_person0.bam"
/home/fer/miniconda3/envs/genomics/bin/samtools index "strspy/input/${run_number}_${chrom}_person1.bam"
echo "Saved: ${datadir}/${run_number}_${chrom}_person0.bam"
echo "Saved: ${datadir}/${run_number}_${chrom}_person1.bam"
echo "Indexed BAM files."