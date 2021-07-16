#!/usr/bin/env bash

# basecaller.sh
# This script runs the guppy_basecaller. Reads should be in the fast5/ folder.
# Output is written to the basecall-latest/ folder.
#
# ARGS:
#   ALGO: algorithm: fast or high accuracy
#   num_callers: number of CPUs to run in parallel
#   cpu_threads_per_caller: number of threads per CPU

ALGO=$1
fast5_folder=$2
output_folder=$3
num_callers=12
cpu_threads_per_caller=1

if [ "$#" -eq 0 ];
then
    echo 'Running the guppy basecaller with the high accuracy algorithm.'
	echo 'Input folder: /home/fer/genomics/fast5'
	echo 'Output folder: /home/fer/genomics/basecall-latest'
    guppy_basecaller --input_path /home/fer/genomics/fast5 --save_path /home/fer/genomics/basecall-latest --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers
elif [ $ALGO == "fast" ];
then
    echo 'Running the guppy basecaller with the fast algorithm.'
	echo 'Input folder:' $fast5_folder 
	echo 'Output folder:' $output_folder
    guppy_basecaller --input_path $fast5_folder --save_path $output_folder --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers -c dna_r10.3_450bps_fast.cfg
elif [ $ALGO == "gpu" ];
then
    echo 'Running the guppy basecaller with the gpu algorithm.'
	echo 'Input folder:' $fast5_folder 
	echo 'Output folder:' $output_folder
    /opt/ont-guppy/bin/guppy_basecaller \
        --input_path $fast5_folder \
        --save_path $output_folder \
        --flowcell FLO-MIN111 --kit SQK-LSK110 \
        --min_qscore 7 -r \
        --records_per_fastq 0 \
        --device 'auto' \
        --chunks_per_runner 512
else
    echo 'Running the guppy basecaller with the high accuracy algorithm.'
	echo 'Input folder: /home/fer/genomics/fast5'
	echo 'Output folder: /home/fer/genomics/basecall-latest'
    guppy_basecaller --input_path /home/fer/genomics/fast5 --save_path /home/fer/genomics/basecall-latest --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers
fi

# Script Templates

# With fast5_out
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --fast5_out --records_per_fastq 0 --cpu_threads_per_caller 8 --num_callers 4

# Fast algorithm: Ran 05-30, completed 05-31
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller 4 --num_callers 16 -c dna_r10.3_450bps_fast.cfg

# High accuracy algorithm: Ran 05-31
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller 4 --num_callers 16
