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
num_callers=16
cpu_threads_per_caller=4

if [ "$#" -eq 0 ];
then
    echo 'Running the guppy basecaller with the high accuracy algorithm.'
    guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/basecall-latest --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers
elif [ $ALGO == "fast" ];
then
    echo 'Running the guppy basecaller with the fast algorithm.'
    guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/basecall-latest --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers -c dna_r10.3_450bps_fast.cfg
else
    echo 'Running the guppy basecaller with the high accuracy algorithm.'
    guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/basecall-latest --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller $cpu_threads_per_caller --num_callers $num_callers
fi

# Script Templates

# With fast5_out
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --fast5_out --records_per_fastq 0 --cpu_threads_per_caller 8 --num_callers 4

# Fast algorithm: Ran 05-30, completed 05-31
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller 4 --num_callers 16 -c dna_r10.3_450bps_fast.cfg

# High accuracy algorithm: Ran 05-31
# guppy_basecaller --input_path /home/fer/nanopore/fast5 --save_path /home/fer/nanopore/BaseCall --flowcell FLO-MIN111 --kit SQK-LSK110 --min_qscore 7 -r --records_per_fastq 0 --cpu_threads_per_caller 4 --num_callers 16
