# Bioliquid Nanopore

This repository contains the Nanopore sequencing stage of the Aretian Bioliquid project.

It currently contains the following pipeline:

0. Basecalling
1. Alignment
2. Quality Control
3. Assembly

# Basecalling

The steps for basecalling fast5 files and converting to fastq are as follows:

- Setup a GPU Virtual Machine
- Install Guppy Basecaller
- Run script

For more details view `0_setup.md`