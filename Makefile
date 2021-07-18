# BIOLIQUID PIPELINE
#
# Example routine:
# 
#   $make extract   run=1 dis=sca
#
# There are three routines that take long:
#
#   $make cluster   run=1 dis=sca  :  10 mins
#   $make strspy    run=1 dis=sca  :  36 mins
#   $make tag_reads run=1 dis=sca  :  2 mins
#
# Routines will only be run if dependencies have been updated. The full list is under 'Targets' below.

# Set variables
chrom          := $(shell bash src/get_chrom.sh $(dis))
run_number     := run$(run)
rootdir        := /mnt/aretian/genomics/nanopore
datadir        := /mnt/aretian/genomics/nanopore/$(run_number)

# Targets
get_ref        := $(rootdir)/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
extract_ref    := $(rootdir)/$(chrom)_selected.fa
extract_reads  := $(datadir)/$(run_number)_$(chrom)_$(dis).sam
remove_gaps    := $(datadir)/$(run_number)_$(chrom)_$(dis)_clean.csv
cluster        := $(datadir)/$(run_number)_$(chrom)_read_clusters.txt
assign         := $(datadir)/$(run_number)_$(chrom)_person0_uniqueids.txt
create_bams    := $(datadir)/strspy/$(dis)/input/$(run_number)_$(chrom)_person0.bam
str_list       := $(rootdir)/hg38.hipstr_reference_full_strs.bed
strspy_config  := $(datadir)/strspy/$(dis)/input/regions/all_strs.bed
strspy         := $(datadir)/strspy/$(dis)/output/Countings/run${run}_${chrom}_person0_strs.txt
add_strs       := $(datadir)/$(run_number)_$(chrom)_person_full.txt
tag_reads      := $(datadir)/$(run_number)_$(chrom)_tagged_reads.csv
boolean_matrix := $(datadir)/$(run_number)_$(chrom)_bool_tagged_reads.csv
str_cluster    := $(datadir)/$(run_number)_$(chrom)_kmeans_clusters.csv
score          := $(datadir)/$(run_number)_$(chrom)_recall_score.csv

all: score

.PHONY: basecall align get_ref extract_ref extract remove_gaps cluster assign create_bams strspy_config strspy add_strs tag_reads boolean_matrix str_cluster score

# BASECALL: run from aretian-genomics-gpu server. First get fast5 files from s3.
.PHONY: basecall
basecall:
	@bash 00_basecaller.sh gpu /home/fer/genomics/fast5 /home/fer/genomics/basecall-latest
	@cat /home/fer/genomics/basecall-latest/pass/*fastq > /home/fer/genomics/basecall-latest/bioliquid_$(run_number).fastq

# ALIGN: Run from inside docker container
# $docker run -d -v /home/fer/genomics:/home/jovyan/work -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian:guppy-gpu
# $docker exec -it bioaretian /bin/bash
# $cd ~/work/bioliquid-nanopore
# $make align run=3
.PHONY: align
align:
	@/opt/ont-guppy/bin/minimap2 -x map-ont -a ~/work/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ~/work/basecall-latest/bioliquid_$(run_number).fastq > ~/work/bioliquid_$(run_number).sam
	@/opt/conda/bin/samtools view -bSh  ~/work/bioliquid_$(run_number).sam          > ~/work/bioliquid_$(run_number)_unsorted.bam
	@/opt/conda/bin/samtools sort -@ 32 ~/work/bioliquid_$(run_number)_unsorted.bam > ~/work/bioliquid_$(run_number).bam
	@/opt/conda/bin/samtools index      ~/work/bioliquid_$(run_number).bam
	@/opt/conda/bin/samtools flagstat   ~/work/bioliquid_$(run_number).bam          > ~/work/bioliquid_$(run_number).bam.flag

.PHONY: get_ref
get_ref: $(get_ref) 
$(get_ref): 0_get_reference.sh
	@bash 0_get_reference.sh $(run) $(dis)

.PHONY: extract_ref
extract_ref: $(extract_ref) 
$(extract_ref): $(get_ref) 0_extract_reference.sh
	@bash 0_extract_reference.sh $(run) $(dis)

.PHONY: extract_reads
extract_reads: $(extract_reads) 
$(extract_reads): 03_extract_reads.sh
	@bash 03_extract_reads.sh $(run) $(dis)

.PHONY: remove_gaps
remove_gaps: $(remove_gaps)
$(remove_gaps): $(extract_ref) $(extract_reads) 04_remove_gaps.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

.PHONY: cluster
cluster: $(cluster)
$(cluster): $(remove_gaps) 05_clustering.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

.PHONY: assign
assign: $(assign)
$(assign): $(cluster) 06_assign_reads.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_assign_reads.py $(run) $(dis)

.PHONY: create_bams
create_bams: $(create_bams)
$(create_bams): $(assign) 07_create_bams.sh
	@bash 07_create_bams.sh $(run) $(dis)

.PHONY: strspy_config
strspy_config: $(strspy_config)
$(strspy_config): $(create_bams) 08_strspy_config.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 08_strspy_config.py $(run) $(dis)

.PHONY: strspy
strspy: $(strspy)
$(strspy): $(strspy_config) 09_strspy.sh
	@bash 09_strspy.sh $(run) $(dis)

.PHONY: clean_strspy
clean_strspy:
	rm -rf $(datadir)/strspy/$(dis)/output

.PHONY: str_list
str_list: $(str_list)
$(str_list): 0_str_list.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 0_str_list.py

.PHONY: add_strs
add_strs: $(add_strs)
$(add_strs): $(str_list) $(strspy) 10_add_strs.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 10_add_strs.py $(run) $(dis)

.PHONY: tag_reads
tag_reads: $(tag_reads)
$(tag_reads): $(add_strs) 11_tag_reads.R
	@/usr/bin/Rscript 11_tag_reads.R $(run) $(dis)

.PHONY: boolean_matrix
boolean_matrix: $(boolean_matrix)
$(boolean_matrix): $(tag_reads) 12_boolean_matrix.R
	@/usr/bin/Rscript 12_boolean_matrix.R $(run) $(dis)

.PHONY: str_cluster
str_cluster: $(str_cluster)
$(str_cluster): $(boolean_matrix) 13_str_clustering.R
	@/usr/bin/Rscript 13_str_clustering.R $(run) $(dis)
    @mkdir -p /home/fer/genomics/bioliquid-nanopore/cluster_plots
	@cp $(datadir)/$(run_number)_$(chrom)_assigned_kmeans_clusters.png /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@cp $(datadir)/$(run_number)_$(chrom)_real_sample_labels.png       /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@echo Copying plots to home folder.
	@mkdir -p /home/fer/genomics/bioliquid-nanopore/cluster_plots
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/$(run_number)_$(chrom)_assigned_kmeans_clusters.png
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/$(run_number)_$(chrom)_real_sample_labels.png 

.PHONY: score
score: $(score)
$(score): $(str_cluster) 14_score.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 14_score.py $(run) $(dis)
