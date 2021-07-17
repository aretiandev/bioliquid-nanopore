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
chr            := $(shell bash get_chrom.sh $(dis))
rootdir        := /mnt/aretian/genomics/nanopore
datadir        := /mnt/aretian/genomics/nanopore/run$(run)
# Targets
#basecall      := $(datadir)/bioliquid_run$(run).fastq
#align
get_ref        := $(rootdir)/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
extract_ref    := $(datadir)/$(chr)_selected.fa
extract_reads  := $(datadir)/run$(run)_$(chr)_$(dis).sam
remove_gaps    := $(datadir)/run$(run)_$(chr)_$(dis)_clean.csv
cluster        := $(datadir)/run$(run)_$(chr)_read_clusters.txt
assign         := $(datadir)/run$(run)_$(chr)_person0_uniqueids.txt
create_bams    := $(datadir)/strspy/$(dis)/input/run$(run)_$(chr)_person0.bam
strspy_config  := $(datadir)/strspy/$(dis)/input/regions/all_strs.bed
strspy         := $(datadir)/strspy/$(dis)/output
strspy_clean   := $(datadir)/run$(run)_$(chr)_person_full.txt
tag_reads      := $(datadir)/run$(run)_$(chr)_tagged_reads.csv
boolean_matrix := $(datadir)/run$(run)_$(chr)_bool_tagged_reads.csv
str_cluster    := $(datadir)/run$(run)_$(chr)_kmeans_clusters.csv
score          := $(datadir)/run$(run)_$(chr)_recall_score.csv

all: score

.PHONY: basecall align get_ref extract_ref extract remove_gaps cluster assign create_bams strspy_config strspy strspy_clean tag_reads boolean_matrix str_cluster score

# BASECALL: run from aretian-genomics-gpu server. First get fast5 files from s3.
.PHONY: basecall
#basecall: $(basecall) 
#$(basecall): 00_basecaller.sh
basecall:
	@bash 00_basecaller.sh gpu /home/fer/genomics/fast5 /home/fer/genomics/basecall-latest
	@cat ~/home/fer/genomics/basecall-latest/pass/*fastq > ~/home/fer/genomics/basecall-latest/bioliquid_run$(run).fastq

# ALIGN: Run from inside docker container
# $docker run -d -v /home/fer/genomics:/home/jovyan/work -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian:guppy-gpu
# $docker exec -it bioaretian /bin/bash
# $cd ~/work/bioliquid-nanopore
# $make align
.PHONY: align
#align: $(align) 
#$(align): $(basecall)
align:
	@minimap2 -x map-ont -a ~/work/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ~/work/basecall-latest/bioliquid_run#(run).fastq > ~/work/bioliquid_run$(run).sam
	@samtools view -bSh ~/work/bioliquid_run$(run).sam > ~/work/bioliquid_run$(run)_unsorted.bam
	@samtools sort -@ 32 ~/work/bioliquid_run$(run)_unsorted.bam > ~/work/bioliquid_run$(run).bam
	@samtools index ~/work/bioliquid_run$(run).bam
	@samtools flagstat ~/work/bioliquid_run$(run).bam > ~/work/bioliquid_run$(run).bam.flag

.PHONY: get_ref
get_ref: $(get_ref) 
$(get_ref): 0_get_reference.sh
	@bash 0_get_reference.sh

.PHONY: extract_ref
extract_ref: $(extract_ref) 
$(extract_ref): $(get_ref) 0_extract_reference.sh
	@bash 0_extract_reference.sh $(run) $(dis)

.PHONY: extract_reads
extract_reads: $(extract_reads) 
$(extract_reads): $(extract_reads) 03_extract_reads.sh
	@bash 03_extract_reads.sh $(run) $(dis)

.PHONY: remove_gaps
remove_gaps: $(remove_gaps)
$(remove_gaps): $(extract) 04_remove_gaps.py
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
	rm -rf $(strspy)

.PHONY: strspy_clean
strspy_clean: $(strspy_clean)
$(strspy_clean): $(strspy) 10_strspy_clean.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 10_strspy_clean.py $(run) $(dis)

.PHONY: tag_reads
tag_reads: $(tag_reads)
$(tag_reads): $(strspy_clean) 11_tag_reads.R
	@/usr/bin/Rscript 11_tag_reads.R $(run) $(dis)

.PHONY: boolean_matrix
boolean_matrix: $(boolean_matrix)
$(boolean_matrix): $(tag_reads) 12_boolean_matrix.R
	@/usr/bin/Rscript 12_boolean_matrix.R $(run) $(dis)

.PHONY: str_cluster
str_cluster: $(str_cluster)
$(str_cluster): $(boolean_matrix) 13_str_clustering.R
	@/usr/bin/Rscript 13_str_clustering.R $(run) $(dis)
	@cp $(datadir)/run$(run)_$(chr)_assigned_kmeans_clusters.png /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@cp $(datadir)/run$(run)_$(chr)_real_sample_labels.png       /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@echo Copying plots to home folder.
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/run$(run)_$(chr)_assigned_kmeans_clusters.png
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/run$(run)_$(chr)_real_sample_labels.png 

.PHONY: score
score: $(score)
$(score): $(str_cluster) 14_score.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 14_score.py $(run) $(dis)
