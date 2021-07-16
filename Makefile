# Bioliquid Pipeline
#
# Example routine:
# 
#   $make extract run=1 dis=sca
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
datadir        := /mnt/aretian/genomics/nanopore/run$(run)
# Targets
extract        := $(datadir)/run$(run)_$(chr)_$(dis).sam
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

.PHONY: extract remove_gaps cluster assign create_bams strspy_config strspy strspy_clean tag_reads boolean_matrix str_cluster score

.PHONY: extract
extract: $(extract)
$(extract):
	@bash 03_extract_location.sh $(run) $(dis)

.PHONY: remove_gaps
remove_gaps: $(remove_gaps)
$(remove_gaps): $(extract)
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

.PHONY: cluster
cluster: $(cluster)
$(cluster): $(remove_gaps)
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

.PHONY: assign
assign: $(assign)
$(assign): $(cluster)
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_assign_reads.py $(run) $(dis)

.PHONY: create_bams
create_bams: $(create_bams)
$(create_bams): $(assign)
	@bash 06_create_bams.sh $(run) $(dis)

.PHONY: strspy_config
strspy_config: $(strspy_config)
$(strspy_config): $(create_bams)
	@/home/fer/miniconda3/envs/genomics/bin/python3 07_strspy_config.py $(run) $(dis)

.PHONY: strspy
strspy: $(strspy)
$(strspy): $(strspy_config)
	@bash 07_strspy.sh $(run) $(dis)

.PHONY: strspy_clean
strspy_clean: $(strspy_clean)
$(strspy_clean): $(strspy)
	@/home/fer/miniconda3/envs/genomics/bin/python3 07_strspy_clean.py $(run) $(dis)

.PHONY: tag_reads
tag_reads: $(tag_reads)
$(tag_reads): $(strspy_clean)
	@/usr/bin/Rscript 08_tag_reads.R $(run) $(dis)

.PHONY: boolean_matrix
boolean_matrix: $(boolean_matrix)
$(boolean_matrix): $(tag_reads)
	@/usr/bin/Rscript 09_boolean_matrix.R $(run) $(dis)

.PHONY: str_cluster
str_cluster: $(str_cluster)
$(str_cluster): $(boolean_matrix)
	@/usr/bin/Rscript 10_str_clustering.R $(run) $(dis)

.PHONY: score
score: $(score)
$(score): $(str_cluster)
	@/home/fer/miniconda3/envs/genomics/bin/python3 11_score.py $(run) $(dis)
