# Bioliquid Pipeline
#
# Example routine: make extract run=1 dis=sca
#
# Instructions:
#
# There are four blocks of code to be run in order:
#
#   $make cluster   run=1 dis=sca  :  10 mins
#   $make strspy    run=1 dis=sca  :  5 mins
#   $make tag_reads run=1 dis=sca  :  2 mins
#   $make score     run=1 dis=sca  :  1 min
#
# Every other step in between will be run automatically. The reason to split the three blocks is that cluster, strspy and tag_reads take a long time to run so you might want to avoid those depending on the case.

# Set variables
chr           := $(shell bash get_chrom.sh $(dis))
datadir       := /mnt/aretian/genomics/nanopore/run$(run)
extract       := $(datadir)/run$(run)_$(chr)_$(dis).sam
remove_gaps   := $(datadir)/run$(run)_$(chr)_$(dis)_clean.csv
cluster       := $(datadir)/run$(run)_$(chr)_read_clusters.txt
assign        := $(datadir)/run$(run)_$(chr)_person0_uniqueids.txt
create_bams   := $(datadir)/strspy/$(dis)/input/run$(run)_$(chr)_person0.bam
strspy_config := $(datadir)/strspy/$(dis)/input/regions/all_strs.bed
strspy        := $(datadir)/strspy/$(dis)/output
strspy_clean  := $(datadir)/run$(run)_$(chr)_person_full.txt
tag_reads     := $(datadir)/run$(run)_$(chr)_tagged_reads.csv

.PHONY: extract remove_gaps cluster assign create_bams strspy_config strspy strspy_clean tag_reads boolean_matrix str_cluster score print

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
$(tag_reads): $(strspy)
	@/usr/bin/Rscript 08_tag_reads.R $(run) $(dis)

boolean_matrix:
	@/usr/bin/Rscript 09_boolean_matrix.R $(run) $(dis)

#str_cluster: boolean_matrix
str_cluster:
	@/usr/bin/Rscript 10_str_clustering.R $(run) $(dis)

#score: str_cluster
score:
	@/home/fer/miniconda3/envs/genomics/bin/python3 11_score.py $(run) $(dis)
