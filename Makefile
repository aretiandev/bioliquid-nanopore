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

.PHONY: extract remove_gaps cluster assign create_bams strspy_config strspy tag_reads boolean_matrix str_cluster score

extract:
	@bash 03_extract_location.sh $(run) $(dis)

#remove_gaps: extract
remove_gaps:
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

#cluster: remove_gaps
cluster: 
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

assign:
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_read_assignment.py $(run) $(dis)

#create_bams: assign
create_bams:
	@bash 06_create_bams.sh $(run) $(dis)

#strspy_config: create_bams
strspy_config:
	@/home/fer/miniconda3/envs/genomics/bin/python3 07_strspy_config.py $(run) $(dis)

#strspy: strspy_config
strspy:
	@bash 07_strspy.sh $(run) $(dis)

tag_reads:
	@/usr/bin/Rscript 08_tag_reads.R $(run) $(dis)

boolean_matrix:
	@/usr/bin/Rscript 09_boolean_matrix.R $(run) $(dis)

#str_cluster: boolean_matrix
str_cluster:
	@/usr/bin/Rscript 10_str_clustering.R $(run) $(dis)

#score: str_cluster
score:
	@/home/fer/miniconda3/envs/genomics/bin/python3 11_score.py $(run) $(dis)
