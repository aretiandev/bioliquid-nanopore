# Bioliquid Pipeline
#
# Instructions:
#
# There are three blocks of code to be run in order:
#
#   $make cluster
#   $make strspy
#   $make str_cluster
#
# Every other step in between will be run automatically. The reason to split the three blocks is that 05_clustering.py and 07_strspy.sh take a long time to run so you might want to avoid those depending on the case

.PHONY: extract remove_gaps cluster assign create_bams strspy_config strspy tag_reads boolean_matrix str_cluster score

extract:
	@bash 03_extract_location.sh $(run) $(dis)

remove_gaps: extract
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

cluster: remove_gaps
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

assign:
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_read_assignment.py $(run) $(dis)

create_bams: assign
	@bash 06_create_bams.sh $(run) $(dis)

strspy_config: create_bams
	@/home/fer/miniconda3/envs/genomics/bin/python3 07_strspy_config.py $(run) $(dis)

strspy: 
	@bash 07_strspy.sh $(run) $(dis)

tag_reads:
	@/usr/bin/Rscript 08_tag_reads.R $(run) $(dis)

boolean_matrix: tag_reads
	@/usr/bin/Rscript 09_boolean_matrix.R $(run) $(dis)

str_cluster: boolean_matrix
	@/usr/bin/Rscript 10_str_clustering.R $(run) $(dis)

score: 
	@/home/fer/miniconda3/envs/genomics/bin/python3 11_score.py $(run) $(dis)
