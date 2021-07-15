.PHONY: extract remove_gaps cluster assign

extract:
	@bash 03_extract_location.sh $(run) $(dis)

remove_gaps: extract
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

cluster: remove_gaps
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

assign:
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_read_assignment.py $(run) $(dis)

create_bams: assign
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_create_bams.py $(run) $(dis)
