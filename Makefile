.PHONY: extract remove_gaps clustering

extract:
	@bash 03_extract_location.sh $(run) $(dis)

remove_gaps: extract
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)

clustering: remove_gaps
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)
