.PHONY: extract
extract:
	@bash 03_extract_location.sh $(run) $(dis)

.PHONY: remove_gaps
remove_gaps: extract
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)
