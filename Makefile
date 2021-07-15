# Makefile

.PHONY: extract
extract:
	bash 03_extract_location.sh $(run) $(dis)
	
.PHONY: 