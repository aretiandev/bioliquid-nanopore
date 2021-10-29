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
dis            ?= sca
chrom          := $(shell bash src/get_chrom.sh $(dis))
chrom_dis      := $(chrom)_$(dis)
run_number     := run$(run)
rootdir        := /mnt/aretian/genomics/nanopore
datadir        := /mnt/aretian/genomics/nanopore/$(run_number)

# Targets
get_ref        := $(rootdir)/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
extract_ref    := $(rootdir)/reference_genome_$(chrom).fa
extract_reads  := $(datadir)/$(run_number)_$(chrom_dis).sam
remove_gaps    := $(datadir)/$(run_number)_$(chrom_dis)_clean.csv
cluster        := $(datadir)/$(run_number)_$(chrom_dis)_read_clusters.txt
assign         := $(datadir)/$(run_number)_$(chrom_dis)_person0_uniqueids.txt
create_bams    := $(datadir)/strspy/$(dis)/input/$(run_number)_$(chrom_dis)_person0.bam
str_list       := $(rootdir)/hg38.hipstr_reference_full_strs.bed
strspy_config  := $(datadir)/strspy/$(dis)/input/regions/all_strs.bed
strspy         := $(datadir)/strspy/$(dis)/output/Countings/run${run}_${chrom_dis}_person0_strs.txt
add_strs       := $(datadir)/$(run_number)_$(chrom_dis)_person_full.txt
tag_reads      := $(datadir)/$(run_number)_$(chrom_dis)_tagged_reads.csv
long_reads     := $(datadir)/$(run_number)_$(chrom_dis)_long_tagged_reads.csv
boolean_matrix := $(datadir)/$(run_number)_$(chrom_dis)_bool_tagged_reads.csv
str_cluster    := $(datadir)/$(run_number)_$(chrom_dis)_kmeans_clusters.csv
score          := $(datadir)/$(run_number)_$(chrom_dis)_recall_score.csv
sample_id      := $(datadir)/$(run_number)_$(chrom_dis)_sample_id_read_shares.csv

.PHONY: all all_log basecall align qc get_ref extract_ref extract_reads clean_extract_reads remove_gaps clean_remove_gaps cluster clean_cluster assign create_bams str_list strspy_config strspy clean_strspy add_strs tag_reads long_reads boolean_matrix str_cluster score sample_id disease_diag

all: disease_diag

all_log:
	@mkdir -p logs
	@echo "BIOLIQUID PIPELINE"
	@echo "Run: $(run), disease: $(dis)"
	@echo ""
	@echo "Log file: logs/$(run_number)_$(dis)_$(shell date '+%Y-%m-%d-%Hh').log"
	@$(MAKE) all run=$(run) dis=$(dis) 2>&1 | tee logs/$(run_number)_$(dis)_$(shell date '+%Y-%m-%d-%Hh').log

# BASECALL: run from GPU server. First get fast5 files into ~/genomics/fast5
# $make basecall run=6
basecall:
	bash 00_basecaller.sh gpu ~/genomics/fast5 ~/genomics/basecall-latest
	cat ~/genomics/basecall-latest/pass/*fastq > ~/genomics/basecall-latest/bioliquid_$(run_number).fastq

# ALIGN: run from GPU server.
# $make align run=6
align:
	/opt/ont-guppy/bin/minimap2 -t 8 -x map-ont -a ~/genomics/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ~/genomics/basecall-latest/bioliquid_$(run_number).fastq > ~/genomics/bioliquid_$(run_number).sam
	samtools view -@ 32 -bSh ~/genomics/bioliquid_$(run_number).sam          > ~/genomics/bioliquid_$(run_number)_unsorted.bam
	samtools sort -@ 32      ~/genomics/bioliquid_$(run_number)_unsorted.bam > ~/genomics/bioliquid_$(run_number).bam
	samtools index -@ 32     ~/genomics/bioliquid_$(run_number).bam
	samtools flagstat        ~/genomics/bioliquid_$(run_number).bam          > ~/genomics/bioliquid_$(run_number).bam.flag

# QUALIMAP:
# $make qc run=6
qc:
	docker run --rm -v /home/fer/genomics:/data pegi3s/qualimap qualimap bamqc -bam /data/bioliquid_$(run_number).bam
	s3cmd put -r ~/genomics/bioliquid_$(run_number)_stats s3://aretian-genomics/nanopore/$(run_number)/
	
get_ref: $(get_ref) 
$(get_ref): 0_get_reference.sh
	@bash 0_get_reference.sh $(run) $(dis)

extract_ref: $(extract_ref) 
$(extract_ref): $(get_ref) 0_extract_reference.sh
	@bash 0_extract_reference.sh $(run) $(dis)

extract_reads: $(extract_reads) 
$(extract_reads): $(extract_ref) 03_extract_reads.sh
	@bash 03_extract_reads.sh $(run) $(dis)

clean_extract_reads: 
	rm $(extract_reads)
    
remove_gaps: $(remove_gaps)
$(remove_gaps): $(extract_ref) $(extract_reads) 04_remove_gaps.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 04_remove_gaps.py $(run) $(dis)
    
clean_remove_gaps: 
	@rm $(remove_gaps)

cluster: $(cluster)
$(cluster): $(remove_gaps) 05_clustering.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 05_clustering.py $(run) $(dis)

clean_cluster:
	@rm $(cluster)

assign: $(assign)
$(assign): $(cluster) 06_assign_reads.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 06_assign_reads.py $(run) $(dis)

create_bams: $(create_bams)
$(create_bams): $(assign) 07_create_bams.sh
	@bash 07_create_bams.sh $(run) $(dis)

strspy_config: $(strspy_config)
$(strspy_config): $(create_bams) 08_strspy_config.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 08_strspy_config.py $(run) $(dis)

strspy: $(strspy)
$(strspy): $(strspy_config) 09_strspy.sh
	@bash 09_strspy.sh $(run) $(dis)

clean_strspy:
	rm -rf $(datadir)/strspy/$(dis)/output

str_list: $(str_list)
$(str_list): 0_str_list.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 0_str_list.py

add_strs: $(add_strs)
$(add_strs): $(str_list) $(strspy) 10_add_strs.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 10_add_strs.py $(run) $(dis)

tag_reads: $(tag_reads)
$(tag_reads): $(add_strs) 11_tag_reads.R
	@/usr/bin/Rscript 11_tag_reads.R $(run) $(dis)

long_reads: $(long_reads)
$(long_reads): $(tag_reads) 12_long_reads.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 12_long_reads.py $(run) $(dis)

boolean_matrix: $(boolean_matrix)
$(boolean_matrix): $(long_reads) 13_boolean_matrix.R
	@/usr/bin/Rscript 13_boolean_matrix.R $(run) $(dis)

str_cluster: $(str_cluster)
$(str_cluster): $(boolean_matrix) 14_str_clustering.R
	@/usr/bin/Rscript 14_str_clustering.R $(run) $(dis)
	@mkdir -p /home/fer/genomics/bioliquid-nanopore/cluster_plots
	@cp $(datadir)/$(run_number)_$(chrom_dis)_assigned_kmeans_clusters.png /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@cp $(datadir)/$(run_number)_$(chrom_dis)_real_sample_labels.png       /home/fer/genomics/bioliquid-nanopore/cluster_plots/
	@echo Copying plots to home folder.
	@mkdir -p /home/fer/genomics/bioliquid-nanopore/cluster_plots
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/$(run_number)_$(chrom_dis)_assigned_kmeans_clusters.png
	@echo Saved: /home/fer/genomics/bioliquid-nanopore/cluster_plots/$(run_number)_$(chrom_dis)_real_sample_labels.png 

score: $(score)
$(score): $(str_cluster) 15_score.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 15_score.py $(run) $(dis)
    
sample_id: $(sample_id)
$(sample_id): $(score) 16_sample_id.py
	@/home/fer/miniconda3/envs/genomics/bin/python3 16_sample_id.py $(run) $(dis)

disease_diag: $(sample_id) 17_disease_diagnostic.R
	@/usr/bin/Rscript 17_disease_diagnostic.R $(run) $(dis)
