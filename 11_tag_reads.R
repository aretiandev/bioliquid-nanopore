# 11 - Tag Reads
#
# This script links each read for person0 and person1 with the presence or absence of STRs. The script is based on 08_tag_reads.ipynb
#
# INPUTS:
#   run_number
#   disease
#   Identified STRs: e.g. {datadir}/strspy/${dis}/output/Countings
#   nanopore reads for location: e.g. run1_chr11_sca.sam
# 
# OUTPUTS:
#   reads without gaps: e.g. run1_chr11_sca_clean.csv
#   reference genome without gaps: e.g. run1_chr11_reference_genome.json
print('')
print('11 - TAG READS')
# Set variables
args = commandArgs(trailingOnly=TRUE)
run_num = args[1]
dis = args[2]
print(paste0("Run: ", run_num, ", disease: ", dis))
print('')
# Load Modules
# -----------------------------------------------------------------------------
print('Loading modules.')
suppressMessages(library(Rsamtools, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
# suppressMessages(library(ggplot2, quietly=TRUE))

# Set Parameters
# -----------------------------------------------------------------------------
# print('Setting parameters.')

run_number=paste0("run",run_num)

dis_params = data.frame('disease' = c('sca', 'cystic', 'spinal1', 'spinal2', 'thal1', 'thal2', 'thal3', 'pompe'),
                    'chr' = c('chr11','chr7','chr5', 'chr5', 'chr16', 'chr16', 'chr16', 'chr17'),
                    'start' = c(5227002, 117559590, 70924941, 70049523, 176680, 172876, 5225464, 0),
                    'end'   = c(5227002, 117559590, 70953015, 70077595, 177522, 173710, 5227071, 0))

chr = dis_params %>% filter(disease == dis) %>% select(chr) %>% pull

datadir = paste0('/mnt/aretian/genomics/nanopore/',run_number)

P0bam <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_', chr,'_','person0.bam')
P0ind <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_', chr,'_','person0.bam.bai')
P1bam <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_', chr,'_','person1.bam')
P1ind <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_', chr,'_','person1.bam.bai')

Person0Bam <- BamFile(file = P0bam, index = P0ind)
Person1Bam <- BamFile(file = P1bam, index = P1ind)
samples <- c(Person0Bam, Person1Bam)

# Load data
# -----------------------------------------------------------------------------
print(paste0('Loading STRspy output: ',datadir,'/',run_number,'_',chr,'_person_full.txt'))

vcf_filepath <- paste0(datadir,'/',run_number,'_',chr,'_person_full.txt')
vcf <- read.table(vcf_filepath, sep = '\t', stringsAsFactors = FALSE)
# x <- readLines(vcf_filepath)

# Set column names
names <- c('name', 'count', 'chr', 'start', 'end', 'motif','str')
colnames(vcf) <- names
vcf$strname <- paste0('str', seq(1:nrow(vcf)))

# Helper functions
# -----------------------------------------------------------------------------
get_samplename <- function(sample){
    samplename <- unlist(strsplit(sample$`.->path`, '/'))[grep('.bam', unlist(strsplit(sample$`.->path`, '/')))]
    samplename <- unlist(strsplit(samplename, '\\.'))[1] # changed to reflect current bam paths
    return(samplename)
}

# init_reads: This function creates a dataframe with all reads that overlap a given str
# inputs: str: str row number in the vcf, sample, vcf
# output: dataframe with reads attached to a given str
init_reads <- function(str, sample, vcf) {
    
    # get sample name
    samplename <- get_samplename(sample)
    
    str_chr <- vcf[str, 'chr']
    str_start <- vcf[str, 'start']
    str_end <- vcf[str, 'end']
    str_id <- vcf[str, 'name']
    str_pattern <- vcf[str, 'str']
    motif <- vcf[str, 'motif']

    # isolate reads covering that area 
    gr <- GRanges(seqnames = chr,
                  ranges = IRanges(start = str_start, end = str_end))
    params <- ScanBamParam(which = gr, what = scanBamWhat())
    aln <- scanBam(sample, param = params)
    
    # handle cases with 0 reads
    if (length(aln[[1]]$pos) == 0){
        reads <- NA
    } else {
        reads <- data.frame('read_id'=aln[[1]]$qname,
                            'str_id' = str_id,
                            'samplename' = samplename,
                            'startpos' = aln[[1]]$pos, 
                            'str_start'= 0,
                            'str_end' = 0,
                            'read' = aln[[1]]$seq,
                            'str' = str_pattern,
                            'motif' = motif,
                            'has_str' = 0)

        reads$str_start <- abs(reads$startpos - str_start)
        reads$str_end <- abs(reads$startpos - str_end)
    }
    
    return(reads)
}

# update reads function updates reads to include str sequences and read categorizations
update_reads <- function(reads, vcf, str) {

    # for each row (read) in reads
    for (read in c(1:nrow(reads))) {
            
        # Find STR in the reads
        if (grepl(reads[read,'motif'], reads[read,'read'])) {
            reads[read,'has_str'] <- 1
        }
    }

    return(reads)
}

# Implementation
# -----------------------------------------------------------------------------
print(paste0('Tagging reads. Saving: ',datadir,'/strs'))
# for each STR
for (str in c(1:nrow(vcf))) { #
                                 
    data <- as_tibble(data.frame('read_id'='',
                                 'str_id' = '',
                                 'samplename' = '',
                                 'startpos' = 0, 
                                 'str_start'= 0,
                                 'str_end' = 0,
                                 'read' = '',
                                 'str' = '',
                                 'motif' = '',
                                 'has_str' = 0))
    data <- data[-1,]
    
    # for each sample
    for (sample in samples) { # samples = c(Person0Bam, Person1Bam)

        # Get all reads that cover the VCF
        reads <- init_reads(str = str, sample = sample, vcf = vcf)
        if (is.null(nrow(reads))) {
#             meets_threshold <- FALSE
            break
        }
        
        # Get reads that contain the STR pattern
        reads <- update_reads(reads = reads, vcf = vcf, str = str)
        
        # remove reads without motif
        reads <- reads[reads$has_str == 1,]
        
        # update summary dataframe
        data <- rbind(data, reads)
    }

#     data <- data %>% mutate(uid = paste0(sample, '-', readgroup))
    data <- data %>% mutate(read_id = paste0(samplename, '-', read_id))
    data <- data %>% select(read_id, str_id, startpos)

    
    # save each str df (only if all samples met the threshold)
#     if (meets_threshold) {
#        write.csv(data, paste0('/mnt/aretian/genomics/data/str_pipeline/output/', disease, '-test/', fid, '/str', str, '.csv'), row.names = FALSE)
#     }
    write.csv(data, paste0(datadir,'/strs/str', str, '.csv'), row.names = FALSE)

    # time / str tracker:
    cat('\r',paste0('Progress: STR ', str, ' out of ',nrow(vcf),'.'))
#     print(paste0('str', str, ': ', Sys.time()))
}

# Concatenate all results
# -----------------------------------------------------------------------------
print("")
print(paste0('Concatenating and removing all individual files.'))
input = paste0(datadir,'/strs/str[[:digit:]]*.csv')
output = paste0(datadir,'/',run_number,'_',chr,'_tagged_reads.csv')

concatenate_cmd = paste0("awk 'FNR==1 && NR!=1{next;}{print}' ", input, ' > ', output)
system(concatenate_cmd)

system(paste0('rm ', input))

print(paste0('Saved: ', output))
