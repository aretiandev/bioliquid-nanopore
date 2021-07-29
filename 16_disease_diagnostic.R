# 16 - Disease Diagnostic
#
# This script performs disease diagnostic. Reads bam files of two persons and extracts relevant reads to analyze.
#
# INPUTS:
#   run_number
#   disease
#   BAM files
# 
# OUTPUTS:
#   relevant reads with the SNP
args = commandArgs(trailingOnly=TRUE)
print('')
print('------------------------------------------------------------------------------------------')
print('16 - DISEASE DIAGNOSTIC (16_disease_diagnostic.R)')
print(paste0("Run: ", args[1], ", disease: ", args[2]))
print('')
print('')

# Set variables
# -----------------------------------------------------------------------------
run_num = args[1]
dis = args[2]

# Load Modules
# -----------------------------------------------------------------------------
print('Loading modules.')
suppressMessages(library(tidyverse, quietly=TRUE))
suppressMessages(library(Rsamtools, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))

# Set Parameters
# -----------------------------------------------------------------------------
# print('Setting parameters.')

run_number=paste0("run",run_num)

dis_params = data.frame('disease' = c('sca', 'cf', 'sma1', 'sma2', 'thal1', 'thal2', 'thal3', 'pompe'),
                    'chr' = c('chr11','chr7','chr5', 'chr5', 'chr16', 'chr16', 'chr11', 'chr17'),
                    'start' = c(5227002, 117559590, 70924941, 70049523, 176680, 172876, 5225464, 25000000),
                    'end'   = c(5227002, 117559590, 70953015, 70077595, 177522, 173710, 5227071, 25000000))

chr       = dis_params %>% filter(disease == dis) %>% select(chr)   %>% pull
snp_start = dis_params %>% filter(disease == dis) %>% select(start) %>% pull
snp_end   = dis_params %>% filter(disease == dis) %>% select(end)   %>% pull

datadir = paste0('/mnt/aretian/genomics/nanopore/',run_number)

# P0bam <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_','person0.bam')
# P0ind <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_','person0.bam.bai')
# P1bam <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_','person1.bam')
# P1ind <- paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_','person1.bam.bai')

# Person0Bam <- BamFile(file = P0bam, index = P0ind)
# Person1Bam <- BamFile(file = P1bam, index = P1ind)
# samples <- c(Person0Bam, Person1Bam)

# Load data
# -----------------------------------------------------------------------------
print('Loading data.')

person0bam <- BamFile(
                file  = paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_person0.bam'    ), 
                index = paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_person0.bam.bai'))

person1bam <- BamFile(
                file  = paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_person1.bam'    ), 
                index = paste0(datadir,'/strspy/',dis,'/input/',run_number,'_',chr,'_',dis,'_person1.bam.bai'))


# Diagnostic and save table
# -----------------------------------------------------------------------------
print('Extracting reads containing the relevant SNP.')

gr <- GRanges(seqnames = chr,
              ranges = IRanges(start = snp_start, end = snp_end))
params <- ScanBamParam(which = gr, what = scanBamWhat())

# Person 0
aln <- scanBam(person0bam, param = params)
reads <- data.frame('startpos' = aln[[1]]$pos, 
                    'read' = aln[[1]]$seq)
reads2 <- reads                                                        %>%
    mutate(length         = nchar(read)                   ,
           snp_start_0    = snp_start - startpos          ,
           snp_end_0      = snp_end   - startpos          ,
           left_padding   = str_sub(read,0,snp_start_0   ),
           right_padding  = str_sub(read,  snp_end_0  +2 ),
           segment        = str_sub(read,  snp_start_0+1  , snp_end_0+1),
           snp_start      = snp_start                     ,
           snp_end        = snp_end                       )            %>%
    select(snp_start, snp_end, left_padding,  segment, right_padding)  %>%
    sample_n(20)

write.csv(reads2, paste0(datadir,'/',run_number,'_',chr,'_',dis,'person0_diagnostic_reads.csv'), row.names = FALSE)
write.csv(reads2, paste0('disease_diagnostic/',run_number,'_',chr,'_',dis,'person0_diagnostic_reads.csv'), row.names = FALSE)

# Person 1
aln <- scanBam(person1bam, param = params)
reads <- data.frame('startpos' = aln[[1]]$pos, 
                    'read' = aln[[1]]$seq)

reads2 <- reads                                                        %>%
    mutate(length         = nchar(read)                   ,
           snp_start_0    = snp_start - startpos          ,
           snp_end_0      = snp_end   - startpos          ,
           left_padding   = str_sub(read,0,snp_start_0   ),
           right_padding  = str_sub(read,  snp_end_0  +2 ),
           segment        = str_sub(read,  snp_start_0+1  , snp_end_0+1),
           snp_start      = snp_start                     ,
           snp_end        = snp_end                       )            %>%
    select(snp_start, snp_end, left_padding,  segment, right_padding)  %>%
    sample_n(20)

write.csv(reads2, paste0(datadir,'/',run_number,'_',chr,'_',dis,'person1_diagnostic_reads.csv'), row.names = FALSE)
write.csv(reads2, paste0('disease_diagnostic/',run_number,'_',chr,'_',dis,'person1_diagnostic_reads.csv'), row.names = FALSE)

print(paste0('Saved: ',datadir,'/',run_number,'_',chr,'_',dis,'person0_diagnostic_reads.csv'))
print(paste0('Saved: ',datadir,'/',run_number,'_',chr,'_',dis,'person1_diagnostic_reads.csv'))
print('Saved copies of each csv in disease_diagnostic/ folder in homedir.')