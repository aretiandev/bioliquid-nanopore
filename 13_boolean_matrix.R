# 13 - Create Boolean Matrix
#
# This scripts creates a matrix of boolean values per read based on str presence. Example:
# 
#  person-read   str1  str2  str3  str4  str5
#  p0-01         0     1     0     1     0
#  p0-02         1     0     0     1     0
#  p1-01         1     0     0     1     0
#
# INPUTS:
#   run_number
#   disease
#   tagged reads: datadir/run_number_chr_tagged_reads.csv')
# 
# OUTPUTS:
#   boolean matrix of tagged reads: e.g. datadir/run_chr_bool_tagged_reads.csv'
args = commandArgs(trailingOnly=TRUE)
print('')
print('------------------------------------------------------------------------------------------')
print('13 - BOOLEAN MATRIX (13_boolean_matrix.R)')
print(paste0("Run: ", args[1], ", disease: ", args[2]))
print('')
print('')

# Set variables
run_num = args[1]
dis = args[2]

# Load Modules
# -----------------------------------------------------------------------------
print('Loading Modules.')
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))

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

# Load Data
# -----------------------------------------------------------------------------
df <- read.csv(paste0(datadir,'/', run_number,'_',chr,'_long_tagged_reads.csv'), stringsAsFactors = FALSE)

# Helper Functions
# -----------------------------------------------------------------------------
# instead of using length (which refers to length of STR) just use 1 to indicate it exists
agg_func <- function(l){return(1)}

unmeltdf <- function(data) {
    df <- data %>% mutate(val = 1)
    df2 <- dcast(df, read_id ~ str_id, fill = 0, value.var = 'val', fun.aggregate = agg_func)
    return(df2)
}

# Implementation
# -----------------------------------------------------------------------------
print('Creating Boolean Matrix.')
print(paste0('Loaded tagged reads: ',datadir,'/', run_number,'_',chr,'_tagged_reads.csv'))
booldf <- unmeltdf(df)

# Save
# -----------------------------------------------------------------------------
write.csv(booldf, paste0(datadir,'/',run_number,'_',chr,'_bool_tagged_reads.csv'), row.names = FALSE)

print(paste0('Saved:               ', datadir,'/',run_number,'_',chr,'_bool_tagged_reads.csv'))
