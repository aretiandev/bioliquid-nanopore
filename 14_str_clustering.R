# 14 - STR Clustering
#
# This script performs clustering of the person-reads based on the boolean matrix.
#
# INPUTS:
#   run_number
#   disease
#   boolean matrix of tagged reads: e.g. datadir/run_chr_bool_tagged_reads.csv'
# 
# OUTPUTS:
#   cluster results: datadir/run_chr_kmeans_clusters.csv'
args = commandArgs(trailingOnly=TRUE)
print('')
print('------------------------------------------------------------------------------------------')
print('14 - STR CLUSTERING (14_str_clustering.R)')
print(paste0("Run: ", args[1], ", disease: ", args[2]))
print('')
print('')

# Set variables
run_num = args[1]
dis = args[2]

# Load Modules
# -----------------------------------------------------------------------------
print('Loading Modules.')
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(factoextra))
suppressMessages(library(FactoMineR))
suppressMessages(library(dplyr))

# Set Parameters
# -----------------------------------------------------------------------------
# print('Setting parameters.')

run_number=paste0("run",run_num)

dis_params = data.frame('disease' = c('sca', 'cf', 'sma1', 'sma2', 'thal1', 'thal2', 'thal3', 'pompe'),
                    'chr' = c('chr11','chr7','chr5', 'chr5', 'chr16', 'chr16', 'chr11', 'chr17'),
                    'start' = c(5227002, 117559590, 70924941, 70049523, 176680, 172876, 5225464, 25000000),
                    'end'   = c(5227002, 117559590, 70953015, 70077595, 177522, 173710, 5227071, 25000000))

chr = dis_params %>% filter(disease == dis) %>% select(chr) %>% pull

datadir = paste0('/mnt/aretian/genomics/nanopore/',run_number)

# Helper functions
# -----------------------------------------------------------------------------
mysplit <- function(x) {
  y <- unlist(strsplit(x, '-'))[1]
  return(y)
}

constfunc <- function(x) {
  if (length(unique(x)) > 1) {return(1)}
  else {return(0)}
}

# Load data
# -----------------------------------------------------------------------------
print(paste0('Loading Boolean Matrix:                 ',datadir,'/',run_number,'_',chr,'_',dis,'_bool_tagged_reads.csv'))

data <- read.csv(paste0(datadir,'/',run_number,'_',chr,'_',dis,'_bool_tagged_reads.csv'), stringsAsFactors = FALSE)
data$sample <- unlist(lapply(data$read_id, mysplit))

data$sample <- as.factor(data$sample)
data_nolabels <- data %>% select(-read_id, -sample)

# remove constant columns
presents_variation <- apply(X = data_nolabels, MARGIN = 2, FUN = constfunc)
presents_variation <- as.logical(presents_variation)

data_nolabels_noconstants <- data_nolabels[,presents_variation]

# Implementation
# -----------------------------------------------------------------------------
#################################################################################
# K-means Clustering
#################################################################################
set.seed(123)

kres <- kmeans(data_nolabels_noconstants, centers = 2)


print(paste0('Generating plot with assigned clusters: ',datadir,'/',run_number,'_',chr,'_',dis,'_assigned_kmeans_clusters.png'))

png(file=paste0(datadir,'/',run_number,'_',chr,'_',dis,'_assigned_kmeans_clusters.png'))
fviz_cluster(kres, data = data_nolabels_noconstants, labelsize = 0,
             ellipse = TRUE, ellipse.type = "convex",
              ellipse.level = 0.95, ellipse.alpha = 0.2,
             main = 'K-Means clustering results: Assigned clusters',
             submain = 'Person0-Person1 merged samples',
#              caption = paste0(mid, ' = mother, ', cid, ' = child'),
             xlab = 'PC1',
             ylab = 'PC2',
             palette = 'Set2') + 
  theme_bw()
dev.off()

kres_lab <- kres # color by real group
kres_lab$cluster <- data$sample

print(paste0('Generating plot with real clusters:     ',datadir,'/',run_number,'_',chr,'_',dis,'_real_sample_labels.png'))

png(file=paste0(datadir,'/',run_number,'_',chr,'_',dis,'_real_sample_labels.png'))
fviz_cluster(kres_lab, data = data_nolabels_noconstants, labelsize = 0,
             ellipse = TRUE, ellipse.type = "convex",
              ellipse.level = 0.95, ellipse.alpha = 0.2,
             main = 'K-Means clustering results: Actual sample labels', 
             submain = 'Person0-Person1 merged samples',
#              caption = paste0(mid, ' = mother, ', cid, ' = child'),
             xlab = 'PC1',
             ylab = 'PC2',
             palette = 'Set1') + 
  theme_bw()
dev.off()

data$kmeans_clusters <- kres$cluster
write.csv(data, paste0(datadir,'/',run_number,'_',chr,'_',dis,'_kmeans_clusters.csv'), row.names = FALSE)
print(paste0('Saved:                                  ',datadir,'/',run_number,'_',chr,'_',dis,'_kmeans_clusters.csv'))
