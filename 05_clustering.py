# 05 - Clustering
#
# This script performs the clustering algorithm. It contains the following steps:
#
# - Padding reads with the reference genome to fill empty spaces
# - Sliding a window over the region of interest to identify reads
#
# The script is based on 05_clustering.ipynb
#
# INPUTS:
#   run_number
#   disease
#   reads without gaps: e.g. run1_chr11_sca_clean.csv
#   reference genome without gaps: e.g. run1_chr11_reference_genome.json
# 
# OUTPUTS:
#   Clustered reads with windows: run1_chr11_read_clusters.txt
print('')
print('05 - CLUSTERING')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import json
import os
import sys
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
encoder = OneHotEncoder()

# Set Variables
# -----------------------------------------------------------------------------
# Input variables
run_num = sys.argv[1]
dis = sys.argv[2]
print(f"Run: {run_num}, disease: {dis}.")

# Setup
run_number=f"run{run_num}"

if dis =="cf":
    chrom = "chr7"
elif dis == "sca":
    chrom = "chr11"
elif dis == "sma":
    chrom = "chr5"
elif dis == "thal":
    chrom = "chr16"
elif dis == "pompe":
    chrom = "chr17"
else:
    print("Disease should be in disease list: cf, sca, sma, thal, pompe.")

chrom_dis=f"{chrom}_{dis}"
datadir=f"/mnt/aretian/genomics/nanopore/{run_number}"
results_file_path = f"{datadir}/{run_number}_{chrom}_read_clusters.txt"

# Import data
# -----------------------------------------------------------------------------
# Import Nanopore reads
print(f'Reading clean Nanopore reads from: {datadir}/{run_number}_{chrom_dis}_clean.csv')
nanopore_reads_path = f'{datadir}/{run_number}_{chrom_dis}_clean.csv'
nanopore_reads = pd.read_csv(nanopore_reads_path)

# Import reference genome for the location of interest without gaps
print(f'Reading Reference Genome from:     {datadir}/{run_number}_{chrom}_reference_genome.json')
with open(f'{datadir}/{run_number}_{chrom}_reference_genome.json', 'r') as f:
    ref_genome_json = json.load(f)
    
ref_genome = ref_genome_json['reference_genome']

# Helper functions
# -----------------------------------------------------------------------------
# Define hyper-parameters
window_width = 5000 # 10K BP
jump_width = 1000 # 1k BP sliding window (ideally jump_width=1)

# check if the start of a sequence is in range
def check_sliding_window(read):
    read_filter = False
    
# #     Slow algorithm
#     window = set(list(range(left_bound,right_bound+1)))
#     read = set(list(range(read['POS'],read['END_POS']+1)))
#     if len(window & read)>0:
#         read_filter = True

    starts_inside = read['POS']>=left_bound and read['POS']<=right_bound
    ends_inside = read['END_POS']>=left_bound and read['END_POS']<=right_bound
    covers_window = read['POS']<left_bound and read['END_POS']>right_bound
    if starts_inside or ends_inside or covers_window:
        read_filter=True
    
    return read_filter

# Input: read --> Series
# Output: nucleotide_sequence_trimmed: String

def left_pad_read(read):
    nucleotide_sequence = list(read['SEQ']) # turn sequence string into a list
    genome_sequence = list(window_ref_genome)
    
    start = read['POS']-left_bound

    if start == 0:
        return read['SEQ'] # read starts on the start of the window
    
    elif start > 0: # sequence starts to the right of the bound
#         print('HERE1')
#         genome_fill = window_ref_genome[:start] # select portion of reference genome to pad the left with
        genome_fill = genome_sequence[:start] # select portion of reference genome to pad the left with
#         print('genome_fill type:', type(genome_fill))
#         print('nucleotide_sequence type:', type(nucleotide_sequence))
        nucleotide_sequence = genome_fill+nucleotide_sequence
        nucleotide_sequence_trimmed = ''.join(nucleotide_sequence) # keep entire sequence
        return nucleotide_sequence_trimmed
    
    elif start < 0: # sequence starts to the left of the bound
#         print('HERE2')
        nucleotide_sequence_trimmed = ''.join(nucleotide_sequence[np.abs(start):])
        return nucleotide_sequence_trimmed

# Input: read --> Series
# Output: nucleotide_sequence_trimmed: String

def right_pad_read(read):
    nucleotide_sequence = list(read['left_padded']) # turn sequence string into a list
    genome_sequence = list(window_ref_genome)

    end = right_bound-read['END_POS']

    if end == 0:
        return read['left_padded'] # read end on the end of the window
    
    elif end > 0: # sequence ends to the left of the bound
#         genome_fill = window_ref_genome[-end:] # select portion of reference genome to pad the right with
        genome_fill = genome_sequence[-end:] # select portion of reference genome to pad the right with
        nucleotide_sequence = nucleotide_sequence+genome_fill
        nucleotide_sequence_trimmed = ''.join(nucleotide_sequence) # keep entire sequence
        return nucleotide_sequence_trimmed
    
    elif end < 0: # sequence ends to the right of the bound
        nucleotide_sequence_trimmed = ''.join(nucleotide_sequence[:end])
        return nucleotide_sequence_trimmed

# Should be updated
NUCLEOTIDE_VOCABULARY = [
    'A','C','G','T','X'
]
        
# Not being used, instead, we use the sklearn one hot encoding
def nucleotide_to_one_hot(nucleotide_sequence):
    to_return = []
    for char in nucleotide_sequence:
        if char in NUCLEOTIDE_VOCABULARY:
            to_append = np.zeros(len(NUCLEOTIDE_VOCABULARY))
            to_append[NUCLEOTIDE_VOCABULARY.index(char)] = 1.
            to_return.append(to_append)
        else:
            raise ValueError('Could not one-hot code character {}'.format(char))
    return np.array(to_return)

#nucleotide_to_one_hot('GTCATACX') # uncomment example to see what the encoding does

# Input: read (Series)

# Text file for results
def write_results_to_file(read, results_file):
    results_file.write(f"{read['ID']},{read['kmeans_cls2']},{read['window_num']}\n")

# Main Loop: padding and clustering
# -----------------------------------------------------------------------------
print("")
print('Performing padding and clustering.')
# User feedback
print(f"Total range: {max(nanopore_reads['END_POS']):,}")
print(f"Window width: {window_width}")
print(f"Jump width: {jump_width}")
total_jumps = round(max(nanopore_reads['END_POS'])/jump_width)
print(f"Number of iterations required: {total_jumps}")
print(f"Reading Nanopore reads from: {nanopore_reads_path}")
print(f"Writing results to file:     {results_file_path}")
print("")
print("0: Found 0 reads in window. Saving empty window.")
print("1: Found 1 read in window. Applying arbitrary cluster 0.")
print("*: Found more than 1 read in window and ran clustering. Saving clusters.")
print("")
print(f"Iterations: ", end="")

results_file = open(results_file_path,"w")
iter = 1
empty_count = 0
threshold = 5
for left_bound in range(min(nanopore_reads['POS']),max(nanopore_reads['END_POS']),jump_width):
# for left_bound in range(950000,970000,jump_width):

#     print(f" {progress:.0f}% ", end="")

    # Print progress in jumps of 5%
    progress = iter/total_jumps*100
    if progress>=threshold:
        print(f'{threshold}%')
        threshold = threshold+5

    iter+=1
#     if round(progress*1000)%50==0:
#         print(f" {round(progress*100):.0f}% ", end="", flush=True)
        
    right_bound = left_bound+window_width
    window_ref_genome = ref_genome[left_bound:left_bound+window_width]
    
    # Identify if each read is in the window: True/False
#     print("Identifying if each read is in the window...")
    nanopore_reads['read_filter'] = nanopore_reads.apply(lambda x: check_sliding_window(x), axis=1)
    
    # Get them
    window_reads = nanopore_reads.loc[nanopore_reads['read_filter']==True] 
#     print(f" {window_reads['ID'].unique()}-", end="")
    window_reads = window_reads.reset_index()
    
    # Display feedback to user every 100 empty windows
#     if empty_count%100 == 0:
#         print(empty_count)

    if len(window_reads) == 0:
        # There are no nanopore reads in this window (next window)
        empty_count += 1
        print(f"0", end="", flush=True)
        
        # Record window with empty read
        data = {'ID':[np.nan],'kmeans_cls2':[np.nan], 'window_num':[iter]}
        empty_window_reads = pd.DataFrame(data)
   
        empty_window_reads.apply(lambda x: write_results_to_file(x, results_file), axis=1) # write results to file
        continue 
        
    elif len(window_reads) == 1:
        # There is only one read in this window
        empty_count += 1
        print(f"1", end="", flush=True)
        
        # Record window and assign arbitrary individual
        window_reads['window_num']=iter
        window_reads['kmeans_cls2']=0
        window_reads.apply(lambda x: write_results_to_file(x, results_file), axis=1) # write results to file
        
        continue
        
    # More than one read in window. Proceeding to padding and clustering

#     Padding
    window_reads['left_padded'] = window_reads.apply(lambda x: left_pad_read(x), axis=1) # fill reference genome on the left of the read
    window_reads['right_padded'] = window_reads.apply(lambda x: right_pad_read(x), axis=1) # fill reference genome on the right of the read
    window_reads['final_padded_read'] = window_reads['right_padded']
    window_reads['FINAL_SEQ_LEN'] = window_reads['final_padded_read'].apply(lambda x: len(x)) # should always be 5k
    
    # not currently used. Should be used if we use the custom one-hot-encoding
#     window_reads['one_hot_read_V1'] = window_reads['final_padded_read'].apply(lambda x: nucleotide_to_one_hot(x).flatten())  # apply one-hot encoding V1

#     One hot encoding
    unique_reads = []
    for index, read in window_reads.iterrows(): # TODO: TRY TO USE LAMBDA FUNCTION IF POSSIBLE
        unique_reads.append(list(read['final_padded_read']))
    X_onehot = encoder.fit_transform(unique_reads).toarray()

    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_onehot)
    
    window_reads['PCA1'] = np.nan
    window_reads['PCA2'] = np.nan
    for ID in window_reads.index:
        window_reads.loc[ID,'PCA1'] = X_pca[ID][0] # TODO: change rounding if desired
        window_reads.loc[ID,'PCA2'] = X_pca[ID][1] # TODO: change rounding if desired
        
    # Standardizing the features
    X = window_reads[['PCA1','PCA2']]
    X = StandardScaler().fit_transform(X)

    # Run Kmeans
    model = KMeans(n_clusters=2, random_state=42)
    cls2 = model.fit(X)
    window_reads['kmeans_cls2'] = cls2.labels_
    window_reads['window_num'] = iter
    
    # Write results to file
    window_reads.apply(lambda x: write_results_to_file(x, results_file), axis=1) # write results to file
    
    print(f"*", end="", flush=True)
    
#     Break after N iterations
#     if iter>100:
#         break
        
results_file.close()
print("")
print("Done padding and clustering.")
print(f"Saved results: {datadir}/{run_number}_{chrom}_read_clusters.txt")
