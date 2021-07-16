# 11 - Score
#
# This script evaluates clustering accuracy.. The script is based on 11_evaluate_clusters.ipynb
#
# INPUTS:
#   run_number
#   disease
#   cluster results: datadir/run_chr_kmeans_clusters.csv'
# 
# OUTPUTS:
#   recall scores: e.g. datadir/run_chr_recall_score.csv
print('')
print('11 - SCORE')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import pandas as pd
import os
import sys
from sklearn.metrics import precision_score, recall_score

# Set variables
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
    location=5227002
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

# Load data
# -----------------------------------------------------------------------------
print(f'Reading clustering results: {datadir}/{run_number}_{chrom}_kmeans_clusters.csv')
df = pd.read_csv(f'{datadir}/{run_number}_{chrom}_kmeans_clusters.csv')
df = df[['sample','kmeans_clusters']]
df.columns = ['true','predicted']

# Clean and manipulate data
df['true'] = df['true'].str[-1:].astype(int)
df['predicted'] = df['predicted']-1

# Calculate Recall Score
# -----------------------------------------------------------------------------
print(f'Saving score results:       {datadir}/{run_number}_{chrom}_recall_score.csv')

recall0 = recall_score(df['true'],df['predicted'], pos_label=0)
print(f"The recall score for cluster 0 is {recall0}")

recall1 = recall_score(df['true'],df['predicted'], pos_label=1)
print(f"The recall score for cluster 1 is {recall1}")

# Save results
with open(f'{datadir}/{run_number}_{chrom}_recall_score.csv', 'w') as f:
    f.write(f"{recall0},{recall1}")
