# 15 - Score
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
import sys
print('')
print('----------------------------------------------------------------------')
print(f'15 - SCORE ({__file__})')
print(f"Run: {sys.argv[1]}, disease: {sys.argv[2]}.")
print('')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import pandas as pd
import os
from sklearn.metrics import precision_score, recall_score

# Set Variables
# ---------------------------------------------------------------------------------------------------
run_num = sys.argv[1]
dis = sys.argv[2]

from src.setup_variables import *
try:
    chrom=dis_data[dis]['chr']
    location=dis_data[dis]['location']
except:
    print("Disease should be in disease list: cf, sca, sma1, sma2, thal1, thal2, thal3, pompe.")
    
run_number=f"run{run_num}"
chrom_dis=f"{chrom}_{dis}"
datadir=f"{rootdir}/{run_number}"

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
recall0 = recall_score(df['true'],df['predicted'], pos_label=0)
print(f"The recall score for cluster 0 is {recall0}")

recall1 = recall_score(df['true'],df['predicted'], pos_label=1)
print(f"The recall score for cluster 1 is {recall1}")

# Save results
with open(f'{datadir}/{run_number}_{chrom}_recall_score.csv', 'w') as f:
    f.write(f"{recall0},{recall1}")
    
print(f'Saved score results:       {datadir}/{run_number}_{chrom}_recall_score.csv')