# 16 - Sample ID
#
# This script performs sample identification.
#
# INPUTS:
#   run_number
#   disease
#   bool_tagged_reads.csv
# 
# OUTPUTS:
#   read numbers and shares and bool_matrix plot
import sys
print('')
print('------------------------------------------------------------------------------------------')
print(f'16 - SAMPLE ID ({__file__})')
print(f"Run: {sys.argv[1]}, disease: {sys.argv[2]}.")
print('')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import pandas as pd
import os
import matplotlib.pyplot as plt

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
df = pd.read_csv(f'{datadir}/{run_number}_{chrom_dis}_bool_tagged_reads.csv')

# Calculate read and STR shares
# -----------------------------------------------------------------------------
# Read shares
print("Calculating mother-child shares and STR shares.")
# Create Dataframe with shares
df['person']=df['read_id'].str.split('_').str[3].str[:7]
df_out = pd.DataFrame({'n_long_reads':df['person'].value_counts()})
df_out['share_long_reads']=df_out['n_long_reads']/df_out['n_long_reads'].sum()
df_out.index = df_out.index.rename('person')
# Save
df_out.to_csv(f'{datadir}/{run_number}_{chrom_dis}_sample_id_read_shares.csv')
df_out.to_csv(f'disease_diagnostic/{run_number}_{chrom_dis}_sample_id_read_shares.csv')
print(f"Saved: {datadir}/{run_number}_{chrom_dis}_sample_id_read_shares.csv")

# STR shares
sums_df = df.groupby('person').sum().sum(axis=1)
sums_df = pd.DataFrame({'n_strs':sums_df})
sums_df['share_strs'] = sums_df['n_strs']/sums_df['n_strs'].sum()
sums_df['is_child']=0
sums_df.loc[sums_df['n_strs']<sums_df['n_strs'].mean(), 'is_child']=1
# Save
sums_df.to_csv(f'{datadir}/{run_number}_{chrom_dis}_sample_id_str_shares.csv')
sums_df.to_csv(f'disease_diagnostic/{run_number}_{chrom_dis}_sample_id_str_shares.csv')
print(f"Saved: {datadir}/{run_number}_{chrom_dis}_sample_id_str_shares.csv")
print("Created a copy in the homedir.")

# Plot boolean matrix
# -----------------------------------------------------------------------------
print('Plotting the boolean matrix.')
df_plot = df.iloc[:,1:-1]
fig, ax = plt.subplots(figsize=(20,5), dpi=100)
myplot = ax.imshow(df_plot.to_numpy(),
          cmap=plt.cm.Blues, 
          aspect='auto',
          interpolation='none')

fig.colorbar(myplot, ax=ax)
ax.set_xlabel('STRs')
ax.set_ylabel('Long Reads')

# Save fig
# -----------------------------------------------------------------------------
plt.savefig(f"cluster_plots/{run_number}_{chrom_dis}_boolean_matrix.png",dpi=100, bbox_inches='tight')
print(f"Saved: cluster_plots/{run_number}_{chrom_dis}_boolean_matrix.png")