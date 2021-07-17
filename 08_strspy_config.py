# 08 - STRspy config
#
# This script creates the necessary config files to run STRspy on the Bioliquid Nanopore data. The script is based on 07_STRspy.ipynb
#
# INPUTS:
#   run_number
#   disease
#   reference genome fasta for location: e.g. chr11_selected.fa
#   Full list of STRs
# 
# OUTPUTS:
#   BED and fasta files needed to run STRspy
print('')
print('08 - STRSPY CONFIG')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
import sys

# Set variables
# -----------------------------------------------------------------------------
# Input variables
run_num = sys.argv[1]
dis = sys.argv[2]
windowwidth=2000000
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
    location=25000000
else:
    print("Disease should be in disease list: cf, sca, sma, thal, pompe.")

chrom_dis=f"{chrom}_{dis}"
rootdir=f"/mnt/aretian/genomics/nanopore"
datadir=f"/mnt/aretian/genomics/nanopore/{run_number}"

# Individual BED and Fasta files
# -----------------------------------------------------------------------------
# Get reference genome
# Read in fasta file: remove line breaks and header
print(f'Reading Reference Genome:           {datadir}/{chrom}_selected.fa')

def read_fasta_genome(fasta_file,chromosome_header):
    clean_data = fasta_file.read().replace("\n", "")
    clean_data = clean_data.replace(chromosome_header,"") # get rid of header

    return clean_data

with open(f'{datadir}/{chrom}_selected.fa') as f: # update path if needed
    ref_genome = read_fasta_genome(f,f'>{chrom}')
    
# print(f"Unique characters: {list(set(ref_genome))}") 
# print(f"Selected chromosome from reference genome is {len(ref_genome)} BP long")

# Load list of STRs
# Load Full STR list
print(f'Reading list of STRs:               {rootdir}/hg38.hipstr_reference.bed')
df = pd.read_csv(f'{rootdir}/hg38.hipstr_reference.bed', sep='\t', header=None)
df.columns=['chr','start','end','NA','repeats','name','unit']

# Select Chrom
df = df.loc[df['chr']==chrom]
# Get columns
df = df[['chr','start','end','name']]

# Save each STR in different BED file
selected_strs = df.loc[(df['start']>location-windowwidth)&(df['end']<location+windowwidth)]

# Loop: create single STR files
print(f'Creating BED files for each STR:    {datadir}/strspy/{dis}/input/db' )

if not os.path.exists(f"{datadir}/strspy/{dis}/input/db"):
	os.makedirs(f"{datadir}/strspy/{dis}/input/db")

for n in range(len(selected_strs)):
# for n in range(3):
    str_out = selected_strs.iloc[[n]]
    str_name = str_out['name'].values[0]
    str_out.to_csv(f"{datadir}/strspy/{dis}/input/db/{str_name}.bed", header=False, index=False, sep='\t')
    
    myfasta = open(f"{datadir}/strspy/{dis}/input/db/{str_name}.fa","w")
    start = str_out['start'].values[0]
    end = str_out['end'].values[0]
    # Extract reads
    padded_str=ref_genome[start-500:end+500]
    # Write to file
    myfasta.write('>')
    myfasta.write(str_name)
    myfasta.write('\n')
    myfasta.write(padded_str)
    myfasta.write('\n')
    myfasta.close()

# Region BED file (all STRs)
# -----------------------------------------------------------------------------
if not os.path.exists(f"{datadir}/strspy/{dis}/input/regions"):
	os.makedirs(f"{datadir}/strspy/{dis}/input/regions")
selected_strs.to_csv(f'{datadir}/strspy/{dis}/input/regions/all_strs.bed', header=False, index=False, sep='\t')
print(f'Saved single BED file for all STRs: {datadir}/strspy/{dis}/input/regions/all_strs.bed' )
