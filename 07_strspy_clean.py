# 07 - STRspy clean
#
# This script edits the output of STRspy to the right format with additional columns. The script is based on 07_STRspy.ipynb
#
# INPUTS:
#   run_number
#   disease
#   reference genome fasta for location: e.g. chr11_selected.fa
#   Full list of STRs
# 
# OUTPUTS:
#   VCF file with STRs present: {datadir}/{run_number}_{chrom}_person_full.txt
print('')
print('07 - STRSPY CLEAN')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
import sys
import subprocess

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
else:
    print("Disease should be in disease list: cf, sca, sma, thal, pompe.")

chrom_dis=f"{chrom}_{dis}"
datadir=f"/mnt/aretian/genomics/nanopore/{run_number}"

# Concatenate all output
# -----------------------------------------------------------------------------
print(f'Concatenating all STRs in {datadir}/strspy/{dis}/output/Countings')
concatenate_cmd = f"""
cd "{datadir}/strspy/{dis}/output/Countings"
cat *person0*Allele_freqs.txt > "{run_number}_{chrom}_person0_strs_raw.txt"
cat *person1*Allele_freqs.txt > "{run_number}_{chrom}_person1_strs_raw.txt"

grep Human_STR "{run_number}_{chrom}_person0_strs_raw.txt" > "{run_number}_{chrom}_person0_strs.txt"
grep Human_STR "{run_number}_{chrom}_person1_strs_raw.txt" > "{run_number}_{chrom}_person1_strs.txt"
"""

subprocess.run(concatenate_cmd, shell=True)

# Add necessary columns
# -----------------------------------------------------------------------------
print('Adding motif and full STR columns.')
# This cell takes ~2 minutes to run
def complete_str_df(person):
    # Load STRspy output
    strspy_df = pd.read_csv(f'{datadir}/strspy/output/Countings/{run_number}_{chrom}_{person}_strs.txt', sep='\t', header=None)
    strspy_df.columns = ['name', 'count', 'normcount']

    # Load Full STR list
    df = pd.read_csv(f'{datadir}/hg38.hipstr_reference.bed', sep='\t', header=None)
    df.columns=['chr','start','end','NA','repeats','name','motif']

    ### Create STR
    def create_str(row):
        motif_len = len(row['motif']) # get length
        # Get Base
        int_repeat = int(np.floor(row['repeats'])) # 9
        base = int_repeat * row['motif']
        # Get Tail and append
        dec_repeat = row['repeats']%1
        nt_to_pull = round(dec_repeat * motif_len)
        tail = row['motif'][:nt_to_pull]
        base = base + tail
        return base

    # Drop nans
    df = df.loc[df['motif'].notnull()]
    df['str'] = df.apply(lambda x: create_str(x), axis = 1)

    # Append it to STRspy output
    output = strspy_df.merge(df, how='left', on='name')

    output = output[['name','count','chr','start','end','motif', 'str']]

    # Save
    return output
    
output0 = complete_str_df('person0')
output1 = complete_str_df('person1')

# Save
# -----------------------------------------------------------------------------
output0.to_csv(f'{datadir}/{run_number}_{chrom}_person0_full.txt', index=None, header=None, sep='\t')
output1.to_csv(f'{datadir}/{run_number}_{chrom}_person1_full.txt', index=None, header=None, sep='\t')
print(f'Saved: {datadir}/{run_number}_{chrom}_person0_full.txt')
print(f'Saved: {datadir}/{run_number}_{chrom}_person1_full.txt')

# Combine person0 and person1 into single vcf file
output = pd.concat([output0,output1])
output=output.sort_values(by='name')
output = output.drop_duplicates(subset=['name'])
output.to_csv(f'{datadir}/{run_number}_{chrom}_person_full.txt', index=None, header=None, sep='\t')
print(f'Saved: {datadir}/{run_number}_{chrom}_person_full.txt')