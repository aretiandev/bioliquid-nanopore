# 10 - STRspy clean
#
# This script edits the output of STRspy to the right format with additional columns. The script is based on 07_STRspy.ipynb
#
# INPUTS:
#   run_number
#   disease
#   Full list of STRs
# 
# OUTPUTS:
#   VCF file with STRs present: {datadir}/{run_number}_{chrom}_person_full.txt
print('')
print('10 - STRSPY CLEAN')

# Load Modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import os
import sys
import subprocess

# Set Variables
# ---------------------------------------------------------------------------------------------------
run_num = sys.argv[1]
dis = sys.argv[2]
print(f"Run: {run_num}, disease: {dis}.")
print('')

from src.setup_variables import *
try:
    chrom=dis_data[dis]['chr']
    location=dis_data[dis]['location']
except:
    print("Disease should be in disease list: cf, sca, sma1, sma2, thal1, thal2, thal3, pompe.")
    
run_number=f"run{run_num}"
chrom_dis=f"{chrom}_{dis}"
datadir=f"{rootdir}/{run_number}"

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
def complete_str_df(person):
    # Load STRspy output
    strspy_df = pd.read_csv(f'{datadir}/strspy/{dis}/output/Countings/{run_number}_{chrom}_{person}_strs.txt', sep='\t', header=None)
    strspy_df.columns = ['name', 'count', 'normcount']

    # Load Full STR list
    df = pd.read_csv(f'{datadir}/hg38.hipstr_reference_full_strs.bed', sep='\t', header=None)
    df.columns=['chr','start','end','NA','repeats','name','motif','str']

    # Append it to STRspy output
    output = strspy_df.merge(df, how='left', on='name')
    output = output[['name','count','chr','start','end','motif', 'str']]

    # Save
    return output
    
print('Person0: Adding motif and full STR columns.')
output0 = complete_str_df('person0')
print('Person1: Adding motif and full STR columns.')
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