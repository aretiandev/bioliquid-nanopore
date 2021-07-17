#!/usr/bin/env bash
# 0 - STR List
#
# This script adds full STRs and motifs to the HipSTR reference file.
#
# INPUTS:
#   run_number
#   disease
#   HipSTR Reference file: {rootdir}/hg38.hipstr_reference.bed
# 
# OUTPUTS:
#   Expanded STR list: {rootdir}/hg38.hipstr_reference_full_strs.bed
print('')
print('0 - STR LIST')
print('')

# Load Modules
# ---------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd

rootdir=f"/mnt/aretian/genomics/nanopore"

# Create STR List
#-------------------------------------------------------------------------------------------
# Create Hipstr reference file with full STRs
# Load Full STR list
print(f"Loading STR list from HipSTR reference: {rootdir}/hg38.hipstr_reference.bed")
df = pd.read_csv(f'{rootdir}/hg38.hipstr_reference.bed', sep='\t', header=None)
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
df.to_csv(f'{rootdir}/hg38.hipstr_reference_full_strs.bed', sep='\t', header=None, index=None)
print(f'Saved: {rootdir}/hg38.hipstr_reference_full_strs.bed')