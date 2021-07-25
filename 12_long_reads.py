# 12 - LONG READS
#
# This script takes the information in tagged_reads.csv and combines the reads into long reads. The script is based on 12_long_reads.ipynb
#
# INPUTS:
#   run_number
#   disease
#   tagged_reads.csv
# 
# OUTPUTS:
#   tagged_long_reads.csv
import sys
print('')
print('----------------------------------------------------------------------')
print(f'12 - LONG READS ({__file__})')
print(f"Run: {sys.argv[1]}, disease: {sys.argv[2]}.")
print('')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

# Set Variables
# -----------------------------------------------------------------------------
run_num = sys.argv[1]
dis = sys.argv[2]

from src.setup_variables import *
try:
    chrom=dis_data[dis]['chr']
    location=dis_data[dis]['location']
except:
    print("Disease should be in disease list: cf, sca, sma1, sma2, thal1, thal2, thal3, pompe.")

# Setup
run_number=f"run{run_num}"
chrom_dis=f"{chrom}_{dis}"
rootdir=f"/mnt/aretian/genomics/nanopore"
datadir=f"{rootdir}/{run_number}"

# Parameters
expected_gap_fraction = 0.5

# Load Data
# -----------------------------------------------------------------------------
df = pd.read_csv(f'{datadir}/{run_number}_{chrom}_tagged_reads.csv')
original_columns = df.columns
df = df.rename(columns={'startpos':'orig_pos'})
shift = min(df['orig_pos'])
df['pos'] = df['orig_pos'] - shift
df['seq_len']= df['read'].apply(lambda x: len(x))
df['end_pos']=df['pos']+df['seq_len']-1
df['sample'] = df['samplename'].str[-1].astype(int)

# Helper functions
# -----------------------------------------------------------------------------

# Descriptive Stats
def descriptive_stats(df):
    total_length = max(df['end_pos'])
    reads_0 = df[df['sample']==0]
    reads_1 = df[df['sample']==1]

    # Number of reads
    n_0=len(reads_0)
    n_1=len(reads_1)

    # Average length
    mean_length = df['seq_len'].mean()
    mean_length_0 = reads_0['seq_len'].mean()
    mean_length_1 = reads_1['seq_len'].mean()

    # Total bins
    expected_gap = mean_length*expected_gap_fraction
    bin_length = round(mean_length + expected_gap)
    n_bins_exact = total_length/bin_length 
    n_bins = round(n_bins_exact)

    # Bins in  person0
    n_long_reads_0 = n_0/n_bins
    # Bins in person1
    n_long_reads_1 = n_1/n_bins

    print(f"Total length of region:     {total_length:,.0f}")
    print(f"Avg read length:               {mean_length:,.2f}")
    print(f"Avg read length + {expected_gap_fraction*100:,.0f}% padd:    {mean_length+expected_gap:,.2f}")
    print(f"Bin length:                    {bin_length:,.0f}")
    print(f"total_length/bin_length:           {n_bins_exact:,.2f}")
    print(f"Number of bins:                    {n_bins:,.2f} --> last bin is {'shorter' if (n_bins>n_bins_exact) else 'longer'}.")
    print("")
    print(f"               Count       Avg Length    N of long reads ")
    print(f"person0        {n_0:,.0f}       {mean_length:,.0f}        {n_long_reads_0:,.0f}")   
    print(f"person1        {n_1:,.0f}       {mean_length_0:,.0f}        {n_long_reads_1:,.0f}")
    print(f"Full sample    {len(df):,.0f}       {mean_length_1:,.0f}        {n_long_reads_0 + n_long_reads_1:,.0f}")
    
# Create bins
def create_bins():
    # Parameters
    bins = []

    last_bin_end = 0
    for i in range(int(np.floor(n_bins))):
        bin_start = last_bin_end + 1
        bin_end   = bin_start + bin_length - 1
        bins.append([bin_start, bin_end])
        last_bin_end = bin_end

    # Gap
    bin_agg_length = n_bins * bin_length
    gap = total_length - bin_agg_length
    print(f"Total length:         {total_length:,.0f}")
    print(f"Bin aggregate length: {bin_agg_length:,.0f}")
    print(f"Difference:              {gap:,.0f}")

    # Last bin length
    last_bin_length = bin_length + gap
    print(f"Bin length:              {bin_length:,.0f}")
    print(f"Last bin length:         {last_bin_length:,.0f}")

    # Adjust last bin length
    bins[-1][1]=bins[-1][0]+last_bin_length-1
    
    return bins

def assign_reads_to_bins(df, selected_person):
# Returns person_reads with the 'bin' column of bin membership
    
    reads_0 = df[df['sample']==0]
    reads_1 = df[df['sample']==1]
    
    print(f"Selected person: person{selected_person}")

    if   selected_person == 0:
        person_reads = reads_0.copy()
    elif selected_person == 1:
        person_reads = reads_1.copy()

    # Count empty bins
    empty_bins = 0
    empty_bins_list = []

    # Initialize bin column with empty lists
    person_reads['bin']=np.empty((len(person_reads), 0)).tolist()

    # Lambda function to append bin membership to list
    def append_bin(read, bin):
        read['bin'].append(bin)
        return read

    # Add bin membership column to person_reads df

    for n, bin in enumerate(bins):

        # Get reads in bin
        read_starts_inside = (person_reads['pos']    >bin[0]) & (person_reads['pos']    <bin[1])
        read_ends_inside   = (person_reads['end_pos']>bin[0]) & (person_reads['end_pos']<bin[1])
        read_covers        = (person_reads['pos']    <bin[0]) & (person_reads['end_pos']>bin[1])

        overlaps_bin = (read_starts_inside | read_ends_inside | read_covers)

        # Count empty bins
        if overlaps_bin.sum() == 0:
            empty_bins = empty_bins + 1
            empty_bins_list.append(n)
#             print(f'Empty bin: {n}')
            continue

        # Append bin number to bin column
        person_reads.loc[overlaps_bin].apply(lambda x: append_bin(x, n), axis=1)

    print(f"Reads with no assigned bin: {person_reads['bin'].isna().sum()}")
    print(f"There are {empty_bins} empty bins: {empty_bins_list}")
    
    return person_reads

def create_long_reads(person_reads):
# Populated the person_reads dataframe with the long_read column

    # List bins to ignore in last stages of the loop
    ignore_bin_list = []
    
    # Initialize long read membership
    person_reads['long_read'] = np.nan
    # Track reads that are already assigned
    person_reads['assigned'] = False

    long_read_number = 0

    while person_reads['assigned'].sum() < len(person_reads):
        # While there are unassigned reads

        n_unassigned = len(person_reads) - person_reads['assigned'].sum()
        progress = (1-n_unassigned/len(person_reads))*100

        for n, bin in enumerate(bins):
            if n in ignore_bin_list:
                continue
        # Run bin loop assigning reads
            print(f"\rLong read: {long_read_number:3.0f}. Bin: {n:3.0f}. Unassigned: {n_unassigned:5.0f}. Progress: {progress:3.0f}%", end="", flush=True)
    #         print(f"Bin number: {n}")

            # Get reads in bin
            bin_reads_boolean = person_reads.apply(lambda x: n in x['bin'], axis=1)

            # Get reads in bin that are not assigned
            bin_reads = person_reads.loc[bin_reads_boolean & ~person_reads['assigned']]

            # Get first read, assign it to long read
            try:
                selected_index = bin_reads.index[0]
            except: # The bin is empty
    #             print(f'Empty bin: {n}')
                if n not in ignore_bin_list:
                    ignore_bin_list.append(n)
                continue

            person_reads.loc[selected_index, 'long_read'] = long_read_number
    #         bin_reads.loc[0, 'long_read'] = long_read_number

            # Record assignment
            person_reads.loc[selected_index, 'assigned'] = True

            # Add check of already assigned

            # Check previous selected read

            #Store end position of last read for next loop

        # Go to the next long read
        long_read_number += 1

    print("")
    print("Done.")
    
    return person_reads
    
def plot_long_reads(person_reads):
    # Get effective number of long reads
    n_long_reads = int(person_reads['long_read'].max() + 1)

    # long_read_binary = list(np.zeros(total_length))
    long_read_collection = []
    print("Identifying points to plot long reads.")
    print(f"Total number of long reads: {n_long_reads}")

    for long_read_n in range(n_long_reads):
        long_read_binary = []

        print(f"\rProcessing long read: {long_read_n}", end="", flush=True)

        selected_long = person_reads.loc[person_reads['long_read']==long_read_n].copy()

        for point in range(0,total_length+1,10000):
            after_start =  selected_long['pos']< point
            before_end  =  point < selected_long['end_pos']
            is_in_read  = (after_start & before_end).any()

            value_to_append = 1 if is_in_read else np.nan
            long_read_binary.append(value_to_append)

        long_read_collection.append(long_read_binary)

    #     if is_in_read:
    #         long_read_binary[point]=1

    print("")
    print("Plotting.")
    # print(f"Length of read to plot: {len(long_read_binary)}")

    plt.subplots(figsize=(20,10))
    x = range(len(long_read_collection[0]))

    # for i in range(len(long_read_collection)):
    for i in range(100):
        long_read_binary = [x * (i+1) for x in long_read_collection[i]]
        plt.scatter(x, long_read_binary)
        
        
# Implementation
# -----------------------------------------------------------------------------
descriptive_stats(df)

bins = create_bins()

selected_person = 0
person_reads0 = assign_reads_to_bins(df, selected_person)
person_reads0 = create_long_reads(person_reads0)
selected_person = 1
person_reads1 = assign_reads_to_bins(df, selected_person)
person_reads1 = create_long_reads(person_reads1)


# Create output 
# -----------------------------------------------------------------------------
person_reads = pd.concat([person_reads0, person_reads1])

# Rename and keep columns as original file
person_reads_out = person_reads.copy()

# Create long read ID
person_reads_out['long_read_id'] = person_reads_out['samplename'] + '-' + person_reads_out['long_read'].astype(int).astype(str)

# Replace new read ID
person_reads_out = person_reads_out.drop(columns=['read_id'])
person_reads_out = person_reads_out.rename(columns={'long_read_id':'read_id'})

# Rename columns as original
person_reads_out = person_reads_out.rename(columns={'orig_pos':'startpos'})
person_reads_out = person_reads_out[original_columns]

# Save
# -----------------------------------------------------------------------------
person_reads_out.to_csv(f'{datadir}/{run_number}_{chrom}_tagged_reads.csv', index=None)