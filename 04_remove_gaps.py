# 04 - Remove Gaps
#
# This script removes gaps in the reads and the genome. It contains the following steps:
#
# - Loading and cleaning read data and reference genome data
# - Visualizing overlap between reads
# - Collapsing the regions where there is empty data
#
# The script is based on 04_remove_gaps.ipynb
#
# INPUTS:
#   run_number
#   disease
#   nanopore reads for location: e.g. run1_chr11_sca.sam
#   reference genome fasta for location: e.g. chr11_selected.fa
# 
# OUTPUTS:
#   reads without gaps: e.g. run1_chr11_sca_clean.csv
#   reference genome without gaps: e.g. run1_chr11_reference_genome.json
import sys
print('')
print('----------------------------------------------------------------------')
print(f'04 - REMOVE GAPS ({__file__})')
print(f"Run: {sys.argv[1]}, disease: {sys.argv[2]}.")
print('')
print('')

# Load Modules
# ---------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import json
import os

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

# Load and Clean Reference Genome
# ---------------------------------------------------------------------------------------------------
print(f'Reading Reference Genome from: {rootdir}/{chrom}_selected.fa')
# Read in fasta file: remove line breaks and header
def read_fasta_genome(fasta_file,chromosome_header):
    clean_data = fasta_file.read().replace("\n", "")
    clean_data = clean_data.replace(chromosome_header,"") # get rid of header
    return clean_data

with open(f'{rootdir}/{chrom}_selected.fa') as f: # update path if needed
    ref_genome = read_fasta_genome(f,f'>{chrom}')

# Load and Clean Reads
# ---------------------------------------------------------------------------------------------------
# Read Sam file
print(f'Reading Nanopore reads from:   {datadir}/{run_number}_{chrom_dis}.sam')
nanopore_reads = pd.read_csv(f'{datadir}/{run_number}_{chrom_dis}.sam',sep='\t',header=None,error_bad_lines=False, warn_bad_lines=False)

# Clean Reads
nanopore_reads = nanopore_reads.iloc[:,:11]
nanopore_reads.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT','TLEN', 'SEQ', 'QUAL']
# Sort and get ID
nanopore_reads = nanopore_reads.sort_values(by='POS',ascending=True) # sort based on starting index of reads
nanopore_reads['ID'] = np.nan
# Get columns of interest
nanopore_reads = nanopore_reads[['ID', 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT','TLEN', 'SEQ', 'QUAL']]
# Drop missing values
nanopore_reads = nanopore_reads.loc[nanopore_reads['SEQ'] != '*'] # drop any rows without a proper nucleotide sequence
# Reset index and save ID
nanopore_reads = nanopore_reads.reset_index()
nanopore_reads['ID'] = nanopore_reads.index
# Calculate length and end position
nanopore_reads['SEQ_LEN'] = nanopore_reads['SEQ'].apply(lambda x: len(x))
nanopore_reads['END_POS'] = nanopore_reads['POS']+nanopore_reads['SEQ_LEN']

# Collapse Gaps
# ---------------------------------------------------------------------------------------------------
# print('Collapsing gaps.')
ref_genome = ref_genome[min(nanopore_reads['POS']):max(nanopore_reads['END_POS'])]

# Reset start and end positions based on location of interest
shift = min(nanopore_reads['POS'])
# archive old positions
nanopore_reads['ORIG_POS'] = nanopore_reads['POS']
nanopore_reads['ORIG_END_POS'] = nanopore_reads['END_POS']
# shift positions
nanopore_reads['POS'] = nanopore_reads['POS']-shift
nanopore_reads['END_POS'] = nanopore_reads['END_POS']-shift
nanopore_reads = nanopore_reads.reset_index() # to make the next step easier
nanopore_reads = nanopore_reads[['ID', 'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR','RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'SEQ_LEN', 'END_POS','ORIG_POS', 'ORIG_END_POS']]

### Remove gaps
confirmed_gaps = [] # Collect gaps in list
max_end_pos = 0 # Store end position of longest read
# Initialize island counter
islandID_counter=0
nanopore_reads["islandID"]=islandID_counter  

for index in nanopore_reads.index:
 
    current_read_end = nanopore_reads.loc[index,'END_POS']
    
    # Current read is redundant, skip
    if current_read_end <= max_end_pos: 
        continue
        
    # This is a "standard" read (not embedded)    
    elif current_read_end > max_end_pos: 
 
        current_read_start = nanopore_reads.loc[index,'POS']
        
        # There is no gap, update max_end_pos
        if current_read_start <= max_end_pos: 
            max_end_pos = current_read_end 
            continue
            
        # There is a gap, remove it  
        elif current_read_start > max_end_pos: 
            gap_start = max_end_pos+1
            gap_end = current_read_start-1
            shift = gap_end - gap_start + 1

            nanopore_reads.loc[index:,'POS'] = nanopore_reads.loc[index:,'POS']-shift
            nanopore_reads.loc[index:,'END_POS'] = nanopore_reads.loc[index:,'END_POS']-shift

            # Calculate gap based on original position
            orig_gap_end   = nanopore_reads.loc[index,'ORIG_POS'] - 1
            orig_gap_start = orig_gap_end - shift + 1
            orig_gap = (orig_gap_start,orig_gap_end)
            confirmed_gaps.append(orig_gap) # Add gap to list
     
            current_read_end = nanopore_reads.loc[index,'END_POS']
            max_end_pos = current_read_end # Current read end should always be larger than max_end_pos
            
            islandID_counter += 1
            nanopore_reads.loc[index:,"islandID"]=islandID_counter 
            
# print(f"There are {len(confirmed_gaps)} gaps.")

# Check if gaps overlap
previous_gap_end = 0
for gap in confirmed_gaps:
    gap_start = gap[0]
    
    if gap_start <= previous_gap_end: # we have a problem
        print(f"Found an intersection in gap {gap}.")
    
    # Update gap_end for next iteration
    previous_gap_end = gap[1]

### Adjust reference_genome based on removed indices
# Reset confirmed_gaps to start at 0
shift = min(nanopore_reads['ORIG_POS'])
orig_confirmed_gaps = [(gap[0]-shift, gap[1]-shift) for gap in confirmed_gaps]

indices_to_remove = []
for gap in orig_confirmed_gaps:
    for i in range(gap[0],gap[1]+1):
        indices_to_remove.append(i)

print(f"Removing {len(confirmed_gaps)} gaps = {len(indices_to_remove):,} bp.")
print(f"Reference genome length before collapsing gaps: {len(ref_genome):,} bp.")
print(f"Final genome length should be:                  {len(ref_genome)-len(indices_to_remove):,} bp.")

final_genome = ''
segment_start = 0

for gap in orig_confirmed_gaps:

    gap_start = gap[0]
    gap_end = gap[1]
    
    # Cut segment until next gap
    new_segment = ref_genome[segment_start:gap_start]
    
    # Append to final genome
    final_genome = final_genome + new_segment
    
    # Start of next segment
    segment_start = gap_end + 1
    
# Add last segment
new_segment = ref_genome[segment_start:]
final_genome = final_genome + new_segment

print(f"Final genome length after collapsing gaps:      {len(final_genome):,} bp.")

# Save output
# ---------------------------------------------------------------------------------------------------
# Save Reference Genome
initial_position = min(nanopore_reads['ORIG_POS'])
final_genome_dict = {'initial_position':initial_position, 'reference_genome': final_genome}
with open(f'{datadir}/{run_number}_{chrom}_reference_genome.json', 'w') as f:
    json.dump(final_genome_dict,f)

# Save Nanopore Reads
nanopore_reads.to_csv(f'{datadir}/{run_number}_{chrom_dis}_clean.csv', index=False)

# Save original confirmed gaps

print(f'Saved: {datadir}/{run_number}_{chrom}_reference_genome.json')
print(f'Saved: {datadir}/{run_number}_{chrom_dis}_clean.csv')