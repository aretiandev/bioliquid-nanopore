# 06 - Read Assignment
#
# This script makes sure the clusters do not flip as the sliding window moves to the right. The script is based on 06_read_assignment.ipynb
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
print('')
print('06 - READ ASSIGNMENT')
print('')

# Load Modules
# -----------------------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

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

# Load data
# -----------------------------------------------------------------------------

# Load read_clusters
print(f'Reading clusters from:             {datadir}/{run_number}_{chrom}_read_clusters.txt')
print(f'Reading clean Nanopore reads from: {datadir}/{run_number}_{chrom_dis}_clean.csv')

results_file = pd.read_csv(f"{datadir}/{run_number}_{chrom}_read_clusters.txt")

results_file.columns = ['ID', 'kmeans_cls2', 'window_num']

# Load Nanopore Reads
nanopore_reads = pd.read_csv(f"{datadir}/{run_number}_{chrom_dis}_clean.csv") # load nanopore reads with padded sequences and other attribute data

# Add column(s)
nanopore_reads['individual'] = np.nan
results_file['individual'] = np.nan

# Helper functions
# -----------------------------------------------------------------------------
def propagate_assignment_1(read):
    print('_', end="", flush=True)
    
    ID = read['ID']
    individual = read['individual']
    
    # Write results
    nanopore_reads.loc[nanopore_reads.ID == ID, 'individual'] = individual # propagate the person assignment to the main reads dataframe
    results_file.loc[results_file['ID'] == ID, 'individual'] = individual # propagate the person assignment to the main reads dataframe
    ID_person_match[ID] = individual # add result to dictionary of results

def propagate_assignment_2(read):
#     print('Inside propagate 2')
    ID = read['ID']
    individual = read['individual']
#     print('Looking for readID in history')

    if ID in ID_person_match.keys():
        # ID has already been assigned 
        # NOTE: same as "if np.isnan(read_result['individual']) == False"
        print("-", end="", flush=True)
        return None
    
#     print('did not find readID in history')
    # Get cluster assignment of read: 0 or 1
#     cluster_num = window_reads.loc[window_reads['ID'] == ID, 'kmeans_cls2'].values[0]
    cluster_num = read['kmeans_cls2']

    # see the individual assigned to IDs in this cluster in this particular window
    cluster_family = window_reads.loc[window_reads['kmeans_cls2'] == cluster_num]
    cluster_family_notna = cluster_family.loc[cluster_family['individual'].notnull()]
    assigned_individual = cluster_family_notna['individual'].unique()
#     print(f'Assigned individual: {assigned_individual}')


    if len(assigned_individual) > 1:
        # this cluster has mixed assignment which is bad
        # we expect that if IDs are already assigned they were in the same cluster in the previous window(s)
        # which means that this issue should theoretically never happen (NOTE: think about/check this logic)
#         print('ERROR: Found a window with person0 and person1 in the same cluster.')
#         print("E", end="")
#         return None
#         print('Multiple assigned individuals in cluster.')

#         # Make a call for this reads with mixed assignment
        mean_individual_to_assign = cluster_family_notna['individual'].mean()
        if mean_individual_to_assign == 0.5:
            if loop_direction=='left':
                individual_to_assign = cluster_family_notna['individual'].iloc[-1]
            else:
                individual_to_assign = cluster_family_notna['individual'].iloc[0]
            print("T", end="", flush=True)
        else:
            individual_to_assign = round(mean_individual_to_assign) 
            print("E", end="", flush=True)
            
#         # Write results
        nanopore_reads.loc[nanopore_reads['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
        results_file.loc[results_file['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
        ID_person_match[ID] = individual_to_assign # add result to dictionary of results
        

    elif len(assigned_individual) == 0:
#         print('All reads in cluster are unassigned.')

        # All reads in the cluster_family are unassigned (NaN)
        # Look at what is the assigned individual of the opposite cluster
        other_cluster_family = window_reads.loc[window_reads['kmeans_cls2'] == (1-cluster_num)]
        other_cluster_family_notna = other_cluster_family.loc[other_cluster_family['individual'].notnull()]
        other_assigned_individual = other_cluster_family_notna['individual'].unique()
#         print(f'New read. Other assigned individual: {other_assigned_individual[0]}')
        print("o", end="", flush=True)
        
        # Assign the opposite individual
        individual_to_assign = 1-other_assigned_individual[0]
        
        # Bulk write results
        for ID in cluster_family['ID'].unique(): # For every read in the same cluster
            if nanopore_reads.loc[nanopore_reads['ID'] == ID, 'individual'].isna().values[0]:
                nanopore_reads.loc[nanopore_reads['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
                results_file.loc[results_file['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
                ID_person_match[ID] = individual_to_assign # add result to dictionary of results
        
    elif len(assigned_individual) == 1:
#         print('Reads in cluster have already been assigned.')

        print("s", end="", flush=True)
        
        # Reads in the same cluster have already been assigned
        individual_to_assign = assigned_individual[0] # check this syntax

        # Bulk write results
        for ID in cluster_family['ID'].unique():
            if nanopore_reads.loc[nanopore_reads['ID'] == ID, 'individual'].isna().values[0]:
                nanopore_reads.loc[nanopore_reads['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
                results_file.loc[results_file['ID'] == ID, 'individual'] = individual_to_assign # propagate the person assignment to the main reads dataframe
                ID_person_match[ID] = individual_to_assign # add result to dictionary of results

# Check Assignment
# -----------------------------------------------------------------------------
print("")
print('Assigning individuals to reads.')

# Setup
ID_person_match = {} #Dictionary to keep track of kmeans cluster and individual assignment
first_iter = True

# Feedback to user
n_windows = results_file['window_num'].nunique()
print(f"Total windows: {n_windows}")
print("_: First read.")
print("-: Individual already assigned.")
print("s: New read. Sibling already assigned.")
print("o: New read. Read in other cluster already assigned.")
print("E: New read. Siblings in both clusters (tie). Assign average rounded to 0 or 1.")
print("T: New read. Siblings in both clusters (tie). Assign individual from oldest read.")
print("")

print(f"Current window:", end="", flush=True)

# Identify idmax of window
window_idxmax = results_file.groupby('window_num').count().idxmax().values[0]
# Loop to the left
loop_direction = 'left'

for window in range(window_idxmax+1,-1,-1):
# for window in results_file['window_num'].unique(): # iterate through instances of the sliding window 
    
    print(f" {window}", end="", flush=True)
    
    window_reads = results_file.loc[results_file['window_num']==window].copy() # get data only from this window
    
#     # If window does not exist, skip
#     if len(window_reads)==0:
#         print("1", end="")
#         continue
    
    # We now have IDs, cluster decisions (0,1) for every read that was included in this window.
    
    # Step 1: We make a decision in the first window
    if first_iter == True:
        
        # Make decision in first window
        window_reads.loc[window_reads['kmeans_cls2'] == 0, 'individual'] = 0
        window_reads.loc[window_reads['kmeans_cls2'] == 1, 'individual'] = 1
        
        window_reads.apply(lambda x: propagate_assignment_1(x), axis=1) # add results to main reads file
        
        first_iter = False
        continue
    
#     print('Stepping inside propagate 2...')
    # Step 2: We perform ID bookeeping and cluster association to manage all following windows
    # Book-keeping:
    # Step 1: make sure IDs match for ones that are already claimed
    # Step 2: if an ID has not been assigned determine its cluster (0,1) assign the associated individual based on already assigned reads in the cluster
    window_reads.apply(lambda x: propagate_assignment_2(x), axis=1) # add results to main reads file
    
# Loop right
loop_direction = 'right'
max(results_file['window_num'])

for window in range(window_idxmax+2,max(results_file['window_num'])+1,1):
    
    print(f" {window}", end="", flush=True)
    
    window_reads = results_file.loc[results_file['window_num']==window].copy() # get data only from this window
    
    # We now have IDs, cluster decisions (0,1) for every read that was included in this window.
    
    # Step 1: We make a decision in the first window
#     if first_iter == True:
        
#         # Make decision in first window
#         window_reads.loc[window_reads['kmeans_cls2'] == 0, 'individual'] = 0
#         window_reads.loc[window_reads['kmeans_cls2'] == 1, 'individual'] = 1
        
#         window_reads.apply(lambda x: propagate_assignment_1(x), axis=1) # add results to main reads file
        
#         first_iter = False
#         continue
        
    # Step 2: We perform ID bookeeping and cluster association to manage all following windows
    # Book-keeping:
    # Step 1: make sure IDs match for ones that are already claimed
    # Step 2: if an ID has not been assigned determine its cluster (0,1) assign the associated individual based on already assigned reads in the cluster
    window_reads.apply(lambda x: propagate_assignment_2(x), axis=1) # add results to main reads file

print("")
print("Done.")

# Get list of IDs
# -----------------------------------------------------------------------------
print("Getting lists of read IDs for each person.")
# Create ID column
# nanopore_reads['uniqueid']=nanopore_reads['QNAME']+'\t'+nanopore_reads['FLAG'].apply(str)+'\t'+nanopore_reads['RNAME'].apply(str)+'\t'+nanopore_reads['POS'].apply(str)
nanopore_reads['uniqueid']=nanopore_reads['QNAME'] # Fix later to replace with unique ids as above

# Write two text files
for person in [0,1]:
    # Get IDs of Person
    personuniqueids = nanopore_reads.loc[nanopore_reads['individual']==person, 'uniqueid']
    # Save to text file
    np.savetxt(f"{datadir}/{run_number}_{chrom}_person{person}_uniqueids.txt", personuniqueids.values , fmt='%s')
    print(f"Saved: {datadir}/{run_number}_{chrom}_person{person}_uniqueids.txt")