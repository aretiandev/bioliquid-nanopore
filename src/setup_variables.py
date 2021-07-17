# SETUP VARIABLES
#
# This script sets up the common variables in all python scripts.
#
# INPUTS:
#   disease
# 
# OUTPUTS:
#   chr
#   location
#   directories: rootdir, datadir

# Set Variables
# ---------------------------------------------------------------------------------------------------
dis_data = {
    'cf'   :{'chr':'chr7','location':''},
    'sca'  :{'chr':'chr11','location':5227002},
    'sma1' :{'chr':'chr5','location':''},
    'sma2' :{'chr':'chr5','location':''},
    'thal1':{'chr':'chr16','location':''},
    'thal2':{'chr':'chr16','location':''},
    'thal3':{'chr':'chr16','location':''},
    'pompe':{'chr':'chr17','location':25000000}
}
rootdir=f"/mnt/aretian/genomics/nanopore"
