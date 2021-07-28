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
import numpy as np

# dis_params = data.frame('disease' = c('sca', 'cystic', 'spinal1', 'spinal2', 'thal1', 'thal2', 'thal3', 'pompe'),
#                     'chr' = c('chr11','chr7','chr5', 'chr5', 'chr16', 'chr16', 'chr16', 'chr17'),
#                     'start' = c(5227002, 117559590, 70924941, 70049523, 176680, 172876, 5225464, 0),
#                     'end'   = c(5227002, 117559590, 70953015, 70077595, 177522, 173710, 5227071, 0))

sma1  = np.mean([70924941, 70953015])
sma2  = np.mean([70049523, 70077595])
thal1 = np.mean([176680  , 177522  ]) 
thal2 = np.mean([172876  , 173710  ])
thal3 = np.mean([5225464 , 5227071 ])
                
dis_data = {
    'cf'   :{'chr':'chr7' ,'location':117559590},
    'sca'  :{'chr':'chr11','location':5227002},
    'sma1' :{'chr':'chr5' ,'location':sma1},
    'sma2' :{'chr':'chr5' ,'location':sma2},
    'thal1':{'chr':'chr16','location':thal1},
    'thal2':{'chr':'chr16','location':thal2},
    'thal3':{'chr':'chr16','location':thal3},
    'pompe':{'chr':'chr17','location':25000000}
}
location_padding=2000000
rootdir=f"/mnt/aretian/genomics/nanopore"
