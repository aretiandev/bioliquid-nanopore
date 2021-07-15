# Remove Gaps

# This script removes gaps in the reads and the genome. It contains the following steps:

# - Loading and cleaning read data and reference genome data
# - Visualizing overlap between reads
# - Collapsing the regions where there is empty data

# Load Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os

# Set Variables
# Input variables
run_number="run1"
chrom="chr11"
dis="sca"
# Setup
chrom_dis=f"{chrom}_{dis}"
datadir=f"/mnt/aretian/genomics/nanopore/{run_number}"

os.environ["run_number"]=run_number
os.environ["chrom_dis"]=chrom_dis
os.environ["datadir"]=datadir

# Load and Clean Reference Genome