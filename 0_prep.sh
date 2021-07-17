# PREP script
# This script should be run before everything else
# It creates the base files that are inputs to the pipeline
# - Get reference genome and index
# - Extract locations of interest and index
# - Add columns to list of STRs

# Set variables
run_number="run1"
chrom="chr7"
dis="sca"
# Setup
chrom_dis="${chrom}_${dis}"
datadir="/mnt/aretian/genomics/nanopore/${run_number}"

# 04_remove_gaps
#-------------------------------------------------------------------------------------------

# Pull reference genome from S3
# !s3cmd get s3://aretian-genomics/nanopore/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# !gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Index reference genome
# !samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# 07_strspy
#-------------------------------------------------------------------------------------------
# Create Hipstr reference file with full STRs
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
df.to_csv(f'{datadir}/hg38.hipstr_reference_full_strs.bed', sep='\t', header=None, index=None)