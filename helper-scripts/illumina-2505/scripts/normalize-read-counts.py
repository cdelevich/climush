import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import warnings, re
import matplotlib.pyplot as plt

## REGEX ###############################################################################################################

READ_SIZE_REGEX = r'(?<=;size=)\d{1,}(?=;)'

########################################################################################################################


## FILE PATHS ##########################################################################################################

## INPUT ##

# main input directories
illumina2505_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/')
illumina2505_data = illumina2505_main / 'data'

# sequence group pre-combining fasta files (zip compressed)
illumina2505_seqgroup_main = illumina2505_data / 'illumina_2505_per-sequence-group.zip'

## OUTPUT ##

########################################################################################################################


## FUNCTIONS ###########################################################################################################

# define a rescaling function that uses min-max normalization to normalize read counts
def min_max_normalize(read_counts):
    '''
    Normalize read counts with min-max scaling.

    :param read_counts: a list or array sequence read counts from
    a sequencing run
    :return: a list or array of normalized read counts
    '''

    ## CHECK INPUT TYPE ##

    # if already an array, do nothing
    if isinstance(read_counts, array):
        pass

    # if a list, create an array
    elif isinstance(read_counts, list):
        read_counts = np.array(read_counts)

    # if neither a list or array, raise error
    else:
        err_msg = f'input must be an array or list of read counts'
        raise TypeError(err_msg)

    ## MIN-MAX NORMALIZE READ COUNTS ##

    normalized_read_counts = (read_counts - np.min(read_counts)) / (np.max(read_counts) - np.min(read_counts))

    ## CHECK VALUE BOUNDS ##

    if all(0 <= normalized_read_counts <= 1):
        return normalized_read_counts
    else:
        normalized_min = np.min(normalized_read_counts)
        normalized_max = np.min(normalized_read_counts)
        err_msg = (f'normalized read count values are not within the bounds of [0,1]\n'
                   f'   [{normalized_min:.2f},{normalized_max:.2f}]')
        raise ValueError(err_msg)

########################################################################################################################


## GET SEQUENCE GROUP READ COUNTS ######################################################################################


## GET SEQUENCE RUN INFO ##

# get the bioinformatics run name
bioinfo_run_name = fasta_file.name.split('_')[1]


## PULL READ COUNTS FROM READ HEADERS ##

# create a dictionary where the key is the read ID, the value is the read count
read_count_dict = {}

# open the fasta file for this sequencing group...
with open(fasta_file, 'rt') as fasta_in:

    # go through each read in the fasta file...
    for record in SeqIO.parse(fasta_in, 'fasta'):

        # get read ID
        read_id_found = re.search()
        if read_id_found:
            read

        # get read count (size)
        size_found = re.search(READ_SIZE_REGEX, record.description)
        if size_found:
            read_count = int(size_found.group(0))
        else:
            warnings.warn(f'The size attribute was not located in read: {}')









# ensure that each input read has an output read count that was normalized
assert len(read_counts_norm) == len(read_counts)



## PLOT SEQUENCE GROUP READ COUNTS #####################################################################################

# plot the distribution of the read counts before and after normalizing
plt, (ax1,ax2) = plt.subplots(nrows=1, ncols=2)

ax1.hist(read_counts)
ax1.set_title('before')

ax2.hist(read_counts_norm)
ax2.set_title('after')

plt.suptitle(f'{bioinfo_run_name}: read count normalization')

savfig_output_file = (fasta_file.parent / f'{bioinfo_run_name}_read-count-norm').with_suffix('.png')
plt.savefig(savfig_output_file)

########################################################################################################################