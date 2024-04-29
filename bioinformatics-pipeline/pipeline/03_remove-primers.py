'''
FIX IDENTIFICATION OF INPUT READS
'''

from mapping import filepath_map as fpm

import argparse, sys
from climush.constants import *
from climush.bioinfo import identify_primers, remove_primers, confirm_no_primers
from climush.utilities import *

settings = get_settings(fpm)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to trim primers from
last_output = [dir for dir in fpm['pipeline-output']['prefiltered']['prefilt02_no-ambig'].glob('*')
               if re.search(f'^{NOAMBIG_PREFIX}', dir.name, re.I)][0]
is_input, illumina_files = check_for_input(file_dir=last_output, seq_platform='illumina')

if is_input:
    trimmed_path = remove_primers(illumina_files, file_map=fpm, platform='illumina', paired_end=True, verbose=False)
    confirm_no_primers(trimmed_path, file_map=fpm, platform='illumina')
else:
    pass

#####################
# PACBIO ############
#####################

# check that there are PacBio reads to trim primers from
last_output = [dir for dir in fpm['pipeline-output']['demultiplexed'].glob('*')
               if re.search(f'^{NOAMBIG_PREFIX}', dir.name, re.I)][0]
is_input, illumina_files = check_for_input(file_dir=last_output, seq_platform='illumina')

# REMOVE AFTER TESTING
is_input = True
pacbio_files = fpm['pipeline-output']['prefiltered'].glob('*.fast*')
####################

if is_input:
    remove_primers(pacbio_files, file_map=fpm, platform='pacbio', paired_end=False)
else:
    pass

#####################
# SANGER ############
#####################

# first check 'nearest' possible directory for sequences
# last_output = [dir for dir in fpm['pipeline-output']['prefiltered'].glob('*') if re.search('^no-ambig', dir.stem, re.I)][0]
# is_input, sanger_files = check_for_input(last_output)
#
# if is_input:
#     remove_primers(sanger_files, file_map=fpm, platform='sanger')
# else:
#     pass

# when all are primers trimmed, continue to next
continue_to_next(__file__, settings)