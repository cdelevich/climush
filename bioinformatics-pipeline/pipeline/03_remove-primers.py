'''
FIX IDENTIFICATION OF INPUT READS
'''

from mapping import filepath_map as fpm

import argparse, sys
from climush.constants import *
from climush.bioinfo import identify_primers, remove_primers
from climush.utilities import *

settings = get_settings(fpm)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to pre-filter
last_output = [dir for dir in fpm['pipeline-output']['prefiltered'].glob('*') if re.search('^no-ambig', dir.stem, re.I)]
is_input, illumina_files = check_for_input(last_output)

if is_input:
    remove_primers(illumina_files, file_map=fpm, platform='custom', paired_end=True)
else:
    pass

#####################
# PACBIO ############
#####################

# first check 'nearest' possible directory for sequences
# last_output = [dir for dir in fpm['pipeline-output']['prefiltered'].glob('*') if re.search('^no-ambig', dir.stem, re.I)][0]
# is_input, pacbio_files = check_for_input(last_output)

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
continue_to_next(Path(__file__), settings)