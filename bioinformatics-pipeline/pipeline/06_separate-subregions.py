from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
from climush.bioinfo import separate_subregions
from climush.utilities import *

settings = get_settings(fpm)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to dereplicate
# last_output = [dir for dir in fpm['pipeline-output']['derep-full-length'].glob('*') if re.search(f'^{DEREP_PREFIX}01', dir.stem, re.I)][0]
# is_input, illumina_files = check_for_input(last_output)
#
# if is_input:
#     separate_subregions(input_files=illumina_files, file_map=fpm)
# else:
#     pass

#####################
# PACBIO ############
#####################

last_output = [dir for dir in fpm['pipeline-output']['derep-full-length'].glob('*') if re.search(f'^{DEREP_PREFIX}01', dir.stem, re.I)][0]
is_input, pacbio_files = check_for_input(last_output)

if is_input:
    separate_subregions(input_files=pacbio_files, file_map=fpm)
else:
    pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)