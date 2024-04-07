from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
from climush.bioinfo import check_chimeras
from climush.utilities import *

settings = import_config_as_dict(file_path=fpm['config']['main'], file_handle='pipeline-settings')
run_name = settings['run_details']['run_name']

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to dereplicate
last_output = [dir for dir in fpm['pipeline-output']['chimera-checked'].glob('*') if re.search(f'^{CHIM}_{run_name}', dir.stem, re.I)][0]
is_input, illumina_files = check_for_input(last_output)

if is_input:
    check_chimeras(input_files=illumina_files, file_map=fpm, ref=None)
else:
    pass

#####################
# PACBIO ############
#####################

# last_output = [dir for dir in fpm['pipeline-output']['derep-full-length'].glob('*') if re.search(f'^{DEREP_PREFIX}01', dir.stem, re.I)][0]
# is_input, pacbio_files = check_for_input(last_output)
#
# if is_input:
#     separate_subregions(input_files=pacbio_files, file_map=fpm)
# else:
#     pass

#####################
# SANGER ############
#####################

# continue to next
continue_to_next(Path(__file__), settings)