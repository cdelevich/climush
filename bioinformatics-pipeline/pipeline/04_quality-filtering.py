from mapping import filepath_map as fpm

import argparse
import sys
import subprocess
from pathlib import Path
##REMOVE AFTER PACKAGE TESTING#######
sys.path.insert(0, str(fpm['package']['main']))
#####################################
from climush.constants import *
from climush.bioinfo import merge_reads, quality_filter
from climush.utilities import *

settings = import_config_as_dict(file_path=fpm['config']['main'], file_handle='pipeline-settings')

#####################
# ILLUMINA ##########
#####################

last_output = [dir for dir in fpm['pipeline-output']['primers-trimmed'].glob('*') if re.search(f'^{TRIMMED_PREFIX}', dir.stem, re.I)][0]
is_input, illumina_files = check_for_input(last_output)

if is_input:
    output = quality_filter(input_files=illumina_files, platform='illumina', file_map=fpm)
    merge_reads(input_files=output, file_map=fpm)
else:
    pass

#####################
# PACBIO ############
#####################

# pacbio_dir = PIPE_OUT_MAIN / Path(f'03_remove-primers/trim_{run_name}')
#
# # check that there are Illumina reads to pre-filter
# if input_files_present(file_path = pacbio_dir):
#     print(f'PacBio reads were detected in {pacbio_dir.parent}. '
#           f'Processing {count_files(pacbio_dir, search_for=SEQ_FILE_GLOB)} samples...\n')
#


#####################
# SANGER ############
#####################

# when all are prefiltered, continue to next
continue_to_next(Path(__file__), settings)