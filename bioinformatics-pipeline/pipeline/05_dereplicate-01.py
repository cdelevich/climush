from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
from climush.bioinfo import dereplicate
from climush.utilities import *

settings = get_settings(fpm)

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to dereplicate
# last_output = [dir for dir in fpm['pipeline-output']['quality-filtered'].glob('*') if re.search(f'^{QUALFILT_PREFIX}', dir.stem, re.I)][0]
# is_input, illumina_files = check_for_input(last_output)
#
# if is_input:
#     dereplicate(input_files=illumina_files, derep_step=1, platform='illumina', file_map=fpm)
# else:
#     pass
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
# dereplicate(input_dir=pacbio_dir, derep_step=1, platform='pacbio', script_name=__file__)

# REMOVE AFTER TESTING
# is_input = True
# pacbio_files = list((fpm['pipeline-output']['primers-trimmed'] / "trim_bioinfo-test").glob('*.fast*'))
# ####################
#
# if is_input:
#     dereplicate(input_files=pacbio_files, derep_step=1, platform='pacbio', file_map=fpm)
#     # dereplicate(pacbio_files, file_map=fpm, platform='pacbio', paired_end=False)
# else:
#     pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
# continue_to_next(__file__, settings)