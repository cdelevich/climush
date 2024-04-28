from mapping import filepath_map as fpm

import argparse, sys, re, json
from pathlib import Path
import pandas as pd
from climush.constants import *
from climush.utilities import *
from climush.bioinfo import demultiplex

settings = get_settings(fpm)

# set variable that will only switch to True if files requiring demultiplexing are detected
files_demuxed = False

#####################
# ILLUMINA ##########
#####################

# check if there are Illumina reads that need to be demultiplexed
is_input, illumina_files = check_for_input(fpm['sequences']['demux'], seq_platform='illumina')

if is_input:
    msg = f'WARNING. Currently, there is no script that can demultiplex Illumina reads. You can continue with the ' \
          f'pipeline, but these {len(illumina_files)} multiplexed Illumina sequencing files will be ignored. Do you ' \
          f'wish to continue without these sequences?'
    prompt_yes_no_quit(message = msg)
else:
    pass

#####################
# PACBIO ############
#####################

# check if there are PacBio reads that need to be demultiplexed
is_input, pacbio_files = check_for_input(fpm['sequences']['demux'], seq_platform='\d{4}')

if is_input:
    files_demuxed = True
    print(f'{len(pacbio_files)} PacBio sequencing files were detected that require demultiplexing...')

    # create output directory for demultiplexed samples
    mkdir_exist_ok(new_dir = fpm['pipeline-output']['demultiplexed'])

    # demultiplex input files
    demultiplex(file_map=fpm, multiplexed_files=pacbio_files)
else:
    pass

#####################
# SANGER ############
#####################

# check if there are Sanger reads that need to be demultiplexed
is_input, sanger_files = check_for_input(fpm['sequences']['demux'], seq_platform='sanger')

if is_input:
    msg = f'WARNING. Currently, there is no script that can demultiplex Sanger reads. You can continue with the ' \
          f'pipeline, but these {len(sanger_files)} multiplexed Sanger sequencing files will be ignored. Do you ' \
          f'wish to continue without these sequences?'
    prompt_yes_no_quit(message = msg)
else:
    pass

#########################

# when all are demultiplexed (if possible), continue to the next script
if not files_demuxed:  # print if no demultiplexing was carried out
    print(f'No sequencing files were detected in the /sequences directory that require demultiplexing.\n')
continue_to_next(Path(__file__), settings)