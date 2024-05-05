from mapping import filepath_map as fpm

import argparse, sys, re, json
from pathlib import Path
import pandas as pd
from climush.constants import *
from climush.utilities import *
from climush.bioinfo import demultiplex

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Demultiplex PacBio reads using PacBio\'s Lima demultiplexing '
                                             'tool.',
                                 epilog='')

# input directory containing the files to demultiplex
parser.add_argument('-i', '--input',
                    default=fpm['sequences']['demux'],
                    help='The path to a directory containing sequencing files to demultiplex. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# output directory to send the demultiplexed reads to
parser.add_argument('-o', '--output',
                    default=fpm['pipeline-output']['demultiplexed'],
                    help='The path to a directory containing sequencing files to demultiplex. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# to run quality filtering and merging separately, for time restriction
parser.add_argument('--merge-from', nargs='?',
                    default=None,
                    help='If running merging separately from quality filtering, provide the path to the '
                         'directory containing the quality-filtered paired sequence files to merge.')

# if option provided, uses this minimum read length filter instead of configuration settings
parser.add_argument('--min-length', nargs='?',
                    default=settings['quality_filtering'],
                    help='The minimum read length permitted to pass the quality filter. Using this command line '
                         'argument will override the settings from your configuration file.')

# if option provided, uses this maximum read length filter instead of configuration settings
parser.add_argument('--max-length', nargs='?',
                    default=settings['quality_filtering'],
                    help='The maximum read length permitted to pass the quality filter. Using this command line '
                         'argument will override the settings from your configuration file.')

# if option provided, uses this maximum expected error filter instead of configuration settings
parser.add_argument('--max-error', nargs='?',
                    default=settings['quality_filtering'],
                    help='The maximum expected error allowed to pass the quality filter. Using this command line '
                         'argument will override the settings from your configuration file.')


args = vars(parser.parse_args())

# parse default or CL arguments
# if an input path is provided, convert to a Path object
if isinstance(args['input'], str):
    input_path = Path(args['input'])
else:
    input_path = args['input']

# set variable that will only switch to True if files requiring demultiplexing are detected
files_demuxed = False

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check if there are Illumina reads that need to be demultiplexed
is_input, illumina_files = check_for_input(fpm['sequences']['demux'], seq_platform=platform)

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
platform = 'pacbio'

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
platform = 'sanger'

# check if there are Sanger reads that need to be demultiplexed
is_input, sanger_files = check_for_input(fpm['sequences']['demux'], seq_platform=platform)

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
    print(f'No sequencing files were detected in the {input_path} that require demultiplexing.\n')
continue_to_next(__file__, settings)