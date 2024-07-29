from mapping import filepath_map as fpm

import argparse, sys, re, json, pathlib
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
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to demultiplex
parser.add_argument('-i', '--input',
                    default=fpm['sequences']['demux'],
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to demultiplex. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# output directory to send the demultiplexed reads to
parser.add_argument('-o', '--output',
                    default=fpm['pipeline-output']['demultiplexed'],
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to demultiplex. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# if option provided, will use this minimum score instead of configuration settings
parser.add_argument('--min-score', nargs='?',
                    default=settings['demultiplex']['min_precision'],
                    help='The minimum precision of detecting barcodes in the multiplexed reads. A minimum score '
                         'of 80 yields a precision of >99.99%. The higher the score, the lower the level of'
                         'contaminant, but this will also lead to a lower post-demultiplexing read yield. The '
                         'default minimum precision score is 93, as recommended by Tedersoo et al. 2021.')

# if option provided, will print out all warnings
parser.add_argument('--verbose',
                    action='store_true',
                    default=settings['automate']['verbose'],
                    help='Whether to print out all errors and warnings in the terminal; if True, will print all '
                         'errors and warnings; if False, will only print out fatal errors.')

args = vars(parser.parse_args())

# parse default or CL arguments

# I don't have a way to accept a non-default minimum precision score yet (--min-score) so print warning if provided one
if not isinstance(args['min_score'], dict):
    default_precision = settings['demultiplex']['min_precision']
    print(f'WARNING. This version of {Path(__file__).name} is not yet configured to accept values for '
          f'--min-score that differ from the default values from the configuration file. Using the default '
          f'values of:\n'
          f'   {default_precision.items()}')

# set variable that will only switch to True if files requiring demultiplexing are detected
files_demuxed = False

#####################
# ILLUMINA ##########
#####################
# platform = 'illumina'
#
# # check if there are Illumina reads that need to be demultiplexed
# is_input, illumina_files = check_for_input(args['input'], seq_platform=platform)
#
# if is_input:
#     msg = f'WARNING. Currently, there is no script that can demultiplex Illumina reads. You can continue with the ' \
#           f'pipeline, but these {len(illumina_files)} multiplexed Illumina sequencing files will be ignored. Do you ' \
#           f'wish to continue without these sequences?'
#     prompt_yes_no_quit(message = msg)
# else:
#     pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check if there are PacBio reads that need to be demultiplexed
is_input, pacbio_files = check_for_input(args['input'], seq_platform=r'\d{4}')

if is_input:
    files_demuxed = True
    print(f'{len(pacbio_files)} PacBio sequencing files were detected that require demultiplexing...')

    # create output directory for demultiplexed samples
    mkdir_exist_ok(new_dir = args['output'])

    # demultiplex input files
    demultiplex(output_dir=args['output'], file_map=fpm, verbose=args['verbose'], multiplexed_files=pacbio_files)
else:
    pass

#####################
# SANGER ############
#####################
# platform = 'sanger'
#
# # check if there are Sanger reads that need to be demultiplexed
# is_input, sanger_files = check_for_input(args['input'], seq_platform=platform)
#
# if is_input:
#     msg = f'WARNING. Currently, there is no script that can demultiplex Sanger reads. You can continue with the ' \
#           f'pipeline, but these {len(sanger_files)} multiplexed Sanger sequencing files will be ignored. Do you ' \
#           f'wish to continue without these sequences?'
#     prompt_yes_no_quit(message = msg)
# else:
#     pass

#########################

# when all are demultiplexed (if possible), continue to the next script
if not files_demuxed:  # print if no demultiplexing was carried out
    print(f'No sequencing files were detected in the {args["input"]} that require demultiplexing.\n')
continue_to_next(__file__, settings)