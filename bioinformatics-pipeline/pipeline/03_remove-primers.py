'''
FIX IDENTIFICATION OF INPUT READS
'''

from mapping import filepath_map as fpm

import argparse, sys, pathlib
from climush.constants import *
from climush.bioinfo import identify_primers, remove_primers, confirm_no_primers
from climush.utilities import *

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Remove primers using cutadapt.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

parser.add_argument('-i', '--input',
                    default=fpm['pipeline-output']['prefiltered']['prefilt02_no-ambig'] / f'no-ambig_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to the sequencing files. Will default to the location that is '
                         'expected with the Docker container\'s native file structure.')

# the action is what will occur if the flag is used
# if check_only flag used, then it will be True; if not used, then False
parser.add_argument('--check-only', action='store_true',
                    help='Boolean (True/False). Whether to double check that the primers are removed from the sequences '
                         'in their forward and reverse complement direction.')

# if option provided, will print out all warnings
parser.add_argument('--verbose',
                    action='store_true',
                    default=settings['automate']['verbose'],
                    help='Whether to print out all errors and warnings in the terminal; if True, will print all '
                         'errors and warnings; if False, will only print out fatal errors.')

args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# check that there are Illumina reads to trim primers from
is_input, illumina_files = check_for_input(
    args['input'],
    config_dict=settings,
    file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
)

if is_input:
    if args['check_only']:
        print(f'Only checking for primers...')
        confirm_no_primers(args['input'], file_map=fpm, platform=platform)
    else:
        print(f'Running cutadapt...\n')
        trimmed_path = remove_primers(illumina_files, file_map=fpm, platform=platform, paired_end=True, verbose=args['verbose'])
        confirm_no_primers(trimmed_path, file_map=fpm, platform=platform)
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# check that there are PacBio reads to trim primers from
is_input, pacbio_files = check_for_input(
    args['input'],
    config_dict=settings,
    file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
)

if is_input:
    remove_primers(input_files=pacbio_files, file_map=fpm, platform=platform, verbose=args['verbose'], paired_end=False)
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# check if there are Sanger sequences to trim primers from
# is_input, sanger_files = check_for_input(
#     args['input'],
#     config_dict=settings,
#     file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
# )
#
# if is_input:
#     remove_primers(sanger_files, file_map=fpm, platform=platform)
# else:
#     pass

# when all are primers trimmed, continue to next
continue_to_next(__file__, settings)