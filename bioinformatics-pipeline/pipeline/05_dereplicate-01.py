from mapping import filepath_map as fpm

import argparse, sys, subprocess, pathlib
from pathlib import Path
from climush.constants import *
from climush.bioinfo import dereplicate
from climush.utilities import *

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Dereplicate full-length reads.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to dereplicate
parser.add_argument('-i', '--input',
                    default=None,  # depends on previous pipeline steps, settings
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to dereplicate. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

args = vars(parser.parse_args())

## REMOVE AFTER TESTING
args = {'input': None}

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# where should the input be located for Illumina reads?
if args['input'] is None:

    # which folder from quality-filtered should be used (i.e., merged or unmerged)
    if settings['quality_filtering'][platform]['merge_reads']:  # if reads were merged
        dir_name = f'{MERGED_PREFIX}_{run_name}'
    else:
        dir_name = f'{QUALFILT_PREFIX}'

    # make full path with dir
    input_path = fpm['pipeline-output']['quality-filtered'] / dir_name

else:

    input_path = args['input']


# check that there are Illumina reads to dereplicate
is_input, illumina_files = check_for_input(file_dir=input_path, seq_platform=platform)

if is_input:
    print(f'Dereplicating {len(illumina_files)} {platform} reads...\n')

    dereplicate(input_files=illumina_files, derep_step=1, platform='illumina', file_map=fpm)
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# where should the input be located for PacBio reads?
if args['input'] is None:
    input_path = fpm['pipeline-output']['primers-trimmed'] / f'trim_{run_name}'
else:
    input_path = Path(args['input'])

is_input, pacbio_files = check_for_input(file_dir=input_path, seq_platform=platform)

if is_input:
    dereplicate(input_files=pacbio_files, derep_step=1, platform=platform, file_map=fpm)
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)