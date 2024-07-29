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
                    default=fpm['pipeline-output']['primers-trimmed'] / f'trim_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to dereplicate. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################
# platform = 'illumina'
#
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
platform = 'pacbio'

is_input, pacbio_files = check_for_input(file_dir=args['input'], seq_platform=platform)

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