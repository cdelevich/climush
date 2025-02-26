from mapping import filepath_map as fpm

import argparse, sys, subprocess, pathlib
from pathlib import Path
from climush.constants import *
from climush.bioinfo import dereplicate
from climush.utilities import *

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Dereplicate separated subregions and full-length ITS reads.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to dereplicate
parser.add_argument('-i', '--input',
                    default=fpm['pipeline-output']['separated-subregions'] / f'itsx_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to dereplicate. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# output directory for the dereplicated files
parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The main output directory in which to create the bioinformatics run subfolder for '
                         'dereplicated sequences.')

args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################

# NOT UPDATED FOR ILLUMINA READS
# check that there are Illumina reads to dereplicate
# is_input, illumina_files = = check_for_input(
#         args['input'],
#         config_dict=settings,
#         file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
#     )
#
# if is_input:
#     dereplicate(input_files=illumina_files, derep_step=2, platform='illumina', file_map=fpm)
# else:
#     pass

#####################
# PACBIO ############
#####################

platform = 'pacbio'

# set the file_ext argument to None; this assumes itsx immediately precedes derep and you want a list of
#   directories, not a list of files (post-itsx, each sample has a directory)
is_input, pacbio_dirs = check_for_input(
    args['input'],
    config_dict=settings,
    file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform],
    file_ext=None
    )

if is_input:
    dereplicate(
        input_files=pacbio_dirs,
        derep_step=2,
        platform=platform,
        file_map=fpm,
        output_path=args['output'],
    )
else:
    pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)