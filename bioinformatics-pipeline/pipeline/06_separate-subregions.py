from mapping import filepath_map as fpm

import argparse, pathlib
from pathlib import Path
from climush.constants import *
from climush.bioinfo import separate_subregions, concat_regions, check_concat_output
from climush.utilities import *

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Identify and separate the ITS/LSU subregions.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to dereplicate
parser.add_argument('-i', '--input',
                    default=fpm['pipeline-output']['derep-full-length'] / f'derep01_{run_name}',
                    type=pathlib.PurePath,
                    help='The path to a directory containing sequencing files to separate. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

parser.add_argument('-c', '--concat-only',
                    action='store_true',
                    help='If this flag is used, ITSx will not run but the ITSx output files will be accessed to '
                         'concatenate subregions.')

# input directory containing the files to dereplicate
parser.add_argument('--concat-in',
                    default=fpm['pipeline-output']['separated-subregions'] / f'itsx_{run_name}',
                    type=pathlib.PurePath,
                    help='The path to a directory containing sequencing files to concatenate. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# parse default or CL arguments
args = vars(parser.parse_args())

# ## REMOVE AFTER TESTING ########################################################################################
# args = {'input': fpm['pipeline-output']['derep-full-length'] / f'derep01_{run_name}',
#         'concat_only': False,
#         'concat_in': fpm['pipeline-output']['separated-subregions'] / f'itsx_{run_name}'}
# ################################################################################################################

#####################
# ILLUMINA ##########
#####################

# check that there are Illumina reads to dereplicate
# last_output = [dir for dir in fpm['pipeline-output']['derep-full-length'].glob('*') if re.search(f'^{DEREP_PREFIX}01', dir.stem, re.I)][0]
# is_input, illumina_files = check_for_input(last_output)
#
# if is_input:
#     separate_subregions(input_files=illumina_files, file_map=fpm)
# else:
#     pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

if args['concat_only']:
    concat_regions(dir_path=concat_path)
    check_concat_output(itsx_dir=concat_path, full_len_dir=input_path, num_bp_compare=50)
else:
    is_input, pacbio_files = check_for_input(file_dir=input_path, seq_platform=platform)
    if is_input:
        itsx_out_path = separate_subregions(input_files=pacbio_files, file_map=fpm)
        concat_regions(dir_path=itsx_out_path)
        check_concat_output(itsx_dir=itsx_out_path, full_len_dir=input_path, num_bp_compare=50)
    else:
        pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)