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
                    type=pathlib.PosixPath,
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
                    type=pathlib.PosixPath,
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

# may want to use this to confirm that reads are ITS1? could check orientation as well?
# I don't think there's a use for concatenating regions for Illumina? at least for CliMush use

#####################
# PACBIO ############
#####################
platform = 'pacbio'

if args['concat_only']:
    for itsx_sample in args['concat_in'].glob('*'):
        print(f'\n{itsx_sample.stem}\n')
        concat_regions(dir_path=itsx_sample, regions_to_concat=['ITS1', '5_8S', 'ITS2'])  # full ITS
        concat_regions(dir_path=itsx_sample, regions_to_concat=['ITS1', '5_8S', 'ITS2', 'LSU'])  # full length read (reoriented)
        # check_concat_output(itsx_dir=concat_path, full_len_dir=args['input'], num_bp_compare=50)
else:
    is_input, pacbio_files = check_for_input(file_dir=args['input'], seq_platform=platform)
    if is_input:
        itsx_out_path = separate_subregions(input_files=pacbio_files, file_map=fpm)
        concat_regions(dir_path=itsx_out_path, regions_to_concat=['ITS1', '5_8S', 'ITS2'])  # full ITS
        concat_regions(dir_path=itsx_out_path, regions_to_concat=['ITS1', '5_8S', 'ITS2', 'LSU'])  # full length read (reoriented)
        # check_concat_output(itsx_dir=itsx_out_path, full_len_dir=args['input'], num_bp_compare=50)
    else:
        pass

#####################
# SANGER ############
#####################


# when all are seqs are dereplicated, continue to next
continue_to_next(__file__, settings)