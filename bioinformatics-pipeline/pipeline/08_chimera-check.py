import argparse, pathlib
from pathlib import Path

from climush.bioinfo import check_chimeras
from climush.constants import *
from climush.utilities import get_settings, check_for_input, continue_to_next

from mapping import filepath_map as fpm

# import settings from the configuration file
settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

# set up command line options
parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Detect and remove chimeras.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to check for chimeras; no default, since depends on platform
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to chimera check. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py. This depends on the sequencing platform.')

# method to use for chimera detection
parser.add_argument('--method',
                    default=settings['chimera_check']['method'],  # depends on previous pipeline steps, settings
                    choices=['denovo', 'reference'],
                    help='The method of chimera detection; either denovo or reference. If denovo, Illumina reads will '
                         'first be denoised using the UNOISE algorithm. If reference-based, sequences will be compared '
                         'against the UNITE reference dataset for UCHIME.')

# optional; whether to keep the chimera reads and write out to chimera output directory
parser.add_argument('--keep-chimers',
                    default=settings['chimera_check']['keep_chimeras'],
                    action='store_true',
                    help='Whether to keep the chimera sequences and write them to their own output directory. If '
                         'flag is used, will write output files for chimera sequences. Files are created regardless '
                         'of whether chimeras are detected in a given sample. If flag is not used, will default to '
                         'whatever the setting is in the bioinformatics configuration file.')

# parse command line options and defaults into a dictionary
args = vars(parser.parse_args())


#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# the default input path for Illumina sequences must be from the denoised sequence output
if args['input'] is None:
    args.update({'input': fpm['pipeline-output']['denoised'] / f'{DENOISE_PREFIX}_{run_name}'})
else:
    pass

# check if there are Illumina reads to check for chimeras
is_input, illumina_files = check_for_input(args['input'], config_dict=settings, seq_platform=platform)

if is_input:
    check_chimeras(input_files=illumina_files, file_map=fpm, method=args['method'], keep_chimeras=args['keep_chimers'])
else:
    pass

#####################
# PACBIO ############
#####################
platform = 'pacbio'

# the default input path for Illumina sequences must be from the denoised sequence output
if args['input'] is None:
    args.update({'input': fpm['pipeline-output']['separated-subregions'] / f'{ITSX_PREFIX}_{run_name}'})
else:
    pass

# check if there are PacBio reads to check for chimeras
is_input, pacbio_files = check_for_input(args['input'], config_dict=settings, seq_platform=platform)

if is_input:
    check_chimeras(input_files=pacbio_files, file_map=fpm, method=args['method'])
else:
    pass

#####################
# SANGER ############
#####################
platform = 'sanger'

# the default input path for Sanger sequences must be from the denoised sequence output (?)
# if args['input'] is None:
#     args.update({'input': fpm['pipeline-output']['denoised'] / f'{DENOISE_PREFIX}_{run_name}'})
# else:
#     pass

# check if there are PacBio reads to check for chimeras
is_input, sanger_files = check_for_input(args['input'], config_dict=settings, seq_platform=platform)

if is_input:
    check_chimeras(input_files=sanger_files, file_map=fpm, ref=args['method'])
else:
    pass

# continue to next
continue_to_next(__file__, settings)