from mapping import filepath_map as fpm

import argparse, pathlib
from pathlib import Path
from climush.constants import *
from climush.bioinfo import denoise
from climush.utilities import *

# import settings, get bioinformatics run name from settings
settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

# ADD COMMAND LINE OPTIONS

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Error correction of Illumina reads using the DENOISE3 algorithm from '
                                             'vsearch.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# optional; input directory
parser.add_argument('-i', '--input',
                    default=fpm['pipeline-output']['derep-full-length'] / f'{DEREP_PREFIX}01_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to a directory containing Illumina sequence files to denoise. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# optional; output directory
parser.add_argument('-o', '--output',
                    default=fpm['pipeline-output']['denoised'] / f'{DENOISE_PREFIX}_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to write denoised Illumina sequences to. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# optional; minimum read count
parser.add_argument('--minsize',
                    default=settings['denoise']['minsize'],
                    type=int,
                    help='The minimum number of identical reads that are required for a centroid sequence to '
                         'be written to the output file.')

# optional; minimum abundance skew
parser.add_argument('--alpha',
                    default=settings['denoise']['alpha'],
                    type=float,
                    help='The minimum abundance skew; the factor by which the centroid sequence read count must '
                         'exceed the target sequence in order for the target to be clustered. The default value '
                         'is 2.0.')

# optional; clustering threshold
parser.add_argument('--clust_threshold',
                    default=settings['denoise']['clust_threshold'],
                    type=float,
                    help='The minimum pairwise alignment identity required to cluster the target sequence to a '
                         'centroid sequence. Must be a value between 0.0 and 1.0, inclusive.')

# parse the command line options, including default settings, into a dictionary
args = vars(parser.parse_args())

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

is_input, illumina_files = check_for_input(
        args['input'],
        config_dict=settings,
        file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
    )

if is_input:
    denoise(input_files=illumina_files, file_map=fpm,
            alpha=args['alpha'], minsize=args['minsize'], clust_threshold=args['clust_threshold'])
else:
    pass


#####################
# PACBIO ############
#####################
# platform = 'pacbio'

# is_input, pacbio_files = check_for_input(
#         args['input'],
#         config_dict=settings,
#         file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
#     )
#
# if is_input:
#     denoise(input_files=pacbio_files, file_map=fpm,
#             alpha=args['alpha'], minsize=args['minsize'], clust_threshold=args['clust_threshold'])
# else:
#     pass

#####################
# SANGER ############
#####################
# platform = 'sanger'

# is_input, sanger_files = check_for_input(
#         args['input'],
#         config_dict=settings,
#         file_identifier=[*SEQ_FILE_PREFIX_DICT[platform], platform]
#     )
#
# if is_input:
#     denoise(input_files=sanger_files, file_map=fpm,
#             alpha=args['alpha'], minsize=args['minsize'], clust_threshold=args['clust_threshold'])
# else:
#     pass


# continue to the next script in the pipeline
continue_to_next(__file__, settings)