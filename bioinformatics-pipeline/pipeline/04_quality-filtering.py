from mapping import filepath_map as fpm

import argparse, sys, subprocess, pathlib
from datetime import datetime
from pathlib import Path
from climush.constants import *
from climush.bioinfo import merge_reads, quality_filter
from climush.utilities import *

settings = get_settings(fpm)
run_name = settings['run_details']['run_name']

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Quality filter reads using VSEARCH.',
                                 epilog='')

# input directory containing the files to quality filter
parser.add_argument('-i', '--input',
                    default=fpm['pipeline-output']['primers-trimmed'] / f'trim_{run_name}',
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to quality filter. If nothing provided, '
                         'will default to the location that is expected in the Docker container\'s native file '
                         'structure, detailed in pipeline/mapping.py.')

# if flag is included, override config settings and use this
parser.add_argument('--merge', action=argparse.BooleanOptionalAction,
                    default=settings['quality_filtering']['illumina']['merge_reads'],
                    help='If you want to override the settings in the configuration file, then include the flag for '
                         '--no-merge if you do not want to merge reads, and --merge if you want to merge reads. '
                         'Otherwise, the setting in the configuration file will be used (default is to merge).')

# to run quality filtering and merging separately, for time restriction
parser.add_argument('--merge-from', nargs='?',
                    default=None,
                    type=pathlib.PosixPath,
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

#####################
# ILLUMINA ##########
#####################
platform = 'illumina'

# if a merge_from path provided, ensure it is a Path object
if args['merge_from'] is None:
    merge_from = None
else:
    merge_from = Path(args['merge_from'])

# FIGURE OUT HOW TO UPDATE THE SETTINGS WITH COMMAND LINE INPUT
print(f'Currently, this is not updated to accomodate changes to the configuration based on command line arguments. '
      f'Working on it... \n')

# look for input in the input path that are illumina reads
is_input, illumina_files = check_for_input(input_path, seq_platform=platform)

if is_input:
    # if the --merge-from flag is not included, go through quality filter process
    if merge_from is None:
        filtered_path = quality_filter(input_files=illumina_files, platform=platform, file_map=fpm)
        is_merge, filtered_files = check_for_input(filtered_path, seq_platform=platform)
    else:
        is_merge, filtered_files = check_for_input(filtered_path, seq_platform=platform)
    # confirm that there are
    if args['merge'] and is_merge:
        print(f'\nMerging reads from {filtered_path}...')
        merge_reads(input_files=filtered_files, file_map=fpm)
    elif args['merge'] and not is_merge:
        print(f'Could not locate files for merging in the path {filtered_path}. Please check the output from '
              f'quality filtering to confirm that {platform.capitalize()} reads are present.\n')
    else:
        print(f'\nNot merging reads; instead using forward reads only.')
else:
    pass
end_time = datetime.now()
run_time = end_time - start_time
print(f'{int(len(filtered_files)/2)} paired-end samples were run through quality filtering and merging '
      f'in {run_time}.\n')
#####################
# PACBIO ############
#####################
platform = 'pacbio'

# pacbio_dir = PIPE_OUT_MAIN / Path(f'03_remove-primers/trim_{run_name}')
#
# # check that there are Illumina reads to pre-filter
# if input_files_present(file_path = pacbio_dir):
#     print(f'PacBio reads were detected in {pacbio_dir.parent}. '
#           f'Processing {count_files(pacbio_dir, search_for=SEQ_FILE_GLOB)} samples...\n')
#


#####################
# SANGER ############
#####################
platform = 'sanger'

# when all are prefiltered, continue to next
continue_to_next(__file__, settings)