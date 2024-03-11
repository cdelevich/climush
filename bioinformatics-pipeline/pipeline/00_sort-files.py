from mapping import filepath_map as fpm

import argparse
import sys
import re
##REMOVE AFTER PACKAGE TESTING#######
from pathlib import Path
sys.path.insert(0, str(fpm['package']['main']))
#####################################
from climush.constants import *
from climush.utilities import *

config_dict = import_config_as_dict(fpm['config']['main'], file_handle='pipeline-settings', config_section='all')

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Identifies sequencing platform of sequencing files '
                                             'and sorts them into their respective pipeline paths. '
                                             'Checks file names for correct naming convention. Confirms '
                                             'that sequences have been demultiplexed, and if not, sends '
                                             'sequences to be demultiplexed.',
                                 epilog='')

parser.add_argument('-sp', '--seq-path', default=fpm['sequences']['main'],
                    help='The path to the sequencing files. Will default to the location that is '
                         'expected with the Docker container\'s native file structure.')

parser.add_argument('-cp', '--config-path', default=fpm['config']['main'],
                    help='The path to the configuration files. Will default to the location that is '
                         'expected with the Docker container\'s native file structure.')

args = vars(parser.parse_args())

seq_path = args['seq_path']
config_path = args['config_path']

## REMOVE AFTER TESTING
# seq_path = fpm['sequences']['main']
# config_path = fpm['config']['main']
######################

sort_input_files(filepath_dict = fpm)

# check which (if any) directories need action before moving to demultiplexing (if necessary)
dir_requires_action = [dir for dir in seq_path.glob('*') if re.search(NEEDS_ACTION_REGEX, dir.stem)]
for dir in dir_requires_action:
    rename_dir = fpm['sequences']['rename']
    if dir.stem == rename_dir.stem:
        print(f'yes, right path')
        num_to_rename = count_files(dir)
        auto_continue = input(f'\n{num_to_rename} sequencing files in the sequencing file set '
                              f'require renaming to match the current file naming convention. If you '
                              f'wish to automatically rename these files, first make sure that '
                              f'there is a file renaming configuration file in the {dir.stem} folder. Do you '
                              f'want to rename these files automatically?[Y/N]\n')
        if re.search(AFFIRM_REGEX, auto_continue, re.I):
            subprocess.run(['python3', fpm['pipeline']['rename']])
        else:
            platform_type = input(f'At minimum, the sequencing files are required to have the sequencing platform '
                                  f'associated with each file. What type of sequencing technology was used to '
                                  f'produce these sequences?[illumina,pacbio,sanger]')
            if re.search('^i', platform_type, re.I):
                add_prefix(dir, prefix='illumina', dest_dir=dir)
            elif re.search('^p', platform_type, re.I):
                add_prefix(dir, prefix='pacbio', dest_dir=dir)
            elif re.search('^s', platform_type, re.I):
                add_prefix(dir, prefix='sanger', dest_dir=dir)
            else:
                print(f'The prefix provided, {platform_type}, does not match any of the available options.\n')
    else:
        num_to_check = count_files(dir)
        print(f'{num_to_check} sequencing files could not be sorted due to their unfamiliar '
              f'naming convention. Please see the directory {dir.stem} to rename and sort '
              f'these files manually. Files may be able to be renamed using the script '
              f'rename_sequence_files.py if a configuration file is provided. If you continue '
              f'in the pipeline without renaming these files, they will be excluded from '
              f'downstream processes.\n')
        platform_type = input(f'At minimum, the sequencing files are required to have the sequencing platform '
                              f'associated with each file. What type of sequencing technology was used to '
                              f'produce these sequences?[illumina,pacbio,sanger]')
        if re.search('^i', platform_type, re.I):
            add_prefix(dir, prefix='illumina', dest_dir=dir)
        elif re.search('^p', platform_type, re.I):
            add_prefix(dir, prefix='pacbio', dest_dir=dir)
        elif re.search('^s', platform_type, re.I):
            add_prefix(dir, prefix='sanger', dest_dir=dir)
        else:
            print(f'The prefix provided, {platform_type}, does not match any of the available options.\n')

        # move back to main sequence folder, then re-sort
        for file in dir.glob('*'):
            move_file(source_path=file, dest_path=mkdir_exist_ok(new_dir=fpm['sequences'][platform_type]))
        filter_empty_files(fpm['sequences']['main'])

# initialize output directory here
mkdir_exist_ok(fpm['pipeline-output']['main'])

# go to next step of the pipeline (demultiplexing) if no other step was chosen from the CLI
continue_to_next(Path(__file__), config_dict)

