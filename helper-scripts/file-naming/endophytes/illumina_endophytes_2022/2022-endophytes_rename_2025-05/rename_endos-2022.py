from pathlib import Path
import argparse, pathlib, re
from climush.utilities import mkdir_exist_ok, add_prefix
from climush.constants import SEQ_FILE_GLOB

parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='rename 2022 endophytes with illumina_ prefix',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)

## FILE PATHS ##

# REQUIRED; input file path to the parent directory containing the sequence files to quality filter
parser.add_argument(
    '-i', '--input',
    default=None,
    required=True,
    type=pathlib.PosixPath,
    action='store',
    help='',
)

parser.add_argument(
    '--old-prefix',
    default='its1',
    required=False,
    type=str,
    action='store',
    help='',
)

parser.add_argument(
    '--new-prefix',
    default='illumina',
    required=False,
    type=str,
    action='store',
    help='',
)

args = vars(parser.parse_args())


# go through each of the sequence files in the input directory
for input_file in args['input'].glob(SEQ_FILE_GLOB):

    # if the file name starts with the prefix you want to replace...
    if input_file.name.startswith(args['old_prefix']):

        # replace this prefix with the new prefix
        add_prefix(
            file_path=input_file,
            prefix=args['new_prefix'],
            dest_dir=args['input'],
            action='rename',
            f_delim='_',
            output_compressed=False,
            replace_prefix=True,
        )

    # if the file name already starts with the prefix you want to use...
    elif input_file.name.startswith(args['new_prefix']):

        # do nothing here, pass over this file without any prefix substitution or addition
        pass

    # if the file name does not start with either the old prefix or the new prefix, assume it has no prefix
    else:

        # rename this file with the new prefix without replacing an existing file name prefix
        add_prefix(
            file_path=input_file,
            prefix=args['new_prefix'],
            dest_dir=args['input'],
            action='rename',
            f_delim='_',
            output_compressed=False,
            replace_prefix=False,
        )

# check the renamed files and check that the R1 / R2 file tags are preceded by a _ not -, and replace - if -
rename_dict = {}

for rename01_file in args['input'].glob(SEQ_FILE_GLOB):

    if re.search(r'.+?(?=-R1)', rename01_file.name, re.I):

        new_filename = re.sub(r'-R1', r'_R1', rename01_file.name)
        new_filepath = rename01_file.parent / new_filename
        rename_dict.update({rename01_file: new_filepath})

    elif re.search(r'.+?(?=-R2)', rename01_file.name, re.I):

        new_filename = re.sub(r'-R2', r'_R2', rename01_file.name)
        new_filepath = rename01_file.parent / new_filename
        rename_dict.update({rename01_file: new_filepath})

    else:
        continue


# go through adn rename what needs renaming
for old_name, new_name in rename_dict.items():
    old_name.rename(new_name)