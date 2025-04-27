import argparse, pathlib, warnings, re
from pathlib import Path

## SCRIPT-WIDE CONSTANTS ###############################################################################################

# if true, the script will use a set of pre-defined command line parameters rather than accept values provided via
#   the command line when running this script
# this is intended for use within the IDE for testing
TEST_MODE = False

########################################################################################################################


## COMMAND LINE OPTIONS ################################################################################################

if TEST_MODE:

    # warn user that the test mode version of this script is being used
    warnings.warn(
        message='This script is running in TEST_MODE, meaning command line options are not being used.',
        category=UserWarning,
    )

    # create a dictionary of pre-set command line options rather than ingesting from command line
    args = {
        'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/02_prefiltered_test/'),
        # 'output': None,
        'wrong_name': 'psf2502_2504',
        'correct_name': 'isp-2504',
    }

else:

    ## INSTANTIATE ARGUMENT PARSER ##
    parser = argparse.ArgumentParser(
        prog=Path(__file__).stem,
        description='Rename bioinformatics output file names with a different bioinformatics run name file/directory '
                    'tab for cases where the run_name value in the bioinformatics configuration file was not updated '
                    'prior to executing pipeline step on the input sequences',
        epilog='This script is part of the CliMush bioinformatics pipeline.',
    )

    ## FILE PATHS ##

    # [REQUIRED] input path to a directory or file that has the incorrect bioinformatics run name
    parser.add_argument(
        '-i', '--input',
        type=pathlib.PosixPath,
        required=True,
        help='The path to the directory containing the files that are incorrectly tagged with the wrong bioinformatics '
             'run name. This directory should be the parent directory of any files or directories that require the '
             'replacement of the incorrect bioinformatics run name with the correct version. This directory name '
             'itself will not be updated.'
    )

    # (optional) output path to the renamed file
    # parser.add_argument(
    #     '-o', '--output',
    #     type=pathlib.PosixPath,
    #     required=False,
    #     default=None,
    #     help='The path to write the corrected file and directory names from the input path to; if nothing provided, it '
    #          'will use the input path (-i / --input) as the output path.'
    # )

    ## FILE / DIRECTORY TAGS ##

    # [REQUIRED] the incorrect run_name tag that is on the files and directories in the input path
    parser.add_argument(
        '--wrong-name',
        type=str,
        required=True,
        help='The incorrect run_name that was originally used to produce the input files and directories in the '
             'input path. This string will be used to search for all files and directories within the input directory '
             'that match this incorrect run_name so that it can be replaced by the correct run_name provided to '
             'this script with the argument --correct-name.'
    )

    # [REQUIRED] the correct run_name tag that will be used to replace the incorrect run_name tag by this script
    parser.add_argument(
        '--correct-name',
        type=str,
        required=True,
        help='The correct run_name to use as the replacement run_name in the files and directories within the input '
             'file path each time the --wrong-name string is encountered.'
    )

    ## PARSE ARGUMENTS INTO DICTIONARY ##

    args = vars(parser.parse_args())

########################################################################################################################


## FIND INCORRECTLY NAMED FILES AND DIRECTORIES ########################################################################

# create an empty list to append the paths of incorrectly named files and directories to
incorrect_paths = []

# keep track of how many files and directories were searched
file_search_count = 0
dir_search_count = 0

# create a simple function that will search and append to the list
def search_and_append(file_or_dir, file_or_dir_root, search_str, match_list):
    '''
    Conduct search for string and add path to input list if found

    :param file_or_dir: a file or directory name, as returned to the filename or dirname
    variable during Path.walk()
    :param file_or_dir_root: the parent directory of the input file or name, as returned
    to the root variable during Path.walk()
    :param search_str: the string to search for in the file or directory name, typically
    the command line argument for --wrong-name
    :param match_list: the list to append the incorrect file path to when a match to the
    search_str is located
    :return: None; update the input list, match_list, if a search result is found, otherwise,
    do nothing
    '''

    if re.search(search_str, file_or_dir, re.I):
        match_list.append(file_or_dir_root / file_or_dir)
    else:
        pass

    return None

# walk from the top of the directories (input path) to bottom looking for incorrectly named files / directories
for root, dirnames, filenames in args['input'].walk(top_down=True):

    ## FILES ##

    # search for the wrong run_bame in the files in this directory
    for file in filenames:
        file_search_count += 1  # add to counter to keep track of the number of files searched
        search_and_append(
            file,
            root,
            args['wrong_name'],
            incorrect_paths,
        )

    ## DIRECTORIES ##

    # search for the wrong run_name in the files in this directory
    for dir in dirnames:
        dir_search_count += 1  # add to counter to keep track of the number of directories searched
        search_and_append(
            dir,
            root,
            args['wrong_name'],
            incorrect_paths,
        )

print(f'The number of files and directories explored during the path walk of the input directory:\n'
      f'   {args["input"]}\n'
      f'   files searched       = {file_search_count}\n'
      f'   directories searched = {dir_search_count}\n'
      f'total matching paths    = {len(incorrect_paths)}\n')

########################################################################################################################


## CREATE CORRECTLY NAMED FILES AND DIRECTORIES ########################################################################

# add the wrong file/directory paths to the dictionary as a key, and add the respective corrected name to this
#  dictionary as its value
renaming_dict = {
    wrong_path:'' for wrong_path in incorrect_paths
}


# iterate through each path that has the incorrect run_name tagged to it...
for wrong_path in renaming_dict:

    # replace the wrong part of the name with the correct version of the run_name to create a corrected path
    corrected_path = Path(re.sub(args['wrong_name'], args['correct_name'], str(wrong_path)))

    # add this corrected path to the dictionary as a value to the wrong path
    renaming_dict.update({wrong_path:corrected_path})

########################################################################################################################


## RENAME INCORRECTLY NAMED FILES AND DIRECTORIES WITH CORRECTED NAMES #################################################

for wrong_path, corrected_path in renaming_dict.items():

    # rename old with new
    wrong_path.replace(corrected_path)

########################################################################################################################


# ## CHECK INSIDE FILES FOR WRONG NAME AND REPLACE #######################################################################
#
# for corrected_path in renaming_dict.values():
#
#     if corrected_path.is_file():
#
#         fixed_text = []
#
#         with open(corrected_path,'rt') as f_in:
#             for line in f_in.readlines():
#                 fixed_text.append(re.sub(args['wrong_name'], args['correct_name'], line))
#
#         corrected_filename = (corrected_path.parent / (corrected_path.stem + '_corrected')).with_suffix(corrected_path.suffix)
#         with open(corrected_filename, 'wt') as f_out:
#             f_out.writelines(fixed_text)
#
#     else:
#         continue
#
# ########################################################################################################################
