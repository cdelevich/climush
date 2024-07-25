import argparse, tomlkit, itertools, os, shutil, re, sys
from pathlib import Path
import pandas as pd
from datetime import datetime
from climush.utilities import prompt_yes_no_quit, is_pathclass, import_filepath, flag_multiple_files, get_settings, mkdir_exist_ok, print_indented_list, prompt_print_options

from climush.constants import CONFIG_FILETYPE, SEQ_FILE_RE, GZIP_REGEX, GZIP_GLOB, MOCK_COMM_RE, NEG_CTRL_RE, UNDET_RE

# start timer for calculating run time
start_time = datetime.now()

# indicate whether to employ 'test mode', used for running interactively in IDE
test_mode = False

# should do this programatically but I'm not
run_date = '2023-10'

# set a variable that will determine if a warning was triggered at some point throughout the script
no_warnings = True

##########################################################
###### ADD TO CLIMUSH PACKAGE ##########################################################################################
##########################################################
print('section 01: add to climush package\n')
# function that will update the neon_domains dictionary if there is a new value for the domain
def update_domain_dict(domain_dict, configuration_dict, site_section='site_name'):
    '''
    Updates the NEON domains dictionary with new site values.

    Looks through the site names given by the configuration
    file submitted by the user and checks these values against
    the values in the neon_domains dictionary. If a site name
    provided in the initialization file is not already a value
    in the neon_domains dictionary, it will add it to the
    dictionary to the corresponding domain.
    :param domain_dict: a nested dictionary where the primary
    key is the state abbreviation, the secondary key is the
    NEON domain, and the value of the NEON domain is an empty
    set.
    :param site_section: the name of the section that contains
    information on the site name used in the old file names;
    defaults to 'site_name', which is the name given in the
    original configuration file template, but can be changed
    if altered in the configuration file.
    :param configuration_dict: dictionary read in from the user-
    provided configuration file; by default, it should be called
    config_dict, but able to change if not.
    :return: None; will update the input dictionary.
    '''

    site_dict = configuration_dict[site_section]

    for state, site in site_dict.items():

        # get the associated NEON domain for this state
        domain = list(domain_dict[state].keys())[0]

        # add the site name to the NEON domain for this state
        if isinstance(site, list):
            for s in site:
                domain_dict[state][domain].add(s)
        else:
            domain_dict[state][domain].add(site)

    return None

# create a dictionary of the original file names and their sizes, as a check after renaming
def file_size_dict(directory, filesize_unit='MB'):
    '''
    Create a dictionary of original file names and their
    corresponding file sizes.

    For each of the original file names, the size of the file
    is calculated and saved into a dictionary. This will serve
    as a test after renaming files that the renamed files are
    corresponding to the correct original files.
    :param directory: the path to the directory containing
    the sequence files to rename; can be absolute or relative to
    the location of this script.
    :param filesize_unit: the memory unit to use for the output,
    default is megabytes (MB); options: bytes (B), kilobytes (KB),
    megabytes (MB), and gigabytes (GB)
    :return: dictionary where keys are names of the original files
    and values are the size of the original files in (bites?)
    '''

    # get path of all original sequence files in the input directory
    filepaths = [file for file in directory.glob(GZIP_GLOB)]

    # get the conversion to use for output file size
    file_size_conversion = {'B': 0, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3}  # conversion relative to bytes

    # only allow filesize_unit input value that is in the file size conversion dictionary
    accepted_unit_regex = '|'.join(list(file_size_conversion.keys()))

    # confirm that file size unit input value is valid; prompt user if not
    if re.search(accepted_unit_regex, filesize_unit, re.I):
        pass
    else:
        print(f'ERROR. The provided file size unit, {filesize_unit}, is not a valid unit. Please chose a valid '
              f'unit from the following options:\n')
        filesize_unit = prompt_print_options(option_list = list(file_size_conversion.keys()))

    # get file conversion and file size for each file in the input dictionary
    conversion = file_size_conversion[filesize_unit]  # conversion value from bytes

    filesize_dict = {file.name: {'file_size': (int(file.stat().st_size) / conversion),
                                 'unit': (filesize_unit.upper())}
                     for file in filepaths}

    return filesize_dict

# before renaming files, create a copy of all files to be renamed and send to zipped subdirectory (for safety)
def copy_original_files(directory, copy_directory='original_files', compress=True):
    '''
    Creates a copy of original sequence files to zipped subdirectory.

    Copies all sequence files in the provided directory
    :param directory: the path to the directory containing the
    sequence files to be renamed.
    :param copy_directory: name or path of the subdirectory to
    copy the original sequence files to; default is a directory
    called 'original_files' within the provided directory
    :param compress: whether to compress the folder containing the
    copies of the original sequence files; default is True and will
    zip the copy directory.
    :return: print confirmation that all files were successfully
    copied to the copy directory
    '''

    # convert input to a Path class object
    if is_pathclass(directory, exit_if_false=False):
        dir_path = directory
    else:
        dir_path = Path(directory)

    # create the output directory (copy) in same directory as the input directory
    if is_pathclass(copy_directory, exit_if_false=False):
        copy_path = copy_directory
    else:
        copy_path = dir_path / copy_directory

    # create the copy directory

    # if the directory already exists, prompt user for input
    if copy_path.is_dir():

        # prompt for user input
        msg = f'A directory with the name {copy_path.name} already exists in: \n'\
              f'\t{dir_path}\n'\
              f'Do you want to delete this directory and continue? If you do not want to delete this directory, '\
              f'then the program will exit.'

        # if response is yes, will continue; if no/quit, will exit here
        prompt_yes_no_quit(message=msg)

        # remove the pre-existing directory that matches the copy directory name
        shutil.rmtree(copy_path)

    # make the new copy path directory
    mkdir_exist_ok(new_dir = copy_path)

    # create a list of files to copy, excluding any non-sequencing files
    files_to_copy = [ file for file in dir_path.glob(GZIP_GLOB) ]

    # copy each sequence file to the new copy directory
    for file in files_to_copy:
        shutil.copy2(file, copy_path)  # to preserve extra metadata, use copy2

    # before archiving copy, get count of sequence files from source and destination
    num_copied_files = len(list(copy_path.glob('*')))
    num_original_files = len(files_to_copy)

    # confirm all files were copied to the copy directory
    if num_copied_files == num_original_files:

        # if all files copied, then compress output
        if compress:
            # zip the copy directory
            shutil.make_archive(copy_path, 'zip', copy_path)

            # remove the unzipped copy directory
            shutil.rmtree(copy_path)

        return print(f'\nA copy of {num_copied_files} original sequence files was successfully made in:\n'
                     f'\t{copy_path}')

    else:
        num_missing = abs(num_copied_files - num_original_files)
        return print(f'ERROR: {num_missing}/{num_original_files} sequence files from the source directory:\n'
                     f'\t{dir_path}\n'
                     f'were not copied to the destination directory:\n'
                     f'\t{copy_path}')

# convert the glob format provided from configuration file to a regex format
def convert_glob_to_regex(glob_str):
    '''
    Convert simple text with asterisk wildcard to regex.

    Config file accepts wildcards in the common form of an
    asterisk (*) to make it more user-friendly. This function
    will take that input and convert to a regular expression.
    If there are multiple strings in the input for this function,
    it will return a list of regex strings, which can be later
    joined to kept separate.
    :param config_input: string or list of strings where asterisk
    wildcard is used; wildcard not necessary, but if not present, will
    build regex for exact match to input (i.e., not part of other
    word).
    :return: list of regex if input list, single regex string if input
    string.
    '''

    if isinstance(glob_str, str):
        pass
    else:
        print(f'ERROR. Input for convert_glob_to_regex() function in the renaming file must be a string, not '
              f'a list. Please fix the code and then rerun.\n')
        print(f'Exiting...')
        # sys.exit()

    glob_str = glob_str.replace(r'.', r'\.')  # if there's a '.' as part of original label, need to escape for regex

    if glob_str.startswith('*'):  # if wildcard at start, search for string at end of label
        result = glob_str.replace(r'*', r'.*') + '$'
    elif glob_str.endswith('*'):  # if wildcard at end, search for string at beginning of label
        result = '^' + glob_str.replace(r'*', r'.*')  # last asterisk, can be anything or nothing after
    else:  # otherwise, search for exact match
        result = '^' + glob_str + '$'

    return result

# extracted nested portion of a nested dictionary
# def extract_nested(nested_dict, depth):
#     # create a depth counter that will keep track of the number of key levels that it iterates through
#     depth_count = 0
#
#     # search through the depths of the dictionary until the desired depth is reached
#     while depth_count < depth:
#         for key, val in nested_dict.items():
#             if isinstance(val, dict):
#                 depth_count += 1
#             else:
#                 if depth_count == 0:  # if the very first set of keys/values explored is determined to not be nested
#                     return print('WARNING. Provided dictionary does not appear to be nested.\n')  # alert user
#
#         # implement a safety check, to prevent any infinite loops
#         if depth_count > 50:
#             return print('WARNING. There may be an error with the function. 50 levels of the dictionary were searched '
#                          'and did not return a result.\n')
#     else:
#         output_dict = {key: val}
#         return output_dict

def extract_nested(nested_dict):
    output_dict = {}

    for key, val in nested_dict.items():
        if isinstance(val, dict):
            output_dict.update(val)  # val itself is a dictionary, so only add this sub-dict to the output dictionary
        else:
            print('WARNING. Provided dictionary does not appear to be nested.\n')

    return output_dict

# NOT ACTUALLY USED IN THIS SCRIPT, BUT MIGHT BE HELPFUL ELSEWHERE?
# function to flatten the values in a dictionary to a single list, recursively
def flatten_dict_values(input_dict, flattened_list=[]):
    '''
    Return values from a simple dictionary (non-nested) as a list of values.
    If values are a mix of lists and/or sets of values, the items within these
    lists/sets will be unpacked.
    :param input_dict: a simple dictionary (non-nested); if using a nested dictionary,
    then subset the dictionary so that it appears non-nested.
    :param flattened_list: list to add the flattened values to; defaults to an empty
    list at the start, but will be a partially filled list upon recursion
    :return: list of values from the input dictionary; always a list, even if input
    values are a set, so that you can observe any values that belong to multiple keys
    '''

    # given this function is recursive, input may not be dict upon recursion

    # if it is the initial input, it will be a dict; get values from the input dict
    if isinstance(input_dict, dict):
        flattened_list = []  # list has to be reset to empty or vals will be added each time it is used
        vals_to_flatten = input_dict.values()
    # if recursed, it will be a set or list; rename to match the object name of dict input
    else:
        vals_to_flatten = input_dict

    # go through each of the values in the input dictionary or recursion set/list
    for val in vals_to_flatten:
        # if the value is a list or set, then recurse to flatten (also added tuple, though I don't use)
        if isinstance(val, list) or isinstance(val, set) or isinstance(val, tuple):
            flatten_dict_values(val, flattened_list=flattened_list)
        # if the value is a dictionary, then the input dictionary is nested, which is not something this function handles
        elif isinstance(val, dict):
            print(f'WARNING. The input dictionary is nested, so it is ambiguous how to flatten the values '
                  f'in this dictionary. Subset the dictionary so that it is non-nested, and then '
                  f'rerun.\n')
        # if not a list or set, then add to flattened list (should mean its a string, float, int, etc.; hopefully)
        else:
            flattened_list.append(val)

    return flattened_list

# search for a value in a dict where values are sets/lists/tuples, and return its key when a match is found
def search_conversion_dict(input_dict, search_val):
    '''
    Search the input dictionary values for a match to the provided string, looking into
    lists, sets, and/or tuples, and return the key corresponding to this value once it
    is located.
    :param input_dict: dictionary to search the values in; should be non-nested
    :param search_val: the value to search for in the dictionary, as a string, float,
    integer, etc. (not a collection-type object like a list)
    :return: the key corresponding to the input search value
    '''

    # look through each key/value pair in the dictionary
    for key, val in input_dict.items():

        # if the value is a collection-type object, then flatten it to get a single value
        if isinstance(val, list) or isinstance(val, set) or isinstance(val, tuple):
            if search_val in val:
                return key
        else:
            if search_val == val:
                return key

    print(f'WARNING. The input search string, {search_val}, was not located as a value in the '
          f'input dictionary.\n')
    return None

########################################################################################################################

# if using interactively in an IDE, set test_mode = True to bypass command line arguments
if test_mode:
    # define path to this script (only works in certain contexts, i.e., my local file system)
    current_script = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/rename_sequence_files.py')

    # confirm that you want to continue before proceeding; will also show if running in terminal if still set to True by accident
    msg = (f'Test mode is currently engaged for {current_script.name}. Command line arguments will not be accessible. '
           f'Would you like to continue?')
    prompt_yes_no_quit(message=msg)

    args = {'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/illumina_soil-litter_2023-10_empty-files'),
            'config': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/climush-bioinfo_file-rename.toml'),
            'output': None,  # use default
            'no_track': False,  # produce the output conversion table
            'no_copy': False}  # save a copy of the original files

# if not using interactively (i.e., running from terminal, etc.), set test_mode = False to activate command line arguments (default)
else:
    current_script = Path(__file__)
    ##########################################################
    ## DEFINE COMMAND LINE OPTIONS ###########################
    ##########################################################

    # command line arguments
    parser = argparse.ArgumentParser(description='Rename CliMush raw sequence files.')

    # input directory containing the sequence files to rename
    parser.add_argument('-i', '--input',  # positional argument
                        help=f'Path to the directory containing the sequences to rename; either an absolute path or '
                             f'relative to the location of this script, which is located in: {Path(__file__).parent}')

    # path to the configuration file to use for renaming the input files
    parser.add_argument('--config',  # optional
                        help = f'Relative or absolute path to the .toml configuration file for this set of sequence files, '
                               f'which should be updated by the user prior to running this script. If a path is not '
                               f'provided, it will look for a file with the extension \'.toml\' in the same directory '
                               f'as the input path provided.')

    # path to use for the renamed output files
    parser.add_argument('--output',  # optional
                        help = f'Relative or absolute path of the directory to contain the renamed sequencing files. By '
                               f'default, this folder will be named after the sequencing run with '
                               f'the prefix \'_raw-reads\'.')

    # flag option to output a conversion table
    parser.add_argument('--no-track',  # optional
                        action = 'store_true',
                        help = f'If this flag is used, a conversion table will be *not* be saved that showing the old file names '
                               f'and their corresponding updated file name. It is advised to always include a conversion table '
                               f'and not use this flag.')

    # flag option to NOT save a compressed copy of the original files
    parser.add_argument('--no-copy',  # optional
                        action = 'store_true',
                        help = f'If this flag is used, no copy of the original file names, compressed as a zip file for '
                               f'the collection, will be saved. This is the safe option if there is not already a copy '
                               f'of the original files elsewhere. However, it does take additional space and time if a '
                               f'compressed copy of the original files is saved.')

    # parse command line arguments into a dictionary
    args = vars(parser.parse_args())  # parse here to facilitate use of 'test_mode'

##########################################################
## PARSE COMMAND LINE ARGUMENTS ##########################
##########################################################
print('section 02: parse command line arguments\n')
## INPUT PATH
# ensure that the input path is a Path object
input_path = import_filepath(args['input'], must_exist=True)


## CONFIGURATION FILE PATH
# check if a configuration filepath was provided
if args['config'] is None:
    # if a config path file was not provided, then search for a config file in the input directory
    config_path = flag_multiple_files(file_path=input_path, search_for=f'*{CONFIG_FORMAT}')
else:
    # if a config path was provided, then ensure the file exists and is a path object
    config_path = import_filepath(args['config'], must_exist=True)


## OUTPUT FILE PATH FOR RENAMED SEQUENCES
# if an output path is not provided...
if args['output'] is None:
    # leave as empty string for now; need to gather info on sequencing run before naming
    rename_path = input_path.parent / f'illumina_soil-litter_{run_date}_file-rename-conversion.csv'
# if an output path is provided...
else:
    # then make sure it is a Path object
    rename_path = Path(args['output'])


## OPTION TO NOT MAKE COPY OF ORIGINAL FILES
# ensure that the --no-copy setting used is what is desired
if args['no_copy']:
    msg = (f'A copy of the files with their original file names will *not* be created, so it'
           f'is advised that you have your own backup of these files with their original file names.'
           f'is advised that you have your own backup of these files with their original file names.'
           f'Do you wish to continue?')
else:
    msg = (f'A copy of the files with their original file names will be created, and may'
           f'be time intensive. Do you wish to continue?')

# if 'yes', returns None and will continue; if 'no' or 'quit', will exit here
prompt_yes_no_quit(message=msg, auto_respond=False)

# save whatever input was received for --no-copy once confirmed
no_copy = args['no_copy']


## OPTION TO NOT MAKE TABLE OF OLD AND NEW FILE NAMES
# ensure that the --no-track setting used is what is desired
if args['no_track']:
    msg = (f'A conversion table showing the original and renamed sequencing file names will *not* be created. It is '
           f'advised that you allow a conversion table to be created. Would you like to continue without creating a'
           f'conversion table? You will not have a record of what the original sequencing files names will be in '
           f'this case.')
    prompt_yes_no_quit(message=msg, auto_respond=False)

# save whatever input was receipt for --no-track once confirmed
no_track = args['no_track']


##########################################################
## IMPORT CONFIGURATION FILE #############################
##########################################################
print('section 03: import configuration file\n')

# read in the configuration file
# not using get_settings function because it won't work and is too complicated for this simple import
with open(config_path, 'rt') as fin:
    config_dict = tomlkit.parse(fin.read())

# define constants
# SEQ_FILE_REGEX now SEQ_FILE_RE
# MOCK_REGEX now MOCK_COMM_RE
# NTC_REGEX now NEG_CTRL_RE
# UNDET_REGEX now UNDET_RE
VALID_ENTRY_REGEX = r'[^(^NA$)|^(^None$)]'  # CHECK HOW USED, NEED TO ADD TO CONSTANTS.PY?

# get the file delimiter from the configuration file
# DELIMITER now file_delim
file_delim = config_dict['delimiter']['delimiter']


##########################################################
## DEFINE REFERENCE DICTIONARIES #########################
##########################################################

# create dict of NEON domains, their associated states, and room to add unique site names from original labels
neon_domains = {'AK': {'D19': set()},
                'AZ': {'D14': set()},
                'CO': {'D13': set()},
                'FL': {'D03': set()},
                'KS': {'D06': set()},
                'MA': {'D01': set()},
                'MN': {'D05': set()},
                'OR': {'D16': set()}}


##########################################################
## GET INFO FROM ORIGINAL FILENAMES ######################
##########################################################
print('section 04: get info from original filenames\n')

# update the neon_domains dictionary to include the site name in the original file names for each domain
update_domain_dict(neon_domains, configuration_dict=config_dict)

# get tuple of the file paths to the original files, with original file names
original_file_paths = tuple(file for file in input_path.glob(GZIP_GLOB))  # use tuple, is immutable

# get a tuple of the original file names without their paths
original_file_names = tuple(file.name for file in original_file_paths)  # use tuple, is immutable

# get the file sizes of each sequencing file before any changes are made
original_file_sizes = file_size_dict(directory=input_path)

# confirm that the original file names exactly match those in the original file size dictionary
assert tuple(original_file_sizes.keys()) == original_file_names, ('ERROR. The original file names do not match the file '
                                                                  'names in the file size dictionary. The final renamed '
                                                                  'file check will be inaccurate if not fixed.')

# create a copy of the original files in a compressed folder, unless --no-copy flag is used
if no_copy:
    pass
else:
    copy_original_files(directory=input_path, copy_directory=(input_path / f'illumina_soil-litter_{run_date}_original-file-names'))


##########################################################
## GET ORIGINAL LABELS USED FOR EACH COMPONENT ###########
##########################################################
print('section 05: get original labels used for each component\n')
# create a dictionary of the file name components, with the new label values and empty sets for the old labels
file_name_conversion = {'seq_platform': {'sanger': set(),
                                         'illumina': set(),
                                         'pacbio': set()
                                         },
                        'compartment': {'litter': set(),
                                       'soil': set(),
                                       'spore': set(),
                                       'sporocarp-f': set(),
                                       'sporocarp-a': set(),
                                       'leaf-sp01': set(),
                                       'leaf-sp02': set(),
                                       'seed-sp01': set(),
                                       'seed-sp02': set(),
                                       'root-sp01': set(),
                                       'root-sp02': set()
                                       },
                        'ecoregion': {'D01': set(),
                                     'D03': set(),
                                     'D05': set(),
                                     'D06': set(),
                                     'D13': set(),
                                     'D14': set(),
                                     'D16': set(),
                                     'D19': set()
                                     },
                        'treatment': {'UC': set(),
                                     'UO': set(),
                                     'UG': set(),
                                     'BC': set(),
                                     'BO': set(),
                                     'BG': set()
                                     },
                        'subplot': {'01': set(),
                                   '02': set(),
                                   '03': set(),
                                   '04': set(),
                                   '05': set(),
                                   '06': set(),
                                   '07': set(),
                                   '08': set(),
                                   '09': set()
                                   }
                        }

# create a dictionary to store the old labels before sorting by updated values
old_label_dict = {k:set() for k in file_name_conversion.keys()}

# add any file names that do not match the label descriptions in the configuration file
invalid_file_names = set()

# for each component in the new file name convention, find the original labels used
for component in file_name_conversion.keys():

    ## FIND LABEL FOR THIS COMPONENT IN THE ORIGINAL FILE NAME

    # from information provided, where is this file component in the original file name?
    original_position = config_dict['component_positions'][component]

    # if this component is in the original file name, position will be tomlkit.items.Integer class
    if isinstance(original_position, int):

        # then go through each file and get the unique labels for this component
        for original_file in original_file_names:

            # get the original label based on the expected position
            original_label = original_file.split(file_delim)[original_position]

            ## CHECK THAT FORMAT OF THIS LABEL MATCHES WHAT IS SPECIFIED IN CONFIGURATION FILE

            # get the accepted label values for this component from the configuration file
            accepted_labels = config_dict['accepted_labels'][component]

            # convert the glob format from the configuration file to regex

            # if acceptable original labels are a list (i.e., multiple possible values)...
            if isinstance(accepted_labels, list):  # e.g., treatment

                # join together possible values
                # sub glob wildcard (*) for regex wildcard (.) w/ zero or more (*) instances of a wildcard
                accepted_labels_regex = '|'.join([convert_glob_to_regex(l) for l in accepted_labels])

            # if acceptable original labels are a dictionary (i.e., changing values depend on compartment type)...
            elif isinstance(accepted_labels, dict):  # e.g., compartment like soil, litter, etc.

                # gather all possible accepted original labels into a list
                accepted_labels_list = []
                for label in accepted_labels.values():
                    if label == "":
                        continue
                    else:
                        accepted_labels_list.append(label)

                # then join these labels together into one regex with same method as for if labels were list
                accepted_labels_regex = '|'.join([convert_glob_to_regex(l) for l in accepted_labels_list])

            # if acceptable original labels are strings...
            elif isinstance(accepted_labels, str):

                # if there is no label in original file name, it will be an empty string; therefore, pass over this
                if len(accepted_labels) == 0:
                    if component == 'ecoregion':
                        accepted_labels_list = []

                        for lab in config_dict['site_name'].values():
                            if isinstance(lab, str):
                                accepted_labels_list.append(lab)
                            elif isinstance(lab, list):
                                accepted_labels_list = accepted_labels_list + lab
                            else:
                                print(f'Type of data entry not recognized.\n')

                        accepted_labels_regex = '|'.join(accepted_labels_list)

                # if there is an acceptable label for this file name component...
                else:

                    # then convert glob format to regex for this string
                    accepted_labels_regex = convert_glob_to_regex(accepted_labels)

            # if the acceptable original label is not a string, dictionary, or list...
            else:

                # print an error message and exit the script
                print(f'ERROR. The accepted label format of the original file name for '
                      f'{component}, provided in the configuration file, is not an accepted data type. Accepted '
                      f'data types for this section of the configuration file are strings, dictionaries, and lists.'
                      f'The provided data type in the configuration file for this component is {type(accepted_labels)}.'
                      f'Please fix this issue in the configuration file, and then rerun the renaming script.\n')
                print(f'Exiting...\n')
                # sys.exit(1)

            # if this is in fact an accepted label, add to old label dictionary

            if re.search(accepted_labels_regex, original_label, re.I):

                # add value(s) to the old dictionary for this component
                old_label_dict[component].add(original_label)

            # otherwise, add to the invalid file names set
            else:

                invalid_file_names.add(original_file)

    # if this component is not in the original file name, and therefore its value is nan, it will be a tomlkit.items.Float class
    elif isinstance(original_position, float):
        continue

    # if it is neither a float or an integer, an incorrect value was provided, so print ERROR and exit
    else:
        print(f'ERROR. The value for the component {component_name} in the file rename configuration file must '
              f'be either an integer, or \'nan\' if this component does not exist in the original file name. The '
              f'provided component value, {original_position}, is of type {type(original_position)}. Please '
              f'update the file rename configuration file and then rerun file renaming.\n')
        # sys.exit(1)


##########################################################
## MATCH ORIGINAL LABELS TO UPDATED LABELS ###############
##########################################################
print('section 06: match original labels to updated labels\n')
# go through each of the updated file name components, and sort the old label values to their corresponding new labels
for component in file_name_conversion.keys():

    # get the old labels that were found in the previous step for this component
    old_labels = old_label_dict[component]

    # if this component is not present in any original file names...
    if len(old_labels) == 0:

        # get the value to use from the missing_components section of the configuration file
        missing_component = config_dict['missing_components'][component].lower()

        # if a value is provided...
        if len(missing_component) > 0:

            # use this value as the 'old' label so that it is clear which new label to add
            file_name_conversion[component][missing_component] = missing_component

        # if no substitute value is provided...
        else:

            # print error message and exit
            print(f'ERROR. Missing a substitute value for the missing component {component} in the original file '
                  f'name. Please update the configuration file so that a default value for this component can be '
                  f'used despite it missing from the original file name.\n')
            # sys.exit(1)

    # if this component was found in any original file names...
    else:

        # go through each original label found, and sort into the corresponding new label
        for label in old_labels:

            # handle sorting of compartments for soil and litter
            if component == 'compartment':
                for new_lab, old_lab in config_dict['accepted_labels'][component].items():  # key in config corresponds to new label
                    # this script currently ignores the wildcards from the config file, so remove here or won't match
                    old_lab = old_lab.replace('*', '')
                    if old_lab == label:
                        file_name_conversion[component][new_lab].add(old_lab)

            # pull the ecoregion values from the neon_domains dictionary
            elif component == 'ecoregion':
                file_name_conversion[component].update(extract_nested(neon_domains))

            # sort treatment (burn/unburn, conifer/oak/grassland)
            elif component == 'treatment':

                # separate burn treatment and habitat from this label
                burn_treatment, habitat = label.split('-')

                # construct new label based on old, and use this to add to correct dictionary key
                new_label = ''

                # burn/unburn treatment in old label is one or two letters; find correct one
                if re.search(r'^U', burn_treatment, re.I):
                    new_label += 'U'
                elif re.search(r'^B', burn_treatment, re.I):
                    new_label += 'B'

                # habitat in old label, like new, is only one letter (matches)
                new_label += habitat

                # use this new label to add entire old label
                file_name_conversion[component][new_label].add(label)

            # remove S prefix from subplot values, add leading zero
            elif component == 'subplot':

                try:

                    # search for a number in the subplot label, excluding 0s as a 0 will be added to the number below
                    subplot_num = re.search(r'[1-9]', label).group(0)

                except AttributeError:
                    print(f'ERROR. The subplot number cannot be located in the original label for subplot for '
                          f'this group of samples. Please confirm that the information in the configuration file '
                          f'for file renaming is consistent with the original file naming convention used. \n')
                    # sys.exit(1)

                # add leading zero (all values are single-digits, no need to check whether leading zero necessary)
                num_with_zero = '0' + subplot_num

                # use this updated label to add original label to conversion dict
                file_name_conversion[component][num_with_zero].add(label)


##########################################################
## CREATE NEW FILE NAMES FOR EACH OLD FILE NAME ##########
##########################################################
print('section 07: create new file names for each old file name\n')
# create a conversion dict for the file paths, important for renaming process
conversion_path_dict = {k:'' for k in original_file_paths}

# create a search regex for finding the read orientation for a file name
#  this will be added to if it doesn't work for the original file names, and a user provides input
#  i.e., it will learn and use this information for future file names as well
read_orientation_regex = r'R1|R2'

# create an empty list to add any original file names that couldn't be renamed using this section of the code
cannot_rename = []

# go through each original file and assemble its new file name and file path
for file in conversion_path_dict.keys():

    # set Boolean to determine whether this file name can be automatically renamed or not
    valid_name = True

    # get just the file name from this file path
    original_file_name = file.stem.replace('.fastq', '')  # stem will only return one of two file extensions, so use replace too

    # create a dictionary for this file to add the updated file name components to
    update_dict = {'seq_platform': '',
                   'compartment': '',
                   'coll_date': '',
                   'ecoregion': '',
                   'treatment': '',
                   'subplot': '',
                   'read_orientation': ''}

    # create a list of the original labels used in the file name for this file
    original_labels = original_file_name.split(file_delim)

    ## FIND THE NEW LABEL, GIVEN THE OLD LABEL FOR EACH FILE NAME COMPONENT

    # go through each component and find the corresponding new label for this original file name
    for component in update_dict.keys():

        # check if the file name component is a key in the component_positions section of the config file
        if component in config_dict['component_positions'].keys():

            # get the position of this component from the original file name
            original_position = config_dict['component_positions'][component]

            # use the position from the config file to get the old label for this file
            if isinstance(original_position, int):
                original_label = original_labels[original_position]
            else:
                continue

            # get a list of old labels for this file name component, to confirm the original label in this list
            old_labels = flatten_dict_values(file_name_conversion[component])

            # add the new label to the update dictionary for this component
            #   (if original label in the list of old labels)
            if original_label in old_labels:
                update_dict[component] = search_conversion_dict(file_name_conversion[component], search_val=original_label)
            else:
                print(f'WARNING. The label {original_label} for {component} in the original file name,\n'
                      f'  {original_file_name}\n'
                      f'was not located in the file name conversion dictionary: {old_labels}. '
                      f'This is likely a coding issue. Check the code to see where this error may have occurred.\n')
                valid_name = False
                no_warnings = False

        # if it is not a key in the component_positions section, it will be handled in the next section, so skip over it
        else:
            continue

    ## CHECK FOR MISSING FILE NAME COMPONENTS FOR THE NEW FILE NAME LABELS

    # once all components are searched for the new label, check for any missing fields
    for comp, label in update_dict.items():

        # if the label is missing after searching the file_name_conversion dictionary...
        if label == '':
            # check the configuration file for any additional information that may be provided for this missing label
            # if this component has a key in the config file under missing_components section...
            if comp in config_dict['missing_components']:
                # check if it has a corresponding value as well
                if config_dict['missing_components'][comp] != '':
                    # if a value is found for this component, use this for the new label
                    update_dict[comp] = config_dict['missing_components'][comp]
                else:
                    # if a value is not found, but there is a place in the config to add a value, then print warning
                    print(f'WARNING. The file name:\n'
                          f'  {original_file_name}\n'
                          f'is missing information for the file name component {comp}. If this information was not '
                          f'included in the original file name, then add this information to the configuration file '
                          f'in the section missing_components, and rerun renaming.\n')
                    valid_name = False
                    no_warnings = False

            # read orientation is expected to be missing, but will be located from the original file name with regex
            elif comp == 'read_orientation':
                # search using standard regex
                try:
                    # search for the read orientation using the regex at the top of this section
                    read_orient = re.search(read_orientation_regex, original_file_name, re.I).group(0)
                # if this doesn't return a match, ask user to add value for match
                except AttributeError:
                    # if a match isn't found using the standard regex, ask user for string of the file name
                    #  that corresponds to the read orientation
                    print(f'WARNING. The read orientation could not be located in the file name:\n'
                          f'  {original_file_name}\n'
                          f'Please enter the read orientation string used in this file name:\n')
                    no_warnings = False
                    # add the provided read orientation regex to the search regex, with or '|' operator
                    read_orientation_regex = read_orientation_regex + '|' + input()
                    read_orient = re.search(read_orientation_regex, original_file_name, re.I).group(0)

                update_dict[comp] = read_orient

            # collection date may be missing (it is from UO samples for soil/litter), so check config section to get
            #   the collection date for this site
            elif comp == 'coll_date':
                # get the name of the original site name used
                original_site_name = original_file_name.split(file_delim)[config_dict['component_positions']['ecoregion']]
                # look up the corresponding state for this site name
                for state, site in config_dict['site_name'].items():
                    # if there is more than one site in this state (i.e., Oregon)
                    if isinstance(site, list):
                        # iterate through the sites associated with this state
                        for i in range(len(site)):
                            # look for the site name that matches the original site name in the file name
                            if site[i] == original_site_name:
                                # use the index value for the match to get the correct collection date for this state
                                update_dict[comp] = config_dict['collection_date'][state][i]
                    # if there is only one site in this state
                    else:
                        # and the site matches the original site name
                        if site == original_site_name:
                            # add the collection date for this site/state to the update dictionary for the new file name
                            update_dict[comp] = config_dict['collection_date'][state]

            # if a file name component is missing that is expected to be present, tell user to update the config file
            #  with a key/value pair for this file name component in the configuration file
            else:
                print(f'WARNING. The file name:\n'
                      f'  {original_file_name}\n'
                      f'is missing information for the file name component {comp}. This information does not currently '
                      f'have a section in the configuration file under the missing_components section. If this component '
                      f'is missing from the original file name, add a key/value pair for this file name component in '
                      f'the missing_components section of the configuration file, using one of the following as the '
                      f'key:')
                print_indented_list(list(update_dict.keys()))
                valid_name = False
                no_warnings = False

    ## USE THE NEW LABELS TO CREATE A NEW FILE NAME AND FILE PATH

    # first check that no warning has been triggered above
    if valid_name:
        # join the labels from the update dictionary; already order in the correct order
        updated_file_name = '_'.join(list(update_dict.values()))

        # create a file path from the new file name; use the same parent directory as the original files
        #  so that the original file names are replaced by the new file names
        updated_file_path = (file.parent / (updated_file_name + '.fastq.gz'))

        conversion_path_dict[file] = updated_file_path

    # if warning was triggered, do not create a new file name and instead add the original file name to error list
    else:
        cannot_rename.append(file)

##########################################################
## RENAME ORIGINAL FILES WITH NEW CONVENTION #############
##########################################################
print('section 08: rename original files with new convention\n')

## CHECK FOR DUPLICATE NEW NAMES
print('section 08.01: check for duplicate new names\n')
# create a list of all new file names, and then create a set from this list (to keep only the unique values)
all_new_names = [new_path for new_path in conversion_path_dict.values() if not new_path == '']
unique_new_names = set(all_new_names)

# if the number of unique values matches the number of all new file names, do nothing here (means no duplicates)
if len(all_new_names) == len(unique_new_names):
    pass
# if there is a difference between the unique values and all file names, then determine which file names are duplicated
else:
    duplicate_new_names = []  # list to add duplicate file names to
    unique_new_names = set()  # recreate this set, but only add to it if it isn't already in this set

    # go through each new file name
    for new_name in all_new_names:

        # if this new name is already in the set of unique new names (meaning it is a duplicate), add to duplicate list
        if new_name in unique_new_names:
            duplicate_new_names.append(new_name)

        # if this new name isn't yet in the set of unique new names, then add it
        else:
            unique_new_names.add(new_name)

    # check if there are duplicate names (there should be, if here)
    if len(duplicate_new_names) > 0:
        print(f'The new file names contain duplicate names for the following files:\n')
        # go through each duplicate file name and print the new name and the original file name
        for old_path, new_path in conversion_path_dict.items():
            if new_path.stem in duplicate_new_names:
                print(f'{new_path.stem}  /   {old_path.stem}')

## TRY TO FIX ANY INVALID NEW NAMES
print(' section 08.02: try to fix invalid new names\n')
for invalid_path in cannot_rename:

    invalid_name = invalid_path.stem.replace('.fastq', '')

    # negative control sequence files
    if re.search(NEG_CTRL_RE, invalid_name, re.I):
        try:
            read_orientation = re.search(read_orientation_regex, invalid_name, re.I).group(0)
            neg_ctrl_num = re.search(r'(?<=ntc)\d', invalid_name, re.I).group(0)
            replacement_file = '_'.join(['illumina', 'soil-litter', run_date, f'NTC-{neg_ctrl_num}', read_orientation])
            replacement_name = invalid_path.parent / (replacement_file + '.fastq' + '.gz')
        except AttributeError:
            replacement_name = 'NA'

    # mock community sequence files
    elif re.search(MOCK_COMM_RE, invalid_name, re.I):
        try:
            read_orientation = re.search(read_orientation_regex, invalid_name, re.I).group(0)
            replacement_file = '_'.join(['illumina', 'soil-litter', run_date, 'mock-community', read_orientation])
            replacement_name = invalid_path.parent / (replacement_file + '.fastq' + '.gz')
        except AttributeError:
            replacement_name = 'NA'

    # undetermined reads, those that couldn't be sorted into samples during demultiplexing
    elif re.search(UNDET_RE, invalid_name, re.I):
        try:
            read_orientation = re.search(read_orientation_regex, invalid_name, re.I).group(0)
            replacement_file = '_'.join(['illumina', 'soil-litter', run_date, 'undetermined', read_orientation])
            replacement_name = invalid_path.parent / (replacement_file + '.fastq' + '.gz')
        except AttributeError:
            replacement_name = 'NA'

    # usually just non-climush samples, but also a catch-all for any other unsortable files
    else:
        new_dir_name = invalid_name.split('_')[0].lower()
        new_dir_parent = mkdir_exist_ok('non-climush_sequences', parent_dir=input_path.parent)  # make dir for non-climush seqs
        new_dir_path = mkdir_exist_ok(new_dir_name, parent_dir=new_dir_parent)  # make subdir for these seqs by group
        replacement_name = new_dir_path / (invalid_name + '.fastq' + '.gz')

    conversion_path_dict[invalid_path] = replacement_name

## RENAME FILES
print('  section 08.03: rename files\n')
renamed_file_counter = 0
old_names = []
new_names = []
not_renamed = []
for old_path, new_path in conversion_path_dict.items():

    # add the old name to the list of old names
    try:
        old_names.append(old_path.stem.replace('.fastq', ''))
    except AttributeError:
        print(f'old path issue: {old_path}\n'
              f'  object type = {type(old_path)}')
        continue

    if is_pathclass(new_path, exit_if_false=False):

        new_stem = new_path.stem.replace('.fastq', '')

        # print(f'Renaming {old_names[-1]} to {new_stem}...\n')

        # add to the counter
        renamed_file_counter += 1

        # rename the file
        old_path.replace(new_path)

        # add new name to the list of new names
        new_names.append(new_stem)

    else:

        # print(f'NOT renaming {old_names[-1]}...\n')

        # if the new_path doesn't exist (an error file), just add an NA value
        new_names.append('NA')

        # add the full path to a list of file paths that were not renamed (hard to get full path from just name later)
        not_renamed.append(old_path)


print(f'{renamed_file_counter} files out of the {len(original_file_names)} input files were renamed or sorted.\n')

# write a summary file of original file names that were not renamed (if any)
if len(not_renamed) > 0:

    # sort these files that couldn't be renamed into their own directory, since looping through anyways
    not_renamed_dir = mkdir_exist_ok('not_renamed', parent_dir=input_path)

    # create a path to write the summary file to, within the new not_renamed/ directory
    not_renamed_summary_path = (not_renamed_dir / 'not_renamed').with_suffix('.txt')

    print(f'WARNING. {renamed_file_counter} files out of the {len(original_file_names)} input files '
          f'were renamed. To see the {len(not_renamed)} files that could not be renamed, see the summary '
          f'file:\n'
          f'  {not_renamed_summary_path}\n'
          f'and the files that were not renamed in the directory:\n'
          f'  {not_renamed_dir}\n')
    no_warnings = False

    # write the unnamed file names to a summary file
    with open(not_renamed_summary_path, 'wt') as fout:

        # create a unified header for the document
        fout.write(f'the following files were not renamed by {__file__}:\n')

        # write the name of the shortened file path, with just the parent directory and the original file name
        for not_renamed_path in not_renamed:

            # write the path of the original file, shortened to just the parent and the original file name
            original_file_name = not_renamed_path.stem.replace('.fastq', '')
            formatted_name = not_renamed_path.parent.name + '/' + original_file_name
            fout.write(f'  {formatted_name}\n')

            # move the file (without renaming the actual file name) to the not_renamed directory
            dest = not_renamed_dir / not_renamed_path.name  # use the name so it will have the .fastq.gz file ext still
            not_renamed_path.rename(dest)

# export summary conversion table, unless --no-track flag invoked
if no_track:
    pass
else:
    conversion_table = pd.DataFrame({'old_name': old_names,
                                     'new_name': new_names})
    conversion_table.to_csv(rename_path, index=False)
    if rename_path.is_file():
        print(f'The conversion table was successfully written to:\n'
              f'  {rename_path}\n')
    else:
        print(f'There were issues saving the conversion table successfully to:\n'
              f'  {rename_path}\n')
        sys.exit(73)  # can't create output file

if no_warnings:
    print(f'SUCCESS! {__file___} has completed.\n')
    sys.exit(0)
else:
    print(f'Please see warnings generated by {Path(__file__).name}.\n')
    sys.exit(0)  # no standard exit code for warnings; only 0 if not strict (same as success) or 1 if strict (same as failure)

##########################################################
## ADDRESS ANY ISSUES OR ERRORS ##########################
##########################################################
#
# # check if
# nonclimush_path = input_path.parent / 'non-climush_sequences'
# if positive_response(fix_error_files):
#     able_to_rename = [f for f in invalid_filenames if re.search(f'{MOCK_REGEX}' + r'|' + f'{NTC_REGEX}' + r'|' + f'{UNDET_REGEX}', f, re.I)]
#     cannot_rename = list(set(invalid_filenames).difference(able_to_rename))
#     rename_sequence_files(able_to_rename)
#     if len(cannot_rename) == 0:
#         print('All invalid file names were successfully renamed.\n')
#     else:
#         move_unnammed(dest_path, cannot_rename)
# else:
#     move_unnammed(dest_path, invalid_filenames)
#
# # check if any files were not renamed because they did not match the
# #   original naming convention as provided in the configuration file
# if len(invalid_filenames) > 0:
#
#     # get examples of filenames that are invalid because they are non-template controls (negative ctrl) or mock communities
#     invalid_ntc = [file.stem for file in invalid_filenames if re.search(NEG_CTRL_RE, file, re.I)]
#     invalid_mock = [file.stem for file in invalid_filenames if re.search(MOCK_COMM_RE, file, re.I)]
#
#     # gather list of invalid examples from ntc and mock, if available
#
#     # if there are no detected ntc or mock files, then presumably nothing can be done here
#     if len(invalid_ntc + invalid_mock) == 0:
#
#         # make a directory to add these invalid file name files to, in the sequencing run's parent directory
#         still_invalid_dir = mkdir_exist_ok(new_dir='renaming_errors', parent=input_path.parent)
#
#         # move any files that could not be renamed to this renaming error directory
#         for invalid_file in invalid_filenames:
#             invalid_file.rename(still_invalid_dir / invalid_file.name)
#
#         # get run time of renaming process
#         end_time = datetime.now()
#         run_time = end_time - start_time
#
#         # print warning message and exit
#         print(f'WARNING. {len(invalid_filenames)} of {len(original_filenames)} original files in:\n'
#               f'\t{input_path}\n'
#               f'were not renamed because they did not follow the original naming convention as provided in the file '
#               f'renaming configuration file. All files that could not be renamed at this point have been moved to a '
#               f'new directory in this sequencing run\'s main directory:\n'
#               f'\t{still_invalid_dir}\n'
#               f'Exiting {Path(__file__).name}...\n'
#               f'  total run time: {str(run_time)}\n')
#         sys.exit(1)
#
#     # if there are any ntc or mock files, prompt user if they would like them renamed, as suggested
#     else:
#         invalid_examples = []
#         if len(invalid_mock) > 0:
#             invalid_examples.append(invalid_mock[0])
#         if len(invalid_ntc) > 0:
#             invalid_examples.append(invalid_ntc[0])
#
#         # compile warning message
#         fix_error_files = (f'WARNING. {len(invalid_filenames)} of {len(original_filenames)} original files in:\n'
#                            f'\t{input_path}\n'
#                            f'were not renamed because they did not follow the original naming convention as provided '
#                            f'in the file renaming configuration file.\n'
#                            f'If these sequencing files belong to the CliMush project (e.g., {invalid_examples}), '
#                            f'they can still be renamed by using the first label of the old file name for these files '
#                            f'in the new label. Do you wish to continue with renaming these files in this way?')
#
#         # prompt user for yes/no/quit input; will continue if 'yes', will quit here if 'no'/'quit'
#         prompt_yes_no_quit(message=fix_error_files)
#
#         ## I DID NOT HAVING ANY FOLLOWING STEPS HERE IN THE ORIGINAL VERSION, NOT SURE OF MY PLAN HERE... #############
#
#         # use same directory as for the files that can't be renamed, aren't climush, can't be otherwise sorted, etc.
#         still_invalid_dir = mkdir_exist_ok(new_dir='renaming_errors', parent=input_path.parent)
#
#         # move any files that could not be renamed to this renaming error directory
#         for invalid_file in invalid_filenames:
#             invalid_file.rename(still_invalid_dir / invalid_file.name)
#
#         # still calculate run time of the renaming process
#         end_time = datetime.now()
#         run_time = end_time - start_time
#
#         # print error statement, explain that this part of this script is incomplete
#         print(f'ERROR. This part of the pipeline in {Path(__file__).name} is incomplete and not yet utilizable. Sorry '
#               f'for any inconveniences. All files that could not be renamed at this point have been moved to a '
#               f'new directory in this sequencing run\'s main directory:\n'
#               f'\t{still_invalid_dir}\n'
#               f'Exiting {Path(__file__).name}...\n'
#               f'  total run time: {str(run_time)}\n')
#         sys.exit(1)
#
# # otherwise, if the script has not yet exited, then everything has run successfully and as expected
# end_time = datetime.now()
# run_time = end_time - start_time
# print(f'SUCCESS. {rename_success_count} sequencing files were successfully updated to the current file naming '
#       f'convention and can be found in the following path:\n'
#       f'\n')
# print(f'Exiting..n\n'
#       f'  total run time: {str(run_time)}\n')
# sys.exit(0)