import argparse, re, pathlib, warnings
from logging import warning
from pathlib import Path
from datetime import datetime
import pandas as pd
from climush.utilities import exit_process, flag_multiple_files

## SET SCRIPT-WIDE CONSTANTS ###########################################################################################

## SCRIPT TEST MODE ##

# if set to True, will bypass the argparse parser and use pre-defined test values instead
TEST_MODE = True


## CONVERSION TABLE GLOB ##

# glob to use to auto-detect the original file name conversion table, if the argparse argument --table isn't used
CONVERT_TABLE_SUFFIX = 'file-name-conversion.'


## FILENAME ERROR REGEX ##

# locates the part of the spore file name that contains the file naming error
FILENAME_ERR_RE = r'(?<=202(2|3)-)\d__'


## TODAYS DATE ##

# get the current date, used to tag the updated file name conversion table output
TODAY = datetime.now().strftime('%Y-%m-%d')

########################################################################################################################


## DEFINE COMMAND LINE ARGUMENTS #######################################################################################

if TEST_MODE:

    # alert the user if the script is being run in testing mode, in case it is accidental
    warning_msg = (f'This script is currently running in test mode and uses pre-defined settings rather than'
                   f'using any provided command line options.')
    warnings.warn(warning_msg, category=Warning)

    # create a default set of input values in place of actual command line input
    args = {
        'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/spores/test-files/illumina_spore_fnames/renamed_seq_test/'),
        'corrections': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/spores/test-files/illumina_spore_error-filenames.txt'),
        'table': None,  # auto-detect the conversion table; default value if nothing provided to this parameter
    }

else:

    ## INSTANTIATE PARSER ##

    parser = argparse.ArgumentParser(
        prog=Path(__file__).stem,
        description='Fix Illumina spore sequence file names that contain extraneous underscores.',
        epilog='This script is part of the CliMush bioinformatics pipeline.',
    )

    ## INPUT ##

    # [REQUIRED] parent file path
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=pathlib.PosixPath,
        help='Path to the parent directory containing the spore sequence files that need to be renamed.'
    )

    # [REQUIRED] file of filenames to correct
    parser.add_argument(
        '--corrections',
        required=True,
        type=pathlib.PosixPath,
        help='Path to the .txt file containing a list of the sequence file names that are incorrect and need to '
             'be updated using this script.'
    )

    # file name conversion table
    parser.add_argument(
        '-t', '--table',
        required=False,
        type=pathlib.PosixPath,
        help='Path to the file name conversion table that will be updated with the corrected updated file name, so '
             'that it can be linked to the original file name after correcting the udpated file name. If one is not '
             'provided, the script will attempt to locate this table within the parent directory of the input sequence '
             'files provided to the -i / --input argument. '
    )

    ## PARSE COMMAND LINE INPUT INTO DICTIONARY ##

    args = vars(parser.parse_args())

########################################################################################################################


## LOCATE + IMPORT FILE NAME CONVERSION TABLE ##########################################################################

# if the path to the conversion table is not provided, locate the file automatically
if args['table'] is None:

    # conversion table should be located in same directory as the sequence files
    convert_table_list = [ file for file in args['input'].glob(f'*{CONVERT_TABLE_SUFFIX}*') ]

    # if only a single match is located (which is ideal), this is the conversion table to use
    if len(convert_table_list) == 1:
        args.update({'table': convert_table_list[0]})

    # if no matches are found, then raise error
    elif len(convert_table_list) == 0:
        err_msg = (f'A conversion table of the original and updated file names matching the substring '
                   f'\'{CONVERT_TABLE_SUFFIX}\' could not be located in the following input directory:\n'
                   f'   {args["input"]}\n'
                   f'Please provide a path to the conversion table using the argument --table for this script and '
                   f'rerun.\n')
        exit_process(err_msg)

    # if many matches are found, then ask user to specify the correct conversion table
    else:
        flag_multiple_files(
            file_path=args['input'],
            search_for=CONVERT_TABLE_SUFFIX,
            auto_respond=False,
        )

# if a path is provided to the conversion table, no search / nothing more
else:
    pass

########################################################################################################################



## CREATE ERROR FILENAME DICTIONARY ####################################################################################

## IMPORT CORRECTION FILE LIST ##

# import the error filenames into an error filenames list
with open(args['corrections'], 'rt') as errors_in:
    error_filenames = [ line.strip() for line in errors_in.readlines() ]

## CREATE ORIGINAL CONVERSION DICTIONARY ##

# import the original conversion table that has the original file names and the error file names
original_convert_table = pd.read_csv(args['table'])

# convert the conversion table to a conversion dictionary
original_convert_dict = {
    original_fname: updated_fname for \
    original_fname, updated_fname in \
    zip(original_convert_table['filename_original'], original_convert_table['filename_updated'])
}

########################################################################################################################


## FIX UPDATED FILE NAMES IN CONVERSION TABLE ##########################################################################

# create a copy of the conversion dictionary to update with the fixed updated file names
updated_convert_dict = original_convert_dict.copy()

# keep track of how many update error file names were located in the original conversion dictionary
found_error_fnames = 0

# iterate through all filenames in the conversion dictionary; use the original to iterate so the copy can be updated
for original_fname, updated_fname in original_convert_dict.items():

    # if the updated file name from the original conversion table matches one of the error file names from
    #   the .txt file of error file names...
    if updated_fname in error_filenames:

        # locate the part of the file name that contains the renaming error
        err_result = re.search(FILENAME_ERR_RE, updated_fname, re.I)

        # check that the regex successfully located the error portion of the filename
        if err_result:

            # if updated file name was in the list of error files and also matched the regex, add to counter
            found_error_fnames += 1

            # add a leading zero to the number and drop the extra underscore
            fixed_filename_portion = '0' + ''.join(list(err_result.group(0))[:-1])

            # replace the error part of the updated file name with this fixed portion
            fixed_filename = re.sub(FILENAME_ERR_RE, fixed_filename_portion, updated_fname)

            # update the copy of the conversion table dictionary so that the original file name associated
            #   with this error updated file name is now associated with the fixed updated file name
            updated_convert_dict.update({original_fname: fixed_filename})

            print(f'error file name     = {updated_fname}\n'
                  f'corrected file name = {fixed_filename}\n')

        # if the regex didn't work, print out the issue
        else:
            print(f'The updated file name from the filename conversion table, {args["table"].name}:\n'
                  f'    {updated_fname}\n'
                  f'was located in the list of error name files, {args["corrections"].name}, but the '
                  f'location of the error in the name could not be detected using the regular expression '
                  f'{FILENAME_ERR_RE}.')

    # continue to the next pair of old and updated file names if the updated file name is not in the
    #   .txt file of error file names
    else:
        continue

# print summary of how many file names were located from the conversion table and fixed
print(f'{found_error_fnames} updated file names from the filename conversion table, {args["table"].name}, were '
      f'matched to an error file in the list of error name files, {args["corrections"].name}, which contained '
      f'{len(error_filenames)} error files to find and fix.\n')

########################################################################################################################


## RENAME CORRECTED UPDATED FILE NAMES #################################################################################

# go through the updated file names conversion dictionary, and rename the error updated files
for original_fname, corrected_fname in updated_convert_dict.items():

    # only rename if they have been changed
    if original_convert_dict[original_fname] in error_filenames:

        # create a path for the error updated file name from the original conversion table, using the common
        #  original file name between the error and fixed updated file names
        error_fpath = args['input'] / original_convert_dict[original_fname]
        corrected_fpath = args['input'] / corrected_fname

        # replace the error updated file name with the corrected updated file name
        error_fpath.replace(corrected_fpath)

    # if there were no corrections made to this updated file name, skip over
    else:
        pass

########################################################################################################################

## WRITE OUT THE UPDATED CONVERSION TABLE ##############################################################################

# create an output path for this conversion table, adding the date to distinguish it from the original
updated_convert_path = (args['input'] / (args['table'].stem + f'_{TODAY}')).with_suffix('.csv')

pd.DataFrame.from_dict(
    data=updated_convert_dict,
    orient='index',
    columns=['filename_updated'],
).reset_index(names='filename_original').to_csv(updated_convert_path, index=False)

########################################################################################################################
