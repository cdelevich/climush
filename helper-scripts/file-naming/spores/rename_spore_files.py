from pathlib import Path
import pathlib, re, argparse
import pandas as pd

from climush.constants import SEQ_FILE_GLOB, ILLUMINA_SEQ_ORIENT_RE, ANY_PLATFORM_REGEX
from climush.utilities import get_settings, mkdir_exist_ok, copy_original_files

## IMPORT BIOINFORMATICS PIPELINE CONFIGURATION SETTINGS ###############################################################

# set a reference directory for a location from which to start looking for the bioinformatics configuration file
ref_dir = Path.cwd().parent.parent

# use get_settings() to locate and read in the bioinformatics configuration file
settings = get_settings(ref_dir)

########################################################################################################################


## DEFINE COMMAND LINE OPTIONS #########################################################################################

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Rename the demultiplexed spore sequence files from '
                                             'the University of Minnesota.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files to rename
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='Path to the directory containing the demultiplexed spore sequence files to rename.')

# output directory for the renamed files
parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='Path to a directory to write the renamed spore sequence files to.')

args = vars(parser.parse_args())

########################################################################################################################


## READ IN SEQUENCE FILES TO RENAME ####################################################################################

# create a list of all files in the input directory that are sequence files
original_file_names = list(args['input'].glob(SEQ_FILE_GLOB))

# how many sequence files need to be renamed
print(f'{len(original_file_names)} sequence files need to be renamed in the directory:\n'
      f'   {str(args["input"])}')

########################################################################################################################


## STORE THE ORIGINAL FILE NAMES IN A DICTIONARY #######################################################################

# create a dictionary where the key is the original file name and the value will be the updated file name
file_name_conversion = {
    og_fname:'' for og_fname in original_file_names
}

########################################################################################################################


## CREATE A DICTIONARY OF ORIGINAL AND NEW LABELS FOR NON-SAMPLE FILES #################################################

nonsample_relabel = {
    'Synmock': 'mock-community',
    'NegCon': 'NTC',
    'PosCon': 'PTC',
}

########################################################################################################################


## CREATE NEW FILE NAMES USING INFO FROM ORIGINAL FILE NAMES ###########################################################

# go through each of the files stored in the file conversion dictionary
for original_fpath in file_name_conversion:

    # get the file extension of the original file name to use for the new file name
    file_ext = ''.join(original_fpath.suffixes)

    # pull only the file name from the original file path and
    #   replace all separators with the correct separator, an underscore
    original_fname = original_fpath.name.replace('-', settings['formatting']['filename_delim'])

    # locate the read orient in the original file name
    read_orient = re.search(ILLUMINA_SEQ_ORIENT_RE, original_fname, re.I).group(0)

    # search for the sequence platform (e.g, illumina) in file name
    seq_platform_result = re.search(ANY_PLATFORM_REGEX, original_fname, re.I)

    # sample sequence files (NOT mock community, pos/neg control) should have the sequence platform
    #   in the original file name
    if seq_platform_result:

        # return the match to the sequence platform from the search of the file name
        seq_platform = seq_platform_result.group(0)

        # all compartment labels should be spore
        compartment_result = re.search(r'(?<=_)spore(?=_)', original_fname, re.I)
        compartment_position = compartment_result.span()[1]  # where does compartment end?
        compartment = compartment_result.group(0)

        # get the NEON ecoregion domain
        neon_domain = re.search(r'(?<=_)D\d{1,2}(?=_)', original_fname, re.I).group(0)

        # get the treatment
        treatment_result = re.search(r'(?<=_)[A-Z]{2}(?=_)', original_fname, re.I)
        treatment_position = treatment_result.span()[1]  # where does treatment end?
        treatment = treatment_result.group(0)

        # get the subplot number
        subplot_start = treatment_position + 1
        subplot_end = subplot_start + 2
        subplot = original_fname[subplot_start:subplot_end]

        # format the date properly

        # year is directly after the compartment
        year_start = compartment_position + 1  # add one to account for file delimiter between
        year_end = year_start + 2  # year is two numbers in original filename
        year_original = original_fname[year_start:year_end]

        # reformat the year by adding '20' to start of string
        year = f'20{year_original}'

        # month is directly after the year
        month_start = year_end + 1  # add one to account for file delimiter between
        month_end = month_start + 2  # month is two numbers in original filename
        month = original_fname[month_start:month_end]

        # combine year and month with hyphen to complete date label
        collection_date = f'{year}-{month}'

        # create new file name with the new labels
        updated_fname = settings['formatting']['filename_delim'].join([
            seq_platform,
            compartment,
            collection_date,
            neon_domain,
            treatment,
            subplot,
            read_orient,
        ])

    # non-sample sequence file (e.g., mock community, pos/neg control) won't have sequence platform in name
    else:

        # the sequence platform for these will all be illumina, so just hard-code that here
        seq_platform = 'illumina'

        # all will be spores as well, so also hard-code here
        compartment = 'spore'

        # use the original-new conversion dict for non-sample files to get the file tag
        nonsample_tag = [updated_tag for original_tag, updated_tag in nonsample_relabel.items() if
                         re.search(original_tag, original_fname, re.I)][0]

        # if the new non-sample tag is PTC, these were numbered in the original file names, so find the associated number
        if nonsample_tag == 'PTC':

            # find the number associated with the original file name tag
            nonsample_tag_num = re.search(r'(?<=PosCon)\d{1}', original_fname, re.I).group(0)

            # add this same number to the updated file name tag, separated with a hyphen
            nonsample_tag = f'{nonsample_tag}-{nonsample_tag_num}'

        # otherwise, the other non-sample tags aren't numbered, so do nothing
        else:
            pass

        # unclear which collection date to add to which non-sample files, so awaiting response from Peter Kennedy
        #  to decide which dates (if any) to use for these non-sample files

        # create new file name with the new labels
        updated_fname = settings['formatting']['filename_delim'].join([
            seq_platform,
            compartment,
            nonsample_tag,
            read_orient,
        ])


    # create a new file path using the output file path and the original file name suffix
    updated_fpath = (args['output'] / updated_fname).with_suffix(file_ext)

    # add this file path to the conversion dictionary as the value to its original file path
    file_name_conversion.update({original_fpath: updated_fpath})

########################################################################################################################


## PRESERVE ORIGINAL FILE NAME INFORMATION #############################################################################


## CREATE A COMPRESSED COPY OF THE ORIGINAL FILES ##

# do this first so that there aren't any new directories or files in the input directory

# copy the original files and write out to a new directory, that will then be compressed
copy_original_files(
    directory=args['input'],
    copy_directory='original_files',
    compress=True,
)


## EXPORT ORIGINAL / UPDATED FILE NAME TABLE ##

# create a dictionary of file names from the dictionary of the file name conversions, which are file paths
file_name_conversion_fnames = {
    original_path.name: updated_path.name
    for original_path, updated_path in
    file_name_conversion.items()
}

# create a dataframe from the conversion dictionary; values will still be file paths, not file names
file_name_conversion_df = pd.DataFrame.from_dict(
    data=file_name_conversion_fnames,
    orient='index',
    columns=['filename_updated'],
).reset_index(names='filename_original')

# create the output file path if it doesn't already exist
mkdir_exist_ok(args['output'])

# create an output path to write this dataframe to the output directory
conversion_df_path = args['output'] / f'illumina_spore_file-name-conversion.csv'

# write this dataframe out to the output path constructed above
file_name_conversion_df.to_csv(conversion_df_path, index=False)


########################################################################################################################


## RENAME ORIGINAL FILES WITH NEW FILE NAMES ###########################################################################

# count how many of the files were successfully renamed
renamed_file_count = 0

# go through the file name conversion dictionary and replace the original file names with the updated file names
for original_fpath, updated_fpath in file_name_conversion.items():

    # rename the original file paths with the updated file paths
    original_fpath.replace(updated_fpath)

    # first confirm that the original file path no longer exists
    if original_fpath.is_file():
        raise OSError(f'The original file path still exists after renaming:\n   {original_fpath}\n')
    else:
        pass

    # then confirm that the update file path does exist
    if updated_fpath.is_file():
        # if it exists, add to the renamed file counter
        renamed_file_count += 1
    else:
        raise OSError(f'The updated file path does not exist after renaming:\n   {updated_fpath}\n')

########################################################################################################################


## SUMMARIZE THE DETAILS OF FILE RENAMING ##############################################################################

# calculate the percent of original files that were successfully renamed
renamed_file_percent = (renamed_file_count / len(original_file_names)) * 100

# print header for renaming summary information
print(f'\nSUMMARY OF FILE RENAMING\n')

# number of successfully renamed sequence files and percent represented
print(f'   {renamed_file_count} out of {len(original_file_names)} ({renamed_file_percent:.2f}%) input sequence '
      f'files were successfully renamed in the input directory:\n'
      f'      {args["input"]}')

# location of archived sequence files with original file names
archive_dir_string = str(args['input'] / 'original_files')
print(f'   A copy of the original file names was archived in the following compressed directory:\n'
      f'      {archive_dir_string}')

# location of the file name conversion table
print(f'   A conversion table of the original file names and their corresponding updated files names '
      f'was exported to the following location:\n'
      f'      {conversion_df_path}')


########################################################################################################################