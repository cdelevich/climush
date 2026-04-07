import argparse, pathlib, re
from pathlib import Path
import pandas as pd

# set a variable for the output file type; plan to always use .csv so I didn't add as an option via argparse
output_filetype = '.csv'

## COMMAND LINE ARGUMENTS ##############################################################################################


## PARSER ##

parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='Filter an OTU table to include only specific taxa.',
    epilog='This script is part of the CliMush bioinformatics pipeline.',
)


## FILE INPUT / OUTPUT ##

# input file [REQUIRED]
parser.add_argument(
    '-i', '--input',
    type=pathlib.PosixPath,
    help='The file path to the original OTU table containing the full record of OTUs.',
)

# establish a mutually exclusive group for the output file name and publication file tag options
outfile_opts = parser.add_mutually_exclusive_group()

# output file [optional]
outfile_opts.add_argument(
    '--fout',
    type=pathlib.PosixPath,
    default=None,
    required=False,
    help='The file path to write the filtered OTU table to; if none provided, the output file '
         'will be written to the same location as the input file and include the filtered taxon '
         'group within the file name.',
)

# publication file tag [optional]
outfile_opts.add_argument(
    '--pub',
    type=str,
    default=None,
    required=False,
    help='Output file prefix that denotes the publication from which the input OTU table originated. '
         'This publication prefix will be added to the beginning of the output filtered OTU table. '
         'This argument is only available if a specific output file path was not provided.',
)


## TAXONOMIC FILTERS ##

# create a list of all taxonomic levels that are options to select for as a filter
taxonomic_levels = [
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]

# create a dictionary using each of these taxonomic levels, where the key is the full string and value is
#   just the first letter of the taxonomic level string, to be later formatted as the argument flag
taxonomic_flags = {tax_str: tax_str[0] for tax_str in taxonomic_levels}

# create a mutually exclusive group so that only one taxonomic level can be used
taxfilter_opts = parser.add_mutually_exclusive_group(
    required=True,  # one filtering group is required (i.e., filtering at one taxonomic level out of the 7)
)

# for each taxonomic level, add an argument
for tax_lvl, tax_flag in taxonomic_flags.items():

    # add an argument for this taxonomic level
    taxfilter_opts.add_argument(
        f'-{tax_flag}', f'--{tax_lvl}',
        type=str,
        default=None,
        help=f'A taxonomic group at the {tax_lvl} level by which to filter the input OTU table.',
    )


## SAVE ARG INPUT TO DICTIONARY ##

args = vars(parser.parse_args())

########################################################################################################################


## GET FILTERING TAXONOMIC LEVEL + GROUP ###############################################################################

# pull out the taxonomic filtering information from the argparse command line input as individual strings
# do this first so that the filtering taxonomic string can be used to assemble the output file path, if needed

# create an empty string for each of the variables; this prevents PyCharm from being mad about an indef variable name
taxfilt_lvl = ''
taxfilt_tax = ''

# go through each of the taxonomic level options from the command line input
for tax_lvl in taxonomic_levels:

    # if this taxonomic level was not used as the filter, continue to the next
    if args[tax_lvl] is None:
        continue

    # if there was a value provided via the command line for this taxonomic level, save the taxon string and level
    #   to a string variable (separately)
    else:
        taxfilt_lvl = tax_lvl
        taxfilt_tax = args[tax_lvl]

########################################################################################################################


## FORMAT OUTPUT FILE PATH #############################################################################################

# if no output path is provided, create a default output file path based on the input file path
if args['fout'] is None:

    # use the parent of the input file as the parent of the output file
    output_parent = args['input'].parent

    # create an output file name using...

    # the input file as the base of the output file
    output_basename = args['input'].stem

    # add delimiter to publication prefix, if present
    if args['pub'] is None:
        pub_prefix = ''
    else:
        # check that a delimiter was not already added, then add
        if args['pub'][-1] == '_':
            pub_prefix = args['pub']
        else:
            pub_prefix = args['pub'] + '_'

    # combine the file path and assemble the output file name from its components; replace None with this path in
    #  the argparse parsed argument dictionary
    args['fout'] = (output_parent / f'{pub_prefix}{output_basename}_{taxfilt_tax}').with_suffix(output_filetype)

# if an output path was provided via the command line, use this value as the output path without altering
else:
    pass

########################################################################################################################


## READ IN INPUT OTU TABLE #############################################################################################

## READ IN FILE ##

# read in the OTU table using the file path provided to argparse command line arguments

# encoding may vary among input files, so first try default (utf-8), then try ascii if that fails
try:
    unfiltered_otu_df = pd.read_csv(args['input'], encoding='utf-8')
except UnicodeDecodeError:
    unfiltered_otu_df = pd.read_csv(args['input'], encoding='ascii')


## LEARN COLUMN NAMES ##

# pull the column names from the OTU table
unfiltered_column_names = unfiltered_otu_df.columns

# create a dictionary that matches the taxonomic levels used by argparse to the column names in the unfiltered table
unfiltered_column_match = {
    tax_lvl:'' for tax_lvl in taxonomic_flags
}

# check that there is a column name for each of the possible taxonomic levels by simply checking for the first letter
for unfilt_col in unfiltered_column_names:

    # get the first letter of this column name
    unfilt_col_startlett = unfilt_col[0]

    # check if this first letter matches a taxonomic level first letter
    for tax_lvl_str, tax_lvl_startlett in taxonomic_flags.items():

        # if it matches, link the taxonomic group from the taxonomic flags to the column name in the unfiltered table
        if unfilt_col_startlett == tax_lvl_startlett:
            unfiltered_column_match[tax_lvl_str] = unfilt_col

        # if it doesn't match, try the next taxonomic level
        else:
            continue

# confirm that each taxonomic level string in the match column dictionary has an associated column name

# create a list of any keys in the column match dictionary that have an empty string as a value
nomatch_lvl = [tax_lvl for tax_lvl in unfiltered_column_match if (unfiltered_column_match[tax_lvl] == '')]

# check if this list of empty string value keys is greater than 0
if len(nomatch_lvl) > 0:
    nomatch_lvl_fmt = ', '.join(nomatch_lvl)
    raise KeyboardInterrupt(
        f'The following {len(nomatch_lvl)} taxonomic levels were not detecting in the input OTU table: \n'
        f'{nomatch_lvl_fmt}.\n'
    )

# if list is empty, then all taxonomic levels were matched to a column name in the original dataframe
else:
    pass

########################################################################################################################


## FILTER INPUT OTU TABLE ##############################################################################################

## FILTER TABLE ##

# filter the original OTU table at the respective taxonomic level (column) by the input taxonomic group
filtered_otu_df = unfiltered_otu_df[unfiltered_otu_df[unfiltered_column_match[taxfilt_lvl]] == taxfilt_tax]

## ASSESS FILTERED OUTPUT ##

# print the number of matches located for this taxonomic group; warn if none were found
tax_match_count = filtered_otu_df.shape[0]

# generic output message
output_msg = (
    f'the {taxfilt_lvl} {taxfilt_tax} were located in the column {unfiltered_column_match[taxfilt_lvl]} '
    f'of the input OTU table:\n'
    f'   {args["input"]}\n'
)

# construct the output message based on the number of matches
if tax_match_count == 0:
    warning_msg = 'No matches were located for ' + output_msg
    raise Warning(warning_msg)
else:
    print(f'{tax_match_count} OTUs matching {output_msg}')

########################################################################################################################


## FETCH FILTERED OTU SEQUENCES FROM GENBANK ###########################################################################



########################################################################################################################


## WRITE OUT FILTERED OTU TABLE ########################################################################################

# use the output argparse argument to write out the filtered OTU table to this output path; do not include row numbers
#  as an additional column
filtered_otu_df.to_csv(args['output'], index=False)

########################################################################################################################