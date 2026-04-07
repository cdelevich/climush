import argparse
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

from mycopull.miscell import import_dataframe


date_suffix = datetime.now().strftime('%Y-%m-%d')
climush_data_out = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/pacbio-seqs/get-dsw-samples/data/output')
dsw_withseqs_out = climush_data_out / f'climush_sporocarp-collection-list_DSW_with-sequences_{date_suffix}.csv'
dsw_withoutseqs_out = climush_data_out / f'climush_sporocarp-collection-list_DSW_without-sequences_{date_suffix}.txt'

## REMOVE AFTER TESTING ##########################
# view more columns from pandas dataframes
pd.set_option('display.max_columns', 15)

args = {
    'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/pacbio-seqs/get-dsw-samples/data/input/PacBio All Runs iNat info_dwnld-2025-09-26.csv'),
    # 'output': None,
    'output': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/pacbio-seqs/get-dsw-samples/data/collection-list_DSW.csv'),
    'sites': ['DSW'],
}

args_fmt = '\n   '.join([f'{k} = {v}' for k, v in args.items()])
print(f'WARNING! Running get_site_collection_list in test mode with the following immutable\n'
      f'command line options:\n'
      f'   {args_fmt}\n')

##################################################

## COMMAND LINE OPTIONS ################################################################################################

## INSTANTIATE PARSER ###########################################################

# create an instance of an argparse Argument parser object
parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='Create a list of sample IDs and sequence file names for '
                'a climush site to use for searching and returning sequence '
                'collections from PacBio climush sequence directories.',
    epilog='This script is part of the climush bioinformatics pipeline.',
)

## ADD ARGUMENTS ################################################################

## I/O ##

# input file from which the sample IDs will be pulled
parser.add_argument(
    '-i', '--input',
    required=True,
    type=Path,
    help='Path to the dataframe from which the sample IDs for a given climush site will be located.',
)

# input file from which the sample IDs will be pulled
parser.add_argument(
    '-o', '--output',
    required=False,
    default=None,
    type=Path,
    help='Path to a file or directory to write the output species list to.',
)

## SITE OPTIONS ##

parser.add_argument(
    '--sites',
    required=False,
    type=str,
    nargs='*',
    help='A site or list of sites in the climush project from which to gather PacBio sporocarp sample IDs.',
)


## PARSE ARGUMENTS ##############################################################

args = vars(parser.parse_args())


## VALIDATE I/O PATHS ###########################################################

## MAKE DEFAULT FILE NAME ##

# create a default file name to use for the output
search_sites_fmt = '-'.join(args['sites'])
date_suffix = datetime.now().strftime('%Y-%m-%d')
default_fileout = f'climush_sporocarp-collection-list_{search_sites_fmt}_{date_suffix}.csv'

## CHECK OUTPUT FILE PATH ##

# if no argument was provided to the output file path option...
if args['output'] is None:

    # use the input file path's parent directory as the directory
    output_dir = args['input'].parent
    output_file = default_fileout

# if some argument provided to the output file path option...
else:

    ## OUTPUT ARG IS DIRECTORY ##

    # if the output file path doesn't have a suffix (file extension), assume it is a directory
    if args['output'].suffix == '':

        # if the directory exists...
        if args['output'].is_dir():

            # do nothing; no need to create an output directory
            pass

        # if the directory doesn't exist...
        else:

            # make this directory
            args['output'].mkdir(exist_ok=True)

        # use this as the output directory with the default file name
        output_dir = args['output']
        output_file = default_fileout

    ## OUTPUT ARG IS FILE ##

    # if the output file path is a path to a file
    else:

        # do nothing here, use only what was given as an argument
        pass



########################################################################################################################


## IMPORT DATA #########################################################################################################

# import the data from the dataframe, throwing an error if the data cannot be formatted as a dataframe
input_df = import_dataframe(
    filepath=args['input'],
)

########################################################################################################################


## SUBSET BY SITE ######################################################################################################

# constant for the column name in the input dataframe that has values that will match the site input
SITE_COL = 'Ecoregion'

# empty list to store row indices of rows of the input dataframe matching the site(s)
wanted_row_nums = []

# for each site in the sites to create a collection list for, add their row indices to the list of rows
for site in args['sites']:

    # cant do with list comprehension; will make a nested list of lists
    wanted_row_nums += input_df[input_df[SITE_COL] == site].index.to_list()

# use the site indices (row numbers) to subset the input dataframe to include only these rows
site_only_df = input_df.iloc[wanted_row_nums,:]

# print summary of the number of matches
sites_fmt = '\n   '.join(args['sites'])
print(f'Out of the {input_df.shape[0]:,} rows of the input dataframe, {args["input"].name},\n'
      f'{len(wanted_row_nums):,} rows were from the following site(s):\n'
      f'   {sites_fmt}\n')

########################################################################################################################

## CLEAN DATAFRAME #####################################################################################################

## DROP COLUMNS #################################################################

# keep only the columns to help me create the climush sequence fasta file name for the wanted site samples
wanted_column_labels = [
    'Run date',
    'Ecoregion',
    'Veg Type',
    'Burned?',
    'sample_id',
    'ITS sequence',
]

# remove any unwanted columns
site_only_dropcols = site_only_df[wanted_column_labels]

## RENAME COLUMNS ###############################################################

# rename the columns that I kept to better labels
rename_wanted_columns = {
    'Run date': 'run_date',
    'Ecoregion': 'site_code',
    'Veg Type': 'habitat',
    'Burned?': 'burned',
    'sample_id': 'collection_num',
    'ITS sequence': 'sequence_its',
}

# rename the remaining columns
site_only_dropcols_renamecols = site_only_dropcols.rename(columns=rename_wanted_columns)

## EXTRACT MISSING RUN DATE #####################################################

# locate where samples are missing a sequence run date
missing_rundate_loc = site_only_dropcols_renamecols[site_only_dropcols_renamecols['run_date'].isna()].index.to_list()
if len(missing_rundate_loc) == 0:
    print(f'No rows of the filtered dataframe are missing a run date.\n')
else:
    print(f'{len(missing_rundate_loc)} rows of the filtered dataframe are missing a run date.\n')


## REFORMAT RUN DATE ############################################################

# all run dates should be reformatted to match: YYYY-MM
run_dates_fmt = []
rdate_regex = re.compile(r'\d{4}-\d{2}')
for original_rdate in site_only_dropcols_renamecols['run_date']:

    fmt_rdate_located = rdate_regex.search(original_rdate)

    if fmt_rdate_located:
        run_dates_fmt.append(fmt_rdate_located.group(0))
    else:
        print(f'Run date not located: {original_rdate}\n')
        run_dates_fmt.append(pd.NA)

site_only_dropcols_renamecols['run_date_fmt'] = run_dates_fmt

# drop unformatted run date column
site_only_dropcols_renamecols_fmtrdate = site_only_dropcols_renamecols.drop(columns=['run_date'])


## COPY FINAL CLEANED DATAFRAME #################################################

site_cleaned = site_only_dropcols_renamecols_fmtrdate.copy()

########################################################################################################################


## SUBSET COLLECTIONS W/ SEQUENCES #####################################################################################

# some collections have sequences associated with them and should not have to be re-processed (I guess?)

# with ITS sequence
site_cleaned_seq = site_cleaned[~site_cleaned['sequence_its'].isna()]

# without ITS sequence
site_cleaned_noseq = site_cleaned[site_cleaned['sequence_its'].isna()]


## SUMMARIZE COLLECTIONS W/ AND W/O SEQUENCES ###################################

## GET COUNTS ##

# total number of collections matching site, after cleaning
total_sitematch_count = site_cleaned.shape[0]

# number with ITS sequence
seq_count = site_cleaned_seq.shape[0]
seq_percent = (seq_count/total_sitematch_count)*100

# number without an ITS sequence
noseq_count = site_cleaned_noseq.shape[0]
noseq_percent = (noseq_count/total_sitematch_count)*100


## PRINT SUMMARY ##

print(f'total: {total_sitematch_count:,}\n'
      f'    with sequence:    {seq_count:,} ({seq_percent:.1f}%)\n'
      f'    without sequence: {noseq_count:,} ({noseq_percent:.1f}%)\n')

## EXPORT COLLECTIONS W/ SEQUENCES ##############################################

site_cleaned_seq.to_csv(dsw_withseqs_out, index=False)

########################################################################################################################


## GET SEQUENCE FILE NAME ##############################################################################################

# common filename features
pb_prefix = 'pacbio'
sporo_ctype = 'sporocarp-f'
filename_delim = '_'

sporocarp_fname = []
for r in range(site_cleaned_noseq.shape[0]):

    # get sequence run date
    seqrun_date = site_cleaned_noseq['run_date_fmt'].iloc[r]

    # get collection number
    collnum = site_cleaned_noseq['collection_num'].iloc[r]

    # combine into sequence file basename (no file extension)
    fname = filename_delim.join([
        pb_prefix,
        sporo_ctype,
        seqrun_date,
        collnum,
    ])

    # append this sequence file basename to the list of sporocarp file names
    sporocarp_fname.append(fname)

# add the sporocarp file name to the dataframe
site_cleaned_noseq['sporocarp_fname'] = sporocarp_fname

# export the sporocarp file names into a single text file (.txt) with one filename per line
with open(dsw_withoutseqs_out, 'w') as fnames_out:
    for fname in sporocarp_fname:
        fnames_out.write(f'{fname}\n')

########################################################################################################################