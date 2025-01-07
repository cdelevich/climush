from datetime import datetime
import pandas as pd
import numpy as np
import re, tomlkit, shutil, argparse, pathlib, sys
from pathlib import Path
from climush.constants import SEQ_FILE_GLOB, GZIP_REGEX, SEQ_FILE_RE, CORRECT_CTAB_PREFIX, ANY_CTAB_PREFIX, ANY_CTAB_CODE
from climush.utilities import strip_file_ext

TEST_MODE = True

## COMMAND LINE OPTIONS ################################################################################################

if TEST_MODE:

    print(f'WARNING. Currently running the renaming script in test mode. This means that command line options are not '
          f'accepted and instead default test values for all parameters are used. If this is undesired behavior, '
          f'then set TEST_MODE = False within the rename_2023_endophytes.py script.')

    args = {
        'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/test-files_2023-sequence-files/'),
        'config': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/illumina_endos_2023-10_file-rename-config.toml'),
        'mapping': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/renaming-mapping-files/from-neda-arnold-lab/its1_endophyte-mapping_2023_combined.csv'),
        'output': None,
        'delim': '-',
        'no_copy': False,
    }

else:

    parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                     description='Rename the 2022 endophyte ITS1 sequence files.',
                                     epilog='This script is part of the climush bioinformatics pipeline.')

    parser.add_argument('-i', '--input',
                        required=False,  # if --dry-run flag is used, an input path to seq files to rename is not needed
                        type=pathlib.PosixPath,
                        help='The path to a directory containing the sequencing files that require renaming.')

    parser.add_argument('-c', '--config',
                        required=True,
                        type=pathlib.PosixPath,
                        help='The path to the .toml configuration file for 2022 endophyte renaming.')

    parser.add_argument('-m', '--mapping',
                        required=True,
                        type=pathlib.PosixPath,
                        help='The path to the CTAB code / bag label mapping file for 2022 endophytes.')

    parser.add_argument('-o', '--output',
                        required=False,
                        type=pathlib.PosixPath,
                        help='The path to the directory to write out renaming summary information.')

    parser.add_argument('-d', '--delim',
                        default='-',
                        help='The delimiter used in the bag labels within the mapping file.')

    parser.add_argument('--no-copy',
                        action='store_true',  # if no arg passed, defaults to False
                        help='By default, a copy of the sequencing files with their original file names will be created '
                             'before the files are renamed. If you want to avoid this behavior and not save a copy of '
                             'the original files, then use the --no-copy flag.')

    # parser.add_argument('--dry-run',
    #                     action='store_true',  # if no arg passed, defaults to False
    #                     help='If you want to test this script on empty files matching the sequence files that need to be '
    #                          'renamed, without renaming the actual sequence files, then use this flag.')

    # store arguments in a dictionary
    args = vars(parser.parse_args())

    # # an input path is not required if doing a dry run; if not doing a dry run, then an input path is required
    # if (args['input'] is None) and (args['dry_run'] is None):
    #     print(f'ERROR. Both input path and dry run cannot be None.')
    #     sys.exit()
    # else:
    #     pass

########################################################################################################################


## CUSTOMIZE SETTINGS ##################################################################################################

## DISPLAY ##################################################################

# see more columns
pd.set_option('display.max_columns', 15)

#############################################################################

## RENAMED FILE PREFIX ######################################################

# update: instead of using the sequencing platform as the prefix (e.g., illumina) use the gene region ('its1')
seq_region_prefix = 'its1'

#############################################################################

## REGEX ####################################################################

# ctab file when using rglob; rglob searches absolute file path so cannot use above prefix
CTAB_RGLOB_RE = r'[MSC]{3,5}\d+?_S\d{1,3}_R(1|2)_001(\.fastq\.gz)$'

#############################################################################

## CUSTOM WARNING HANDLING ##################################################

# as I work through issues, I can silence certain warnings that I've since addressed, but will want to warn me in future
suppress_warnings = {
    'duplicated_ctabs': False,       # duplicate CTAB codes found
    'missing_coll_date': False,      # missing collection date found, fixed following warning
    'nonuniform_coll_date': False,  # a subplot has various collection dates (YYYY-MM) associated with it
    'multiple_domains': False,      # a state abbreviation has multiple NEON domains associated with it
    'missing_old_site': False,      # a NEON domain doesn't have an old climush site name associated with it
    'unrecog_old_site': False,      # an old site name from a sample ID in the dataframe isn't recognized within NEON convert dict
}

#############################################################################

########################################################################################################################


## SET FILE PATHS ######################################################################################################

## OUTPUT ###################################################################

## SUMMARY OUTPUT DIRECTORY ##

# create a file and directory prefix for output files based on the endophyte sampling year
mapping_file_yr = re.search('2023|2022', args['mapping'].name).group(0)
output_file_prefix = f'{mapping_file_yr}-endophytes'

# if no output directory is provided via the command line, create the default
if args['output'] is None:

    # the default summary directory will be at the same level as the input directory
    args['output'] = args['input'].parent / f'{output_file_prefix}_rename-summary'

    # after creating path, make this directory based on the path
    args['output'].mkdir(exist_ok=True)

else:
    pass


## UNLINKED CTAB CODE TABLE ##

# name for the .csv file containing a filtered df with only CTAB codes that don't match a file name
date_today = datetime.today().strftime('%Y-%m-%d')
unlinked_ctab_path = args['output'] / f'{output_file_prefix}_unlinked-ctab-codes_{date_today}.csv'


## ORIGINAL SEQUENCE FILE COPIES ##

# define a path for copy of original sequence files, unless --no-copy is specified
if args['no_copy']:
    pass
else:
    original_copy_path = args['input'] / 'original-files'
    original_copy_path.mkdir(exist_ok=True)


## FILE NAME CONVERSION TABLE ##

# name for the .csv file that has the original file names and their associated new file names
convert_table_path = args['output'] / f'{output_file_prefix}_file-name-conversion_{date_today}.csv'

#############################################################################

########################################################################################################################


# ## CREATE TEST FILES ###################################################################################################
#
# # only make a dummy directory and dummy files if the --dry-run flag is used
# if args['dry_run']:
#
#     # create a parent directory in which to make the dry run files
#     file_rename_path = args['input'].parent / f'{output_file_prefix}_dry-run-files'
#     file_rename_path.mkdir(exist_ok=True)
#
#     # look for a file list in the input path parent directory that matches the current sample year (here, 2023)
#     test_file_list = [f for f in args['input'].parent.glob(f'*{mapping_file_yr}*.txt')]
#
#     # if only a single file was found, then create the dummy files
#     if len(test_file_list) == 1:
#
#         # go through each file name in the demux endophyte Illumina seq file names, and create an empty file for each
#         dummy_file_paths = []
#         with open(test_file_list[0], 'rt') as fin:
#             for f_name in fin.readlines():
#                 # only add paths if the file is a sequence file (fasta/fastq/fastq.gz)
#                 if re.search(SEQ_FILE_RE, f_name) or re.search(GZIP_REGEX, f_name):
#                     dummy_file_paths.append(file_rename_path / f_name.strip())
#                 else:
#                     pass
#
#         # go through each of these file paths and create empty files for each
#         for dummy_file in dummy_file_paths:
#             with open(dummy_file, 'w') as fout:
#                 fout.write('')
#
#     elif len(test_file_list) == 0:
#
#         # error; no files were found
#         print(f'error. no files were automatically detected that might be the sequence file name .txt file.')
#         sys.exit()
#
#     else:
#
#         # error; multiple possible files found
#         print(f'error. {len(test_file_list)} files were automatically detected that might be the sequence file name .txt file:\n'
#               f'   {test_file_list}')
#         sys.exit()
#
# else:
#     pass
#
# ########################################################################################################################


## READ IN SAMPLE MAPPING FILE #########################################################################################

# read in .csv file; 2023 mapping file is a single tab
endo_sampling_df = pd.read_csv(args['mapping'])

########################################################################################################################


## CLEAN ENDOPHYTE MAPPING FILE ########################################################################################

# take only the useful columns in this dataframe (need to look to know what column names are)
endo_df_original = endo_sampling_df[['Code', 'BAG LABEL']]

# make a copy of the original input dataframe before making changes
endo_df = endo_df_original.copy()

# rename the column names to be more clear
endo_df.rename(columns={'Code': 'sample_id_az',
                        'BAG LABEL': 'sample_id_or'},
               inplace=True)

# fix any typos with the CTAB code prefix, which should be MSC (all here are MCSC, which doesn't match seq files)
endo_df['sample_id_az'] = [re.sub(ANY_CTAB_PREFIX, CORRECT_CTAB_PREFIX, c) for c in endo_df['sample_id_az']]

# strip any white space across any sample ID column
for col in endo_df.columns:

    # look for word 'sample' to identify if this column is a sample column
    sample_col_found = re.search(r'sample', col, re.I)

    # if this column is a sample column, then strip white space from all row values
    if sample_col_found:
        endo_df[col] = endo_df[col].apply(lambda row: row.strip())

    # if this column wasn't determined to be a sample ID column, pass over it
    else:
        continue

# check to see if any of the CTAB codes are duplicated
duplicate_ctab = endo_df[endo_df['sample_id_az'].duplicated(keep=False)].sort_values('sample_id_az')  # keep=False to return all duplicate rows

# raise an error if there are duplicate CTAB codes found (unless warning is suppressed)
if (duplicate_ctab.shape[0] == 0) or (suppress_warnings['duplicated_ctabs']):
    pass
else:
    raise KeyboardInterrupt(
        f'{duplicate_ctab.shape[0]} CTAB IDs are duplicated ({duplicate_ctab.shape[0] / 2} pairs):\n'
        f'   {set(duplicate_ctab["sample_id_az"].to_list())}')

########################################################################################################################


## SEQUENCE FILES MISSING FROM TABLE ###################################################################################

# get a list of all files in the input directory that contain the CTAB prefix MSC
ctab_file_names = [file.name for file in args['input'].glob(SEQ_FILE_GLOB) if re.search(CTAB_RGLOB_RE, file.name, re.I)]
nonctab_file_names = [file.name for file in args['input'].glob(SEQ_FILE_GLOB) if not file.name in ctab_file_names]
nonctab_file_names.sort()

# quick check that every file has a pair by making sure number of files is even
if (len(ctab_file_names) % 2 == 0) and (len(nonctab_file_names) % 2 == 0):
    pass
elif (len(ctab_file_names) % 2 == 0):
    print(f'Non-CTAB file names do not have an even number of files.')
elif (len(nonctab_file_names) % 2 == 0):
    print(f'CTAB file names do not have an even number of files.')
else:
    print(f'Neither CTAB nor non-CTAB file names have an even number of files.')

# print counts for non- and ctab file names
print(f'number of paired sequence files:\n'
      f'   CTAB files     = {len(ctab_file_names)} ({len(ctab_file_names)/2:.0f} samples)\n'
      f'   non-CTAB files = {len(nonctab_file_names)}   ({len(nonctab_file_names)/2:.0f} samples)\n')

# pull out just the CTAB codes from the ctab file names
seq_file_ctabs = [re.search(ANY_CTAB_CODE, filename).group(0) for filename in ctab_file_names]

########################################################################################################################


## GET PER-SITE BAG LABEL EXAMPLE ######################################################################################

# create a set of the unique site prefixes used for the bag labels
unique_site_labels = set(label.split('-')[0] for label in endo_df['sample_id_or'])

# for each of the unique site prefixes, pull a random bag label for each to use as an example for renaming
per_site_label_example = {site:'' for site in unique_site_labels}
for site in per_site_label_example:
    per_site_label_example.update({site:endo_df[endo_df['sample_id_or'].str.startswith(site)]['sample_id_or'].iloc[0]})
    print(f'{per_site_label_example[site]}')

########################################################################################################################


## INVESTIGATE PER-SITE LABEL FORMATS ##################################################################################

# get the position and label formats for each site
with open(args['config'], 'rb') as fin:
    rename_settings = tomlkit.loads(fin.read())
    rename_settings.remove('title')

# look at the different labels used within a site
def unique_labels(label_list, site, label_position, label_delim=args['delim']):

    # create a set of unique labels encountered in the provided position of the bag label
    unique_label_set = set()

    # create a list to add any bag labels that don't follow the format of other bag labels for this site
    label_errors = []

    # iterate through each bag label provided...
    for label in label_list:

        # only add the labels that match the site provided
        if label.startswith(site):

            # try to find this label in the provided position
            try:
                unique_label_set.add(label.split(label_delim)[label_position])
            # if the index is out of bounds, then this bag label likely doesn't follow the site's naming convention
            except IndexError:
                # add the full bag label to the list of bag label errors
                label_errors.append(label)
        else:
            continue

    return unique_label_set, label_errors

# test on BNZ first; used when creating the .toml configuration file
unique_labels(endo_df['sample_id_or'], site='BNZ', label_position=4)

# iterate through the config file to see the format used for the treatment label
for site in rename_settings.keys():

    # print the site being viewed
    print(f'\n{site}')

    # what are the unique labels for treatment - burn history?
    print('  burn history:')
    unique_burnhist, errors_burnhist = unique_labels(endo_df['sample_id_or'], site=site,
                                                     label_position=rename_settings[site]['treatment']['burnhist']['position'])
    print(f'    labels = {unique_burnhist}\n'
          f'    errors = {errors_burnhist}')

    # what are the unique labels for treatment - habitat?
    print(f'  habitat:')
    unique_habitat, errors_habitat = unique_labels(endo_df['sample_id_or'], site=site,
                                                   label_position=rename_settings[site]['treatment']['habitat']['position'])
    print(f'    labels = {unique_habitat}\n'
          f'    errors = {errors_habitat}')

########################################################################################################################


## UPDATE OLD CLIMUSH SAMPLE ID ########################################################################################

# create a list to add the new column names (i.e., updated file name components columns); easier when combining at end
update_filename_cols = []

# add a column of repeat values for the gene region; dumb but easiest solution right now
endo_df['sequence_region'] = seq_region_prefix
update_filename_cols.append('sequence_region')


## COMPARTMENT #########################################################################################################

## ISOLATION SOURCE #########################################################

# leaf (L), root, or seed

# what type of isolation sources are in this group of samples
isolate_sources = set(s.split('-')[-1] for s in endo_df['sample_id_or'])
# print(f'isolation sources in dataset: {isolate_sources}')

# get the isolation source, convert to full string from single letter
isolate_source_dict = {
    'L': 'leaf',
    'S': 'seed',  # assuming this part; no seeds in test data
    'R': 'root',
}

# make sure that all isolate sources in df are accounted for in isolate_source_dict
for iso_src in isolate_sources:
    if iso_src in isolate_source_dict:
        pass
    else:
        # no suppression for this warning because it'll cause error below and no 'placeholder' solution (or I guess I could but meh)
        raise KeyboardInterrupt(f'Isolation source {iso_src} is not an expected value.')

# will do conversion at same time as isolation species so I only have to loop through items once in this section

#############################################################################

## ISOLATION SPECIES ########################################################

# this I'll handle when looping through, nothing needs to be done prior

#############################################################################

isolation_source_col = []
for old_id in endo_df['sample_id_or']:

    # get the old species and old isolation source from the old sample ID
    old_sp, old_src = old_id.split('-')[-2:]

    # fix species from Sp.1 to sp01; add leading zero
    new_sp = ''.join([old_sp.split('.')[0].lower(), '0', old_sp.split('.')[-1]])

    # spell out source from single letter
    new_src = isolate_source_dict[old_src]

    # join together, with - delim
    new_isol = '-'.join([new_src, new_sp])

    # add this to the updated list of isolation sources
    isolation_source_col.append(new_isol)

# update the df with this new column
endo_df['isolation_source'] = isolation_source_col
update_filename_cols.append('isolation_source')

########################################################################################################################

### COLLECTION DATE ####################################################################################################

# string format the collection date column, so that it is formatted correctly for the new naming convention
formatted_dates = []
date_found_counter = 0
missing_date_rows = []
for i,date in endo_df['collection_date'].items():  # iterate through row index and row value for this column
    if isinstance(date, datetime):
        date_found_counter += 1
        formatted_dates.append(date.strftime('%Y-%m'))
    else:
        # add row index to a list; will use to determine what val is in this column that isn't a datetime val
        missing_date_rows.append(i)
        # add 'nan' to the list, so it will match same dim as the df
        formatted_dates.append(np.nan)

# formatted dates is the full length of df, date found counter is how many are actual dates
date_count_diff = len(formatted_dates) - date_found_counter
if (date_count_diff == 0) or (suppress_warnings['missing_coll_date']):
    pass
else:
    print(f'Total number of rows checked:     {len(formatted_dates)}\n'
          f'Number of datetime objects found: {date_found_counter}')
    # use .loc, not .iloc, because row indices have not been reset so don't reflect actual positions
    missing_date_df = endo_df.loc[missing_date_rows,:]
    warning_msg = (f'{date_count_diff} row(s) missing a collection date:\n'
                   f'       {missing_date_df}')
    raise KeyboardInterrupt(warning_msg)

# replace the current collection_date column with the formatted collection dates
# do this prior to finding the missing collection date so you can pull the reformatted date to use
endo_df['collection_date'] = formatted_dates
update_filename_cols.append('collection_date')

# for the missing date value, look at other samples from this subplot and use their collection date
missing_date_sample = endo_df.loc[missing_date_rows,:]['sample_id_or'].iloc[0]  # get the sample ID of missing date sample
missing_subplot_info = '-'.join(missing_date_sample.split('-')[:4])  # get just the subplot info from this missing date sample

# make sure to exclude the problematic samples with na in the collection date, too
common_subplot_colldates = endo_df[endo_df['sample_id_or'].str.match(f'^{missing_subplot_info}') & endo_df['collection_date'].notna()]['collection_date']

# TRIED TO EXTRACT EXACT VALUES FOR FLEXIBLE USE IN FUTURE BUT EXTRACT NAN != IN-TABLE NAN (BOTH FLOATS?)
# missing_current_date = endo_df.loc[missing_date_rows,:]['collection_date'].iloc[0]  # get the current value used for this date (usually NaN)
# common_subplot_colldates = endo_df[endo_df['sample_id_or'].str.match(f'^{missing_subplot_info}') & endo_df['collection_date'] != missing_current_date]['collection_date']


# check if they all have the same coll dates
if len(set(common_subplot_colldates)) == 1:
    # if they're all the same, use this to populate the collection_date column for this sample
    # endo_df = endo_df.replace({'collection_date': missing_current_date}, common_subplot_colldates.iloc[0])
    endo_df = endo_df.replace({'collection_date': pd.NA}, common_subplot_colldates.iloc[0])
elif suppress_warnings['nonuniform_coll_date']:
    pass
else:
    # if they aren't all the same, that's weird and figure out why
    subplot_df = endo_df[endo_df['sample_id_or'].str.match(f'^{missing_subplot_info}')]
    raise KeyboardInterrupt(f'Collection dates for this subplots are not all uniform:\n'
                            f'{subplot_df}')


########################################################################################################################


### ECOREGION ##########################################################################################################

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

neon_domains = {'AK': {'D19': set()},
                'AZ': {'D14': set()},
                'CO': {'D13': set()},
                'FL': {'D03': set()},
                'KS': {'D06': set()},
                'MA': {'D01': set()},
                'MN': {'D05': set()},
                'OR': {'D16': set()}}

update_domain_dict(neon_domains, configuration_dict=rename_settings)

# unnest the neon_domains dict, since the state abbreviation isn't used anywhere
neon_domain_convert = {}
for state_abbr in neon_domains:

    # get the neon domain
    n_dom_list = list(neon_domains[state_abbr].keys())

    # confirm only one domain associated with this state (should always be 1)
    if (len(n_dom_list) > 1) or (suppress_warnings['multiple_domains']):
        raise KeyboardInterrupt(f'The state {state_abbr} has multiple NEON domains associated with it.')
    else:
        n_dom = n_dom_list[0]

    # get the old site name or names (these are sets, convert to list)
    old_site_list = list(neon_domains[state_abbr][n_dom])

    # if there are multiple sites for this domain, update the convert dict with multiple entries
    if len(old_site_list) > 1:
        for site in old_site_list:
            neon_domain_convert.update({site: n_dom})
    elif len(old_site_list) == 1:
        neon_domain_convert.update({old_site_list[0]: n_dom})
    elif suppress_warnings['missing_old_site']:
        pass
    else:
        raise KeyboardInterrupt(f'No old CliMush site names were recovered for the NEON domain {n_dom}.')

# use the neon domains to replace the old site names with the neon domains
neon_domain_col = []
for old_id in endo_df['sample_id_or']:

    # get the site name from the old sample ID
    old_site = old_id.split('-')[0]

    # if the site is in the conversion dict (should be) use it to find the NEON domain, add to list
    if old_site in neon_domain_convert:
        neon_domain_col.append(neon_domain_convert[old_site])

    # if the site is not in the convert dict...
    else:
        # check if its because its SRE, which is written SRER in mapping df
        if old_site == 'SRER':
            neon_domain_col.append(neon_domain_convert['SRE'])

        # otherwise, it's some unrecognized old site name
        else:
            if suppress_warnings['unrecog_old_site']:
                # if warning suppressed, add a placeholder value to the list
                neon_domain_col.append(pd.NA)
            else:
                raise KeyboardInterrupt(f'The old site name {old_site} was not recognized.')

# add the neon domains as a column in the mapping df
endo_df['neon_domain'] = neon_domain_col
update_filename_cols.append('neon_domain')

########################################################################################################################


## TREATMENT ###########################################################################################################

def get_treatment(bag_label, config_dict, label_delim=args['delim']):

    # create a list of the original bag label components
    label_list = bag_label.split(label_delim)

    # get the site ID of this bag label, which is always the first part of the bag label for the 2023 endophytes
    site_id = label_list[0]

    # find the original labels
    og_burnhist_label = label_list[config_dict[site_id]['treatment']['burnhist']['position']]
    og_habitat_label = label_list[config_dict[site_id]['treatment']['habitat']['position']]

    # translate the original labels to new labels
    new_burnhist_label = config_dict[site_id]['treatment']['burnhist']['format'][og_burnhist_label]
    new_habitat_label = config_dict[site_id]['treatment']['habitat']['format'][og_habitat_label]

    # return the combined burn history and habitat the treatment labels
    return new_burnhist_label + new_habitat_label

# test function

new_treatment_dict = {}
error_labels = []
for old_label in endo_df['sample_id_or']:

    # skip over SRER and BNZ
    if not (old_label.startswith('SRER') or old_label.startswith('BNZ')):

        try:
            # get the new label by using the get_treatment function
            new_label = get_treatment(bag_label=old_label, config_dict=rename_settings)

            # update the dictionary with the old label as the key and the new label as the value
            new_treatment_dict.update({old_label: new_label})

        except tomlkit.exceptions.NonExistentKey:
            error_labels.append(old_label)


########################################################################################################################


## SUBPLOT #############################################################################################################

# what are the unique subplots in old sample names?
old_subplots = list(set(s.split('-')[-3] for s in endo_df['sample_id_or']))  # make list, easier to check if sorted
old_subplots.sort()
# print(f'old subplot labels = {old_subplots}')

# everything looks like I expected it to, so do simple pass through to get number and format w/ leading zero
subplot_col = []
for samp in endo_df['sample_id_or']:
    # get just the number from the old subplot label
    old_subplot_num = samp.split('-')[-3].replace('S','')  # use replace in case there's a 2-digit subplot num

    # only add a leading zero if its a single-digit number
    if int(old_subplot_num) < 10:
        subplot_col.append('0' + old_subplot_num)
    else:
        subplot_col.append(old_subplot_num)

# add formatted subplot values to a new column
endo_df['subplot'] = subplot_col
update_filename_cols.append('subplot')

########################################################################################################################


## CREATE NEW SAMPLE IDS ###############################################################################################

endo_df['updated_sample_id'] = endo_df[update_filename_cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

########################################################################################################################


## RENAME SEQUENCE FILES ###############################################################################################

## LINK CTAB CODE TO UPDATED CLIMUSH ID #####################################

# create a conversion dictionary of old and new file names (no file extensions, just IDs)
file_rename_dict = {old_name: new_name for old_name, new_name in zip(endo_df['sample_id_az'], endo_df['updated_sample_id'])}

#############################################################################

## LINK OLD AND NEW FILE PATHS ##############################################

# create a set to add file paths to when the file isn't a sample file (e.g., a mock community file)
#   use sets because there will be two files per file base name (r1/r2)
nonsample_files = set()

# create a set to add file paths to when the file name starts with the CTAB prefix (i.e., MSC) but doesn't match
#   a CTAB code in from the mapping dataframe
unknown_sample_files = set()

# create a dictionary to add old file paths and new file paths as key/value pairs; to be used for renaming files
filepath_rename_dict = {}

# keep track of which CTAB codes are linked to a file path; will want to check that no CTAB codes were unlinked
matched_ctab_codes = []

# go through each of the old file paths, creating new file paths (w/ new file names) for each
for filepath_old in file_rename_path.glob(SEQ_FILE_GLOB):

    # get the file stem (i.e., file name without the file extension)
    filename_old = re.sub(r'-R(1|2)', '', strip_file_ext(filepath_old))
    # get the file extension (suffix) to use for new file path
    filesuffix_old = filepath_old.suffix
    # get the read orientation of the file to add to new file path
    read_orient = re.search(r'R(1|2)', strip_file_ext(filepath_old)).group(0)

    # if the old file name is in the file rename dict...
    if filename_old in file_rename_dict:

        # get the new file name; add the read-orient to the end of it (R1/R2)
        filename_new = '_'.join([file_rename_dict[filename_old], read_orient])

        # create a new file path with the new file name, using the same parent directory path and file extension as the old filename
        filepath_new = (filepath_old.parent / filename_new).with_suffix(filesuffix_old)

        # add the old file path and new file path as key/value pair in the file path rename dict
        filepath_rename_dict.update({filepath_old: filepath_new})

        # add the file name (i.e., ctab code) to a list that tracks which CTAB codes from the df were matched to a file
        matched_ctab_codes.append(filename_old)

    # if the old file name isn't in the file name rename dict...
    else:

        # there are some actual samples that are there, but not in table?
        if filename_old.startswith(CORRECT_CTAB_PREFIX):
            unknown_sample_files.add(filename_old)

        # otherwise, it's likely mock communities and ctrl samples
        else:
            # add it to a list for now, will handle once I see what's there
            nonsample_files.add(filename_old)


#############################################################################


## RENAME SEQUENCE FILES ####################################################

# keep track of how many files were successfully renamed
rename_count = 0

# go through each old and new file path pair...
for old_path, new_path in filepath_rename_dict.items():

    # if using safe mode, copy this file to the original file copy directory
    if args['no_copy']:
        pass
    else:
        copy_dst = original_copy_path / old_path.name
        shutil.copy(old_path, copy_dst)

    # rename the old path to the new path
    old_path.rename(new_path)

    # only add to rename count if the new path is confirmed to be an existing file
    if new_path.is_file():
        rename_count += 1
    else:
        pass

print(f'{rename_count} total sequence files were successfully renamed.\n')

#############################################################################

## EXPORT FILE NAME CONVERSION TABLE ########################################

# get the name of the old and new files and create a conversion table for these pairs
original_file_names = [og_file.name for og_file in filepath_rename_dict.keys()]
new_file_names = [new_file.name for new_file in filepath_rename_dict.values()]
filename_convert_df = pd.DataFrame(data={'original_names': original_file_names,
                                         'new_names': new_file_names})

# write this table out to the conversion table path
filename_convert_df.to_csv(convert_table_path, index=False)

#############################################################################

########################################################################################################################


## ADDRESS CTAB ID AND FILENAME DISCREPENCIES ##########################################################################

# some errors don't arise until you can't link the CTAB code in the file name to the one in the mapping dataframe

# how many unknown sample files are there?
unknown_sample_files_list = list(unknown_sample_files)
unknown_sample_files_list.sort()
print(f'There are {len(unknown_sample_files)} unrecognized sample files: {unknown_sample_files_list}\n')

# check original version of the df for these samples (before anything was done to table)
og_samp_check = [code for code in endo_df_original['CTAB Codes'] if code in unknown_sample_files]
print(f'There are {len(og_samp_check)} of these unrecognized sample names are in the original table.\n')

# for any unrecognized sample files, strip te letters and compare just the numbers; often issue with letters
unknown_sample_no = [ re.search(r'(\d+)', s).group(0) for s in unknown_sample_files ]
ctab_sample_no = [ re.search(r'(\d+)', c).group(0) for c in endo_df['sample_id_az'] if pd.notna(c)]


## I THINK THIS IS NOW UNNECESSARY SINCE I UPDATED THE REGEX FOR CORRECTING CTAB PREFIXES ###########
# no_prefix_typo = []
# unknown_sample_list = list(unknown_sample_files)  # need to be able to subscript, can't w/ set
# caught_sample_errors = {  # add an errors that can automatically fixed to the dictionary
#     'dataframe': {},  # errors in the dataframe that need to be corrected
#     'filename': {},  # errors in the sequence file name that need to be corrected
# }
# for u, unknown_no in enumerate(unknown_sample_no):
#
#     # check if the number part of the seq sample name is in the CTAB ID of cleaned df
#     if unknown_no in ctab_sample_no:
#
#         # create a var; used to determine whether to print comparison (no fix found) or not (fix found)
#         fix_found = False
#
#         # if it is found, get the index of the matching number in the CTAB ID df list (corresponds to row num)
#         ctab_row_num = ctab_sample_no.index(unknown_no)
#
#         # use that index number to find the full CTAB ID whose number matches
#         ctab_id = endo_df['sample_id_az'].iloc[ctab_row_num]
#
#         # get the full seq file ID whose number matches
#         unknown_id = unknown_sample_list[u]
#
#         # compare if lengths of sample IDs are same (probably prefix error; otherwise some samples with extra A/B at end)
#         if len(unknown_id) == len(ctab_id):
#             ctab_prefix = re.search('^[A-Z]{3}', ctab_id).group(0)
#             unknown_prefix = re.search('^[A-Z]{3}', unknown_id).group(0)
#
#             # if they both have the correct CTAB prefix, same number, same length; then still confused on issue
#             if (ctab_prefix == CORRECT_CTAB_PREFIX) and (unknown_prefix == CORRECT_CTAB_PREFIX):
#                 no_prefix_typo.append(unknown_id)
#
#             # if only one of the has a correct CTAB prefix...
#             elif (ctab_prefix == CORRECT_CTAB_PREFIX) or (unknown_prefix == CORRECT_CTAB_PREFIX):
#
#                 # switch boolean flag to True; solution found, no need to show comparison
#                 fix_found = True
#
#                 # correct the CTAB prefix for the one that is NOT correct
#                 if ctab_prefix == CORRECT_CTAB_PREFIX:
#                     # correct the CTAB prefix, and add the error/fixed version as key/value pair in fixed dict
#                     unknown_id_fixed = re.sub(f'^{unknown_prefix}', CORRECT_CTAB_PREFIX, unknown_id)
#                     caught_sample_errors['filename'].update({unknown_id: unknown_id_fixed})
#                 else:
#                     # correct the filename prefix, and add the error/fixed version as key/value pair in fixed dict
#                     ctab_id_fixed = re.sub(f'^{ctab_prefix}', CORRECT_CTAB_PREFIX, ctab_id)
#                     caught_sample_errors['dataframe'].update({ctab_id: ctab_id_fixed})
#
#             # if neither of them have the correct CTAB prefix, but same number and length; then still confused on issue
#             else:
#                 no_prefix_typo.append(unknown_id)
#
#         if fix_found:
#             continue
#         else:
#             # show comparison if a solution wasn't found
#             print(f'The number of the following two samples IDs match but the prefix does not:\n'
#                   f'   CTAB label:  {ctab_id}\n'
#                   f'   sequence ID: {unknown_id}\n')
#
#     # if the number part of the seq sample name is still not located
#     else:
#         # add the full unknown sample ID into a 'still unknown' list
#         no_prefix_typo.append(unknown_sample_list[u])


########################################################################################################################


## UNLINKED CTAB CODES #################################################################################################

# which CTAB codes from the dataframe were not matched to a sequence file name?
unlinked_ctab_codes = set(endo_df['sample_id_az']).difference(set(matched_ctab_codes))
print(f'{len(unlinked_ctab_codes)} CTAB codes were not matched to a sequencing file.\n')

# check the other tabs to see if these CTAB codes are recorded elsewhere

# functions to use for searching dataframe tabs
# check for exact matches
def search_exact_matches(search_vals, dataframe, compare_column, return_matches=True):

    match_in_comparison = [c for c in search_vals if c in dataframe[compare_column]]
    if len(match_in_comparison) == 0:
        print(f'No matches were found in the {compare_column} column.\n')
        return None
    else:
        if return_matches:
            return match_in_comparison
        else:
            print(f'{len(match_in_comparison)} matches were located in the \'{compare_column}\' column.\n')
            return None

# ctab processing - details which samples were extracted for sequencing
ctab_process_df = endo_sampling_tabs['CTAB_processing']

# collection processing - prep of samples for culture and sequencing
# this tab has header so that first 6 rows are not column names, but 7 is
# no, i'm not fucking checking this, the header has completely fucked it up
# coll_process_df = endo_sampling_tabs['Collections_processing'].iloc[7:,:]

# collection processing df doesn't have a CTAB code (has col name but no values in this column); look up the old
#   climush sample ID in the mapping df, use these to search
unlinked_old_climush_names = endo_df[endo_df['sample_id_az'].isin(unlinked_ctab_codes)]['sample_id_or'].to_list()
if len(unlinked_old_climush_names) == len(unlinked_ctab_codes):
    pass
else:
    raise KeyboardInterrupt(f'The number of old climush samples IDs does not match the number of unlinked CTAB codes.\n')

# search for exact matches
search_exact_matches(unlinked_ctab_codes, ctab_process_df, 'CTAB Codes MSC####', return_matches=False)
# search_exact_matches(unlinked_old_climush_names, coll_process_df, 'BAG LABEL', return_matches=False)

# filter the original table to only include these unlinked CTAB samples and save to file
unlinked_ctab_df = endo_df[endo_df['sample_id_az'].isin(unlinked_ctab_codes)]
with open(unlinked_ctab_path, 'w') as csv_out:
    unlinked_ctab_df.to_csv(csv_out, index=True, index_label='original_df_row')


########################################################################################################################"