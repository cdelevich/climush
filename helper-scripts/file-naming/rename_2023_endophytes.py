from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import re, tomlkit, shutil
from climush.constants import SEQ_FILE_GLOB, GZIP_REGEX, SEQ_FILE_RE
from climush.utilities import strip_file_ext

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

# lots of errors with the CTAB code prefix, define var with correct one to compare
CORRECT_CTAB_PREFIX = r'MSC'
# regex for correcting any prefix error (captures many combinations of 'M', 'S', or 'C' in CTAB ID)
ANY_CTAB_PREFIX = r'^[MSC]+?(?=\d)'

#############################################################################

## CUSTOM WARNING HANDLING ##################################################

# as I work through issues, I can silence certain warnings that I've since addressed, but will want to warn me in future
suppress_warnings = {
    'duplicated_ctabs': True,       # duplicate CTAB codes found
    'missing_coll_date': True,      # missing collection date found, fixed following warning
    'nonuniform_coll_date': False,  # a subplot has various collection dates (YYYY-MM) associated with it
    'multiple_domains': False,      # a state abbreviation has multiple NEON domains associated with it
    'missing_old_site': False,      # a NEON domain doesn't have an old climush site name associated with it
    'unrecog_old_site': False,      # an old site name from a sample ID in the dataframe isn't recognized within NEON convert dict
}

#############################################################################


## SAFE MODE ################################################################

# whether to create a copy of the original sequencing files before renaming
MAKE_COPY = True

#############################################################################


########################################################################################################################


## SET FILE PATHS ######################################################################################################

## INPUT ####################################################################

# main path to the file naming helper scripts
file_naming_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/')

# path to the directory containing the sequencing files to rename; if doing an empty-file test, leave as empty string
file_rename_path = ''

# endophyte sampling dataframe from Google drive
sampling_main_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/renaming-mapping-files/from-neda-arnold-lab/')
endo_sampling_path = 'its1_endophyte-mapping_2023.csv'

# when testing issue with 2022 v 2023 mapping file, add prefix to missing ctab summary file
mapping_file_yr = re.search('2023|2022', endo_sampling_path.name).group(0)
mapping_file_prefix = f'{mapping_file_yr}-sample-sheet'
print(f'Sample sheet = {mapping_file_yr}')

# text file that contains the name of the files that need to be renamed; one file name per line
# test_name_path = file_naming_path / 'illumina_endophytes_2022_file-list.txt'
test_name_path = file_naming_path / 'its1-endophytes-2023-1_original-file-names.txt'

# when testing issue with 2022 v 2023 mapping file, add prefix to missing ctab summary file
test_file_yr = re.search('2023|2022', test_name_path.name).group(0)
test_file_prefix = f'{test_file_yr}-seq-files'
test_file_dir = f'test-files_{test_file_yr}-sequence-files'

# file renaming configuration file; for whatever reason, the endo toml file wouldn't load correctly
config_main_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/')
rename_settings_path = config_main_path / 'illumina_endos_2023-10_file-rename-config.toml'

## TEST INPUT PATHS #########################################################

# only make a dummy directory and dummy files if there is no value provided for the file_rename_path variable
if file_rename_path == '':

    # print which group of sequence files are being used to created these empty test files
    print(f'Sequence files = {test_file_yr}\n')

    # create a parent directory in which to make the empty dummy files
    file_rename_path = file_naming_path / test_file_dir
    file_rename_path.mkdir(exist_ok=True)

    # go through each file name in the demux endophyte Illumina seq file names, and create an empty file for each
    dummy_file_paths = []
    with open(test_name_path, 'rt') as fin:
        for f_name in fin.readlines():
            # only add paths if the file is a sequence file (fasta/fastq/fastq.gz)
            if re.search(SEQ_FILE_RE, f_name) or re.search(GZIP_REGEX, f_name):
                dummy_file_paths.append(file_rename_path / f_name.strip())
            else:
                pass

    # go through each of these file paths and create empty files for each
    for dummy_file in dummy_file_paths:
        with open(dummy_file, 'w') as fout:
            fout.write('')

else:
    pass

#############################################################################

## OUTPUT ###################################################################

# create a path to write out any summary information for this test
summary_info_path = file_naming_path / 'endophyte-renaming-summary'
summary_info_path.mkdir(exist_ok=True)

# name for the .csv file containing a filtered df with only CTAB codes that don't match a file name
date_today = datetime.today().strftime('%Y-%m-%d')
unlinked_ctab_path = summary_info_path / f'unlinked-ctab-codes_{test_file_prefix}_{mapping_file_prefix}_{date_today}.csv'

# define a path for copy of original sequence files if using safe mode
if MAKE_COPY:
    original_copy_path = file_rename_path / 'original-files'
    original_copy_path.mkdir(exist_ok=True)
else:
    pass

########################################################################################################################


## READ IN SAMPLE MAPPING FILE #########################################################################################

# read in each tab of the excel file; can do this by setting sheet_name=None
# this will output a dictionary, where the key is the name of the tab and the value is the dataframe within this tab
endo_sampling_tabs = pd.read_excel(endo_sampling_path, sheet_name=None)

# figure out which tab has the CTAB codes in them

# look for a column name containing the string 'CTAB codes' (w/ or w/o 's') and that has values in its rows
for tab in endo_sampling_tabs:
    tab_table = endo_sampling_tabs[tab]
    for col in tab_table.columns:
        if re.search(r'ctab code', col, re.I):
            row_vals = set(tab_table[col].to_list())
            if len(row_vals) > 1:
                endo_sampling_df = tab_table
                break
            else:
                continue
        else:
            continue

########################################################################################################################


## CLEAN ENDOPHYTE MAPPING FILE ########################################################################################

# take only the useful columns in this dataframe (need to look to know what column names are)
endo_df_original = endo_sampling_df[['CTAB Codes', 'Sample code', 'Collection month']]

# drop any rows where there is no CTAB code
endo_df = endo_df_original[~endo_df_original['CTAB Codes'].isna()]

# drop the one sample that doesn't have a CTAB Code, instead has '?'
endo_df = endo_df[endo_df['CTAB Codes'] != '?']

# rename the column names to be more clear
endo_df.rename(columns={'CTAB Codes': 'sample_id_az',
                        'Sample code': 'sample_id_or',
                        'Collection month': 'collection_date'},
               inplace=True)

# fix any typos with the CTAB code prefix, which should be MSC
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

# duplicated CTAB codes were found; after looking in CTAB processing tab, I found note to explain this; addressing here

# iterate through the row indices of the original dataframe for these duplicates
for d in duplicate_ctab.index:

    old_ctab_id = endo_df['sample_id_az'].loc[d]

    # I only have a solution from the notes for the 007# sample numbers
    if re.search(r'7\d$', old_ctab_id):

        # if the old climush ID starts with KON, add an 'A' suffix to it in the full table
        if re.search(r'^KON', endo_df['sample_id_or'].loc[d], re.I):

            endo_df.loc[d, 'sample_id_az'] = old_ctab_id + 'A0'

        # if the old climush ID does not start with KON (i.e., BNZ sample), do nothing to ctab ID
        else:
            continue

    # if not a 007# sample number, I don't have solution
    else:
        continue

# re-check for any remaining duplicates and to ensure the fixes above were properly implemented
duplicate_ctab = endo_df[endo_df['sample_id_az'].duplicated(keep=False)].sort_values('sample_id_az')  # keep=False to return all duplicate rows

if (duplicate_ctab.shape[0] == 0) or (suppress_warnings['duplicated_ctabs']):
    pass
else:
    raise KeyboardInterrupt(
        f'{duplicate_ctab.shape[0]} CTAB IDs are duplicated ({duplicate_ctab.shape[0] / 2} pairs):\n'
        f'   {set(duplicate_ctab["sample_id_az"].to_list())}')



# one sample later found to lack a collection date; also missing info elsewhere; check if it has a seq file
problem_sample01 = 'MSC0474'
problem_sample_found = False
for seq_file in file_rename_path.glob(SEQ_FILE_GLOB):
    if seq_file.name.startswith(problem_sample01):
        problem_sample_found = True
        problem_sample_name = seq_file.name
        break
    else:
        continue

search_msg = f'The sample ID {problem_sample01} '
if problem_sample_found:
    search_msg += (f'was linked to the following sequence file:\n'
                   f'   {problem_sample_name}\n')
else:
    search_msg += f'was not linked to a sequencing file.\n'
print(search_msg)

# because this sample has a sequencing file, don't drop it or anything; keep in df

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

# get the name of the sites used for each state from the renaming .toml doc
# if you have any empty values for a key, this won't load; so either put empty string or delete unused sections
with open(rename_settings_path, 'rb') as fin:
    rename_settings = tomlkit.loads(fin.read())

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

## HABITAT TYPE #############################################################

# what are the unique habitat types in old sample names?
old_habitats = set(s.split('-')[1] for s in endo_df['sample_id_or'])
# print(f'old habitat labels = {old_habitats}')

# check if any unaccepted values for habitat
accepted_habitats = ['C', 'O', 'G']
error_habitat = [h for h in old_habitats if not h in accepted_habitats]

def check_treatment_errors(expected_vals, positions, sample_col):

    # what values are observed in the old sample IDs?
    old_obs = set(s.split('-')[positions[0]] for s in sample_col)

    # are any old values that were observed *not* in the expected sample IDs?
    old_err = [h for h in old_obs if not h in expected_vals]

    # if errors found, look at which sample IDs cause this
    if len(old_err) > 0:

        # list of samples with non-accepted habitats
        err_samples = [s for s in sample_col if s.split('-')[positions[0]] in old_err]

        # check if the burn history part of the sample ID contains an accepted habitat
        remaining_err_samples = []
        swapped_position = []
        for err_samp in err_samples:
            # check for an accepted value in a possible location provided in function
            possible_val = err_samp.split('-')[positions[1]]
            if possible_val in expected_vals:
                swapped_position.append(err_samp)
            else:
                remaining_err_samples.append(err_samp)

        return remaining_err_samples, swapped_position

    else:
        return None


hab_err, hab_swap = check_treatment_errors(expected_vals=accepted_habitats,
                                           positions=[1,2],
                                           sample_col=endo_df['sample_id_or'])

def correct_treatment_errors(expected_vals, positions, sample_col):

    corrected_treatment = []
    uncorrected_treatment = []
    for samp in sample_col:

        # get the habitat (or what should be habitat)
        treatment = samp.split('-')[positions[0]]

        # if this is an accepted habitat, add to the list
        if treatment in expected_vals:
            corrected_treatment.append(treatment)

        # if this is not an accepted habitat...
        else:

            # address a specific issue for habitats in NWT samples
            if (positions[0] == 1) and (treatment == 'CO'):
                # if its the NWT samples that have CO instead of C, add C to corrected list
                corrected_treatment.append('C')
            else:
                treatment_swapped = samp.split('-')[positions[1]]
                if treatment_swapped in expected_vals:
                    corrected_treatment.append(treatment_swapped)
                # all errors should be caught by now?
                else:
                    print(f'Unknown error encountered with treatment: {treatment}')
                    uncorrected_treatment.append(treatment)

    # make sure there's the right number of habitats
    if len(corrected_treatment) == len(sample_col):
        return corrected_treatment
    else:
        print(f'The number of corrected habitats does not match the shape of the dataframe.')
        return None

habitats_corrected = correct_treatment_errors(expected_vals=accepted_habitats,
                                              positions=[1,2],
                                              sample_col=endo_df['sample_id_or'])

#############################################################################


## BURN HISTORY #############################################################

# what are the unique burn history types in old sample names?
old_burnhist = set(s.split('-')[2] for s in endo_df['sample_id_or'])
# print(f'old burn history labels = {old_burnhist}')

# check if any unexpected values for burn history
accepted_old_burnhist = ['B', 'UB']
error_burnhist = [bh for bh in old_burnhist if not bh in accepted_old_burnhist]

brn_err, brn_swap = check_treatment_errors(expected_vals=accepted_old_burnhist,
                                           positions=[2,1],
                                           sample_col=endo_df['sample_id_or'])

# no remaining burn errors after checking for position swaps, only position swap errors
burnhist_corrected = correct_treatment_errors(expected_vals=accepted_old_burnhist,
                                              positions=[2,1],
                                              sample_col=endo_df['sample_id_or'])

# burnhist_corrected is NOT corrected for renaming, just corrected to have consistent old names
# still need to change to updated burn history format (do when combining)

#############################################################################


## COMBINE HABITAT + BURN HISTORY ###########################################

# habitat format doesn't change in updated sample ID, so take as it is
# burn history does change, from B/UB to B/U; i.e., just use first letter only
treatment_col = [(brn[0] + hab) for hab,brn in zip(habitats_corrected, burnhist_corrected)]

# add treatment as new column in df
endo_df['treatment'] = treatment_col
update_filename_cols.append('treatment')

#############################################################################

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
    if MAKE_COPY:
        copy_dst = original_copy_path / old_path.name
        shutil.copy(old_path, copy_dst)
    else:
        pass

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