from pathlib import Path
import pandas as pd
from datetime import datetime

pd.set_option('display.max_columns', 10)
pd.set_option('display.max_rows', 10)

## FILE PATHS ##########################################################################################################

# PROJECT PATHS #

# path to the main directory for the ITS variant analysis
itsvar_main_path = Path('/Users/carolyndelevich/main/github_repos/climush/data-analysis/its-variants/')

# path to the data/ subdirectory in the ITS variant analysis directory
itsvar_data_path = itsvar_main_path / 'data'

# INPUT PATHS #

# path to the sequence identification verification table downloaded from Google Drive on 2025-06-22
verified_ids_path = itsvar_data_path / 'pacbio_sporocarps_2023a_genbank-unite-results_2025-06-22.csv'

# path to the table with the number of ITS variants for some sporocarps downloaded form Google Drive on 2025-06-22
itsvar_count_path = itsvar_data_path / 'PacBio All Runs iNat info_2025-06-22.csv'

# OUTPUT PATHS #

# output path for the cleaned file output produced by this script
date_suffix = datetime.now().strftime('%Y-%m-%d')
cleaned_df_out_path = itsvar_data_path / f'pacbio_its-variants_{date_suffix}.csv'

########################################################################################################################


## IMPORT DATA #########################################################################################################

# import the verified sequence ID table
verified_ids_original = pd.read_csv(verified_ids_path)

# import the ITS variant count table
itsvar_counts_original = pd.read_csv(itsvar_count_path)

# create a copy of each original dataframe, to which subsequent filtering and cleaning will be carried out
verified_ids_cleaned = verified_ids_original.copy()
itsvar_counts_cleaned = itsvar_counts_original.copy()

########################################################################################################################


## CLEAN + FILTER COLUMNS ##############################################################################################


# STRIP WHITESPACE FROM COLUMN LABELS #

# remove any leading or trailing whitespace from the column labels of the dataframes
verified_ids_colstrip = {col:col.strip() for col in verified_ids_cleaned}
itsvar_counts_colstrip = {col:col.strip() for col in itsvar_counts_cleaned}

# rename columns with the stripped column labels
verified_ids_cleaned.rename(columns=verified_ids_colstrip, inplace=True)
itsvar_counts_cleaned.rename(columns=itsvar_counts_colstrip, inplace=True)


# KEEP ONLY USEFUL COLUMNS #

# create a list of the columns from each dataframe that I want to keep
verified_ids_cols = [
    'Sample#', 'climush_plot', 'genbank_accession',
    'notes for Carolyn', 'fullITS_seq',
]
itsvar_counts_cols = [
    'sample_id', 'Nr dominant variants', 'Nr Copies of Dom Var',
    'scientific_name red means verified','ITS sequence', 'ITS sequence #2', 'LSU sequence',
    'GenBank Accession ITS', 'GenBank Accession LSU',
]

# filter dataframe columns
verified_ids_cleaned = verified_ids_cleaned[verified_ids_cols]
itsvar_counts_cleaned = itsvar_counts_cleaned[itsvar_counts_cols]


# RENAME REMAINING COLUMNS #

# rename the column containing the sample ID so the two dataframes match; they will later be merged
verified_ids_cleaned.rename(columns={'Sample#':'sample_id'}, inplace=True)

# rename species ID column
itsvar_counts_cleaned.rename(
    columns={
        'Nr dominant variants': 'its_variant_count',
        'Nr Copies of Dom Var': 'its_variant_top_copies',
        'scientific_name red means verified': 'identification',
        'ITS sequence': 'sequence_its_01',
        'ITS sequence #2': 'sequence_its_02',
        'LSU sequence': 'sequence_lsu',
        'GenBank Accession ITS': 'genbank_ref_its',
    },
    inplace=True,
)


# DROP COLUMNS W/O ANY VALUES #

# which columns in the verified IDs dataframe have only NaN rows?
verified_ids_allnull = [ col for col in verified_ids_cleaned.columns if all(verified_ids_cleaned[col].isna())]
itsvar_counts_allnull = [ col for col in itsvar_counts_cleaned.columns if all(itsvar_counts_cleaned[col].isna())]

# drop any columns that exclusively have NaN values
print(f'Columns removed because they only had null (NaN) values:')

# verified IDs dataframe
if len(verified_ids_allnull) > 0:
    verified_ids_cleaned.drop(axis=1, columns=verified_ids_allnull, inplace=True)
    print(f'   verified IDs dataframe = {len(verified_ids_allnull)}\n'
          f'      - {"\n  - ".join(verified_ids_allnull)}')
else:
    print(f'   verified IDs dataframe = None')

# ITS variant dataframe
if len(itsvar_counts_allnull) > 0:
    itsvar_counts_cleaned.drop(axis=1, columns=itsvar_counts_allnull, inplace=True)
    print(f'   ITS variant dataframe = {len(itsvar_counts_allnull)}\n'
          f'      - {"\n   - ".join(itsvar_counts_allnull)}')
else:
    print(f'   ITS variant dataframe = None')

########################################################################################################################


## FILTER ROWS ######################################################################################################

# NUMBER OF DOMINANT VARIANTS #

# the important columns in the itsvar_counts df is the ones related to variants; filter any rows out that don't have
#  data in at least one of these columns
itsvar_counts_varcols = ['its_variant_count','its_variant_top_copies']
temp1 = itsvar_counts_cleaned[(~itsvar_counts_cleaned[itsvar_counts_varcols[0]].isna())]
temp2 = itsvar_counts_cleaned[(~itsvar_counts_cleaned[itsvar_counts_varcols[1]].isna())]
print(f'The {itsvar_counts_varcols[0]} column has {itsvar_counts_cleaned.shape[0]-temp1.shape[0]} missing values.\n'
      f'The {itsvar_counts_varcols[1]} column has {itsvar_counts_cleaned.shape[0]-temp2.shape[0]} missing values.\n'
      )
print(f'There are {temp1["its_variant_count"].isna().sum()} instances where a '
      f'value is in the {itsvar_counts_varcols[0]} column but not the {itsvar_counts_varcols[1]} column.\n'
      f'There are {temp1["its_variant_top_copies"].isna().sum()} instances where a '
      f'value is in the {itsvar_counts_varcols[1]} column but not the {itsvar_counts_varcols[0]} column.\n'
      )

# filter by Nr dominant variants
itsvar_counts_cleaned = temp1.copy()

# ANY MISSING SAMPLE IDS #

verified_ids_cleaned = verified_ids_cleaned[~verified_ids_cleaned['sample_id'].isna()]
itsvar_counts_cleaned = itsvar_counts_cleaned[~itsvar_counts_cleaned['sample_id'].isna()]

# SAMPLE IDS AS STRINGS #
verified_ids_cleaned = verified_ids_cleaned.astype({'sample_id':str})
itsvar_counts_cleaned = itsvar_counts_cleaned.astype({'sample_id':str})

########################################################################################################################


## MERGE DATAFRAMES ####################################################################################################

# CHECK OVERLAP IN SAMPLE_IDS #

# check how many sample IDs overlap between the two dataframes before merging
merge_col = 'sample_id'
itsvar_counts_unique_sampid = set(itsvar_counts_cleaned[merge_col])
verified_ids_unique_sampid = set(verified_ids_cleaned[merge_col])
overlapping_sampid = len(itsvar_counts_unique_sampid.intersection(verified_ids_unique_sampid))
print(f'Number of unique sample IDs in each table:\n'
      f'   ITS variant count dataframe = {len(itsvar_counts_unique_sampid)}\n'
      f'   ID verification dataframe = {len(verified_ids_unique_sampid)}\n'
      f'Number of overlapping sample IDs in both dataframes:\n'
      f'   {overlapping_sampid}\n')

# if the sample IDs don't overlap at all, don't merge
if overlapping_sampid == 0:

    print(f'Not merging dataframes...\n')
    df_merged = False

    # make a copy of the wanted dataframe to keep names consistent if or if not merging
    output_df = itsvar_counts_cleaned.copy()

# if there is any overlap, then merge the two tables
else:

    print(f'Merging dataframes...\n')
    df_merged = True

    # left join the two dataframes
    output_df = itsvar_counts_cleaned.merge(
        verified_ids_cleaned,
        on='sample_id',
        how='left',
    )

    # check number of rows; as left join, should be same number of rows as itsvar_counts_cleaned
    assert output_df.shape[0] == itsvar_counts_cleaned.shape[0], 'Issue w/ merging dataframes.'

########################################################################################################################


## CLEAN DOMINANT VARIANT COUNT COLUMN #################################################################################

# FIND INVALID DOMINANT VARIANT COUNTS #

# determine which values in the its_variant_count column of the merged dataframe are non-numeric
valid_itsvar_counts = []
invalid_itsvar_counts = []

for itsvar_count in output_df['its_variant_count']:
    try:
        itsvar_count_int = int(itsvar_count)
        valid_itsvar_counts.append(itsvar_count)
    except ValueError:
        invalid_itsvar_counts.append(itsvar_count)

valid_itsvar_percent = (len(valid_itsvar_counts) / output_df.shape[0])*100
invalid_itsvar_percent = (len(invalid_itsvar_counts) / output_df.shape[0])*100

invalid_itsvar_counts_fmt = '\n     '.join(invalid_itsvar_counts)

print(f'Out of the {output_df.shape[0]} ITS variant counts in the merged dataframe,\n'
      f'   {len(valid_itsvar_counts)} ({valid_itsvar_percent:.2f}%) were valid integers\n'
      f'   {len(invalid_itsvar_counts)} ({invalid_itsvar_percent:.2f}%) were not valid integers:\n'
      f'     {invalid_itsvar_counts_fmt}')


# DROP NON-NUMERIC DOMINANT VARIANT COUNTS #

output_df_cleaned = output_df[output_df['its_variant_count'].isin(valid_itsvar_counts)].astype({'its_variant_count':int})


########################################################################################################################


## EXPORT MERGED + CLEANED DATAFRAME ###################################################################################

# use the output file path for the cleaned dataframe as the location to write out the cleaned merged dataframe
output_df_cleaned.to_csv(cleaned_df_out_path, index=False)

########################################################################################################################
