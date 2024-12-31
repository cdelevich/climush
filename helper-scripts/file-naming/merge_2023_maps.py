from pathlib import Path
import pandas as pd

## SET FILE PATHS ######################################################################################################

## INPUT ##

# path to the main directory containing the CTAB mapping files to concatenate
mapping_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/renaming-mapping-files/from-neda-arnold-lab/')

# path to the oregon-only mapping file
or_path = mapping_main / 'Arad_Macrosystems_MSC600-MSC719_Carolyn_23Dec2024.xlsx'

# path to the original mapping file that excludes oregon samples
og_path = mapping_main / 'its1_endophyte-mapping_2023.csv'

## OUTPUT ##

# path to the concatenated output file
concat_path = mapping_main / 'its1_endophyte-mapping_2023_combined.csv'

########################################################################################################################


## IMPORT DATAFRAMES ###################################################################################################

## OR TABLE ##

# or samples are in Excel (.xlsx) format; take only relevant columns
or_df = pd.read_excel(or_path)[['CTAB Codes MSC####', 'BAG LABEL']]

# must have same column name as the original table to later concatenate
or_df.rename(columns={'CTAB Codes MSC####': 'Code'}, inplace=True)

## ORIGINAL TABLE ##

# other samples in .csv format; take only relevant columns
og_df = pd.read_csv(og_path)[['Code', 'BAG LABEL']]

# confirm column names are the same
assert list(or_df.columns) == list(og_df.columns), ('Column names of the two tables are not the same and therefore will cause '
                                                    'issues when trying to concatenate them.')

########################################################################################################################


## CONCATENATE TABLES ##################################################################################################

# concatenate the two dataframes into one table
combined_df = pd.concat([og_df, or_df])

# calculate the number of unique CTAB codes in each of the source and combined dataframes
original_ctab_count = len(og_df['Code'].unique())
oregon_ctab_count = len(or_df['Code'].unique())
combined_ctab_count = len(combined_df['Code'].unique())
print(f'unique CTAB code counts:\n'
      f'   combined mapping file = {combined_ctab_count}\n'
      f'   original mapping file = {original_ctab_count}\n'
      f'   OR mapping file       = {oregon_ctab_count}\n')

# because the numbers aren't what is expected, also calculate the number of unique climush sample IDs in each dataframe
original_sid_count = len(og_df['BAG LABEL'].unique())
oregon_sid_count = len(or_df['BAG LABEL'].unique())
combined_sid_count = len(combined_df['BAG LABEL'].unique())
print(f'unique climush sample ID counts:\n'
      f'   combined mapping file = {combined_sid_count}\n'
      f'   original mapping file = {original_sid_count}\n'
      f'   OR mapping file       = {oregon_sid_count}\n')

assert (original_ctab_count + oregon_ctab_count) == combined_ctab_count, ('Error during concatenation - number of unique '
                                                                          'CTAB codes in the two source files does not '
                                                                          'equal the number in the combined file.')

# which climush sample ID is duplicated in the original mapping file?
dup_climush_sample = og_df[og_df['BAG LABEL'].duplicated()]['BAG LABEL'].iloc[0]
print(f'The duplicated climush sample ID in the original mapping file is:\n'
      f'   {dup_climush_sample}')

dup_ctab_codes = og_df[og_df['BAG LABEL'] == dup_climush_sample]['Code'].to_list()
print(f'The CTAB codes associated with the duplicated climush sample ID in the original mapping file are:\n'
      f'   {dup_ctab_codes}')

########################################################################################################################


## EXPORT COMBINED DATAFRAME ###########################################################################################

# write out concatenated dataframe, regardless of errors/issues, to output file path
combined_df.to_csv(concat_path, index=False)

########################################################################################################################
