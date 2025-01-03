from pathlib import Path
import pandas as pd
import re
from climush.constants import SEQ_FILE_GLOB, GZIP_REGEX, SEQ_FILE_RE, CORRECT_CTAB_PREFIX, ANY_CTAB_PREFIX, ANY_CTAB_CODE

## SET FILE PATHS ######################################################################################################

## INPUT ##

# path to the main directory containing the CTAB mapping files
mapping_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/config/renaming-mapping-files/from-neda-arnold-lab/')

# path to the concatenated mapping file, which should contain all known/expected CTAB codes
concat_path = mapping_main / 'its1_endophyte-mapping_2023_combined.csv'

# path to the main directory for file naming
file_naming_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/')

# path to the directory containing the blank sequence files, originally indended for testing file renaming
endo_seqfile_dir = file_naming_main / 'test-files_2023-sequence-files/'

## OUTPUT ##

# path to the file containing the CTAB codes and climush sample IDs in the mapping file but missing in sequence files
# knew I only needed the one output path because I came back after investigating to add path here
missing_from_seqfile_out = file_naming_main / 'its1_endophtyes_2023_missing-from-seqfiles.csv'

########################################################################################################################

## GET CTAB LIST FROM SEQUENCE FILE NAMES ##############################################################################

# create an empty set to add CTAB/slant codes to
# use a set (not list) because R1/R2 files will otherwise result in duplicate CTAB/slant codes
seqfile_ctab_codes = set()

# keep track of number of files, might be useful later?
total_seq_files = 0
ctab_seq_files = 0

# go through each of the sequence files that match CTAB-type codes
for seq_file in endo_seqfile_dir.glob(SEQ_FILE_GLOB):

    # add to the seqfile counter
    total_seq_files += 1

    # check the file for a match to a CTAB-type file prefix (MSC/MSCSC/etc.)
    ctab_match = re.search(ANY_CTAB_CODE, seq_file.name, re.I)

    # if this file does match a CTAB-type file prefix...
    if ctab_match:

        # add to the CTAB match counter
        ctab_seq_files += 1

        # add the CTAB code to the CTAB/slant code set
        seqfile_ctab_codes.add(ctab_match.group(0))

    # if this file doesn't appear to match a CTAB code, pass over it
    else:
        continue


########################################################################################################################


## GET CTAB LIST FROM MAPPING FILE #####################################################################################

# read in the CTAB mapping file created by combining the original and oregon mapping files
mapping_df = pd.read_csv(concat_path)

# replace the slant codes (MSCS) with CTAB codes (MSC) in the Code column
# sequence files are in the CTAB code (MSC) format, so this is required in order to properly compare
mapping_df['ctab_code'] = [re.sub(ANY_CTAB_PREFIX, CORRECT_CTAB_PREFIX, code) for code in mapping_df['Code']]

# create a set of unique CTAB codes from the 'Code' column of the dataframe
# want a set, rather than creating list from .unique(), so that I can use set-type comparisons
mapping_ctab_codes = set(mapping_df['ctab_code'].to_list())

########################################################################################################################


## COMPARE MAPPING FILE AND SEQUENCE FILE CTAB CODES ###################################################################

## COUNTS ##

# compare just the number of unique CTAB codes from the sequence files and mapping file
seqfile_unique_count = len(seqfile_ctab_codes)
mapping_unique_count = len(mapping_ctab_codes)
print(f'number of unique CTAB codes:\n'
      f'   sequence files = {seqfile_unique_count}\n'
      f'   mapping file   = {mapping_unique_count}\n')


## VALUES ##

# compare the actual values of these CTAB codes from the mapping and sequence files

# in the mapping file, but missing from the sequence files
missing_from_seqfile = mapping_ctab_codes.difference(seqfile_ctab_codes)
print(f'missing from the sequence files:\n'
      f'   number of missing CTAB codes = {len(missing_from_seqfile)}\n'
      f'   missing CTAB codes = {missing_from_seqfile}\n')

# in the sequence files, but missing from the mapping file
missing_from_mapping = seqfile_ctab_codes.difference(mapping_ctab_codes)
print(f'missing from the mapping file:\n'
      f'   number of missing CTAB codes = {len(missing_from_mapping)}\n'
      f'   missing CTAB codes = {missing_from_mapping}\n')

########################################################################################################################


## SUMMARIZE MISSING INFORMATION FOR OUTPUT ############################################################################

# CTAB codes were only missing from the sequence files, none from the mapping file

# create a small dataframe of the missing CTAB codes and associated climush sample IDs
missing_from_seqfile_df = mapping_df[mapping_df['ctab_code'].isin(missing_from_seqfile)]

# rename the Code and climush ID columns for clarity; rearrange column
missing_df_out = missing_from_seqfile_df.rename(columns={'Code': 'slant_code',
                                                         'BAG LABEL': 'old_climush_id',})[['ctab_code', 'slant_code', 'old_climush_id']]


# write out this summary dataframe to a .csv file
missing_df_out.to_csv(missing_from_seqfile_out, index=False)

########################################################################################################################