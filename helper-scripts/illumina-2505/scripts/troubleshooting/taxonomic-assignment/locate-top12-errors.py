from pathlib import Path
import re
import pandas as pd
import numpy as np
from datetime import datetime
from Bio import SeqIO
from climush.constants import READ_COUNT_OG_RE
from climush.utilities import sort_taxonomy_info

OTU_ID_RE = r'(?<=;otu=).+?(?=;size)'
OTU_COLUMN_LABEL = '#OTU ID'
TAX_STR_RE = r'(?<=\s).+'

suppress_colerr = True

pd.set_option('display.max_columns', 10)
pd.set_option('display.max_rows', 10)

## FILE PATHS ##########################################################################################################

## INPUT ##

# main directory path for the combined illumina sequence run, illumina-2505
illumina2505_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505')

# directory to the data files pulled from CliMush Sequences (Globus) for illumina-2505
illumina2505_data = illumina2505_main / 'data'

# path to clustered sequences from illumina-2505
illumina2505_data_clustered = illumina2505_data / 'illumina-2505_clustered_2025-05-22'

# path to the OTU table for illumina-2505
illumina2505_otutab_path = illumina2505_data_clustered / 'illumina-2505_98-clusters_otu-table.txt'

# path to the OTU .fasta sequence file for illumina-2505
illumina2505_otufast_path = illumina2505_data_clustered / 'illumina-2505_98-clusters.fasta'

# path to the taxonomic assignment output files produced by amptk taxonomy for illumina-2505
illumina2505_data_taxonomy = illumina2505_data / 'illumina-2505_taxonomy'

# path to the OTU .fasta sequence file WITH TAXONOMY for illumina-2505
illumina2505_otufast_wtax_path = illumina2505_data_taxonomy / 'illumina-2505_98-clusters.otus.taxonomy.fa'

## OUTPUT ##

# path to the main troubleshooting directory within the illumina-2505 directory
illumina2505_troubleshoot = illumina2505_main / 'scripts' / 'troubleshooting'

# path to the troubleshooting directory for taxonomic assignment
illumina2505_troubleshoot_tax = illumina2505_troubleshoot / 'taxonomic-assignment'

# create an output directory for this troubleshooting task
illumina2505_troubleshoot_tax_pathout = illumina2505_troubleshoot_tax / 'output'
illumina2505_troubleshoot_tax_pathout.mkdir(exist_ok=True)

# date file suffix to include in output files
output_date_suffix = datetime.now().strftime('%Y-%m-%d')

# top n OTUs basename (n value appended when writing out below)
topn_otus_clust_basename = 'illumina-2505_clust_top'
topn_otus_tax_basename = 'illumina-2505_tax_top'

########################################################################################################################


## IMPORT DATA #########################################################################################################

## COMMENTED OUT CODE BECAUSE THE TABLE ESPECIALLY TAKES LONG TO LOAD ##################################################

## OTU TABLE ##

# import the OTU table as a pd dataframe
# illumina2505_otutab = pd.read_table(illumina2505_otutab_path, delimiter='\t')


## OTU FASTA ##

# I DON'T KNOW HOW TO GO FROM A SEQRECORD DICT TO WRITING SO I DID NOT CONTINUE WITH THIS METHOD

# import the OTU .fasta file as a dictionary
# with open(illumina2505_otufast_path, 'r') as fast_in:
#     illumina2505_otufast = SeqIO.to_dict(SeqIO.parse(fast_in, format='fasta'))

########################################################################################################################


## FILTER FASTA TO TOP N OTUS ##########################################################################################

# create a function so this can be done for both the clustered and taxonomy-assigned .fasta files
def pull_topn_seqrecords(top_n, full_seqpath, output_dir, output_basename, return_seqrecords=False, regex={'read_count':READ_COUNT_OG_RE}):


    ## PULL TOP N SEQUENCE RECORDS ##

    # gather the sequence records of the top n OTUs into a list
    topn_otu_seqrecords = []

    # iterate through each record (OTU) in the .fasta file and pull top n
    with open(full_seqpath, 'r') as fast_in:
        otu_record_count = 0  # counter to compare against number of wanted records, top_n
        for otu_record in SeqIO.parse(fast_in, format='fasta'):

            # if top n have not yet been selected...
            if otu_record_count < top_n:

                # add this record to the list of top n OTU sequence records
                topn_otu_seqrecords.append(otu_record)

                # if this is the last record, then break here
                if otu_record_count == (top_n-1):
                    break
                # otherwise, add to the counter and continue
                else:
                    otu_record_count += 1
                    continue

    if len(topn_otu_seqrecords) == top_n:
        print(f'SUCCESS. {top_n} OTU sequence records were successfully pulled from the full OTU .fasta file,'
              f'{full_seqpath.name}.\n')
    else:
        err_msg = (f'{len(topn_otu_seqrecords)} OTU sequence records were pulled from the full OTU .fasta file,'
                   f'{full_seqpath.name}, when {top_n} should have been pulled.')
        raise KeyboardInterrupt(err_msg)


    ## CONFIRM DESCENDING READ COUNT ORDER ##

    # of these n OTUs, confirm they are sorted by decreasing size (read count)
    previous_read_size = 1e10  # use very large number, will compare to last to ensure it is larger
    assessment_counter = 0     # to confirm all OTUs were iterated through without error
    for otu_record in topn_otu_seqrecords:

        # add one to the assessment counter
        assessment_counter += 1

        # get the size value from the read header
        size_found = re.search(regex['read_count'], otu_record.id)
        if size_found:
            read_size = int(size_found.group(0))
        else:
            print(f'size not found for:\n'
                  f'   {otu_record.id}')
            break

        # compare current read_size to previous
        # current read size should be less than the previous one, if in descending read count order
        if read_size < previous_read_size:
            # update the previous_read_size to this current read size
            previous_read_size = read_size
        else:
            # print message that the read size isn't decreasing
            err_msg = (f'The current read count, {read_size}, is not less than the previous OTU\'s read count, '
                       f'{previous_read_size}. This indicates that the OTU .fasta file, {full_seqpath.name}, '
                       f'is not sorted by decreasing read count.\n')
            print(err_msg)
            break

    # check that all the OTUs were worked through
    if assessment_counter == top_n:
        print(f'SUCCESS. OTU .fasta file, {full_seqpath.name}, was sorted by decreasing read count abundance.\n')
    else:
        print(f'ERROR.')


    ## EXPORT TOP N OTU RECORDS AS .FASTA OUTPUT FILE ##

    # create the output file date suffix
    output_date_suffix = datetime.now().strftime('%Y-%m-%d')

    # write out the top n OTUs to a separate .fasta file
    topn_fasta_filename = f'{output_basename}{top_n}_{output_date_suffix}.fasta'
    topn_fasta_pathout = output_dir / topn_fasta_filename
    with open(topn_fasta_pathout, 'w') as fasta_out:
        SeqIO.write(topn_otu_seqrecords, fasta_out, 'fasta')

    # confirm file exists
    if topn_fasta_pathout.is_file():
        print(f'SUCCESS. The top {top_n} OTU .fasta output file, {topn_fasta_filename}, was created.\n'
              f'   {topn_fasta_pathout}\n')
    # if not, raise error
    else:
        err_msg = f'There was an error writing the top {top_n} OTU .fasta output file, {topn_fasta_filename}, to a file.'
        raise OSError(err_msg)

    # return the list of sequence records if wanted
    if return_seqrecords:
        return topn_otu_seqrecords
    # otherwise return None
    else:
        return None

# for the clustered .fasta file
pull_topn_seqrecords(
    top_n=12,
    full_seqpath=illumina2505_otufast_path,
    output_dir=illumina2505_troubleshoot_tax_pathout,
    output_basename=topn_otus_clust_basename,
    return_seqrecords=False,
)

# for the taxonomy .fasta file
topn_tax_seqrecords = pull_topn_seqrecords(
    top_n=12,
    full_seqpath=illumina2505_otufast_wtax_path,
    output_dir=illumina2505_troubleshoot_tax_pathout,
    output_basename=topn_otus_tax_basename,
    return_seqrecords=True,
)

########################################################################################################################


## FILTER OTU TABLE BY TOP N OTUS ######################################################################################

# get the OTU IDs from the top n OTU sequence records
topn_otu_ids = [ re.search(OTU_ID_RE, otu_record.id).group(0) for otu_record in topn_otu_seqrecords ]

# filter the OTU ID column of the OTU table to only include these OTUs
topn_otu_taxtab = illumina2505_otutab[illumina2505_otutab[OTU_COLUMN_LABEL].isin(topn_otu_ids)].reset_index(drop=True)

# confirm the correct number of OTUs are included in this subsetted dataframe
if topn_otu_taxtab.shape[0] == top_n:
    print(f'SUCCESS. The top {top_n} OTUs were successfully pulled from the OTU table.\n')
else:

    # get the OTU IDs that are included in this subsetted dataframe
    otus_pulled = topn_otu_taxtab[OTU_COLUMN_LABEL]

    # limit the amount of OTU IDs that print in case of error resulting in large number
    ALLOWED_PRINT_AMT = 15

    # format the OTU IDs that were included for print-out
    if len(otus_pulled) <= ALLOWED_PRINT_AMT:
        otus_pulled_fmt = '\n   '.join(otus_pulled)
    else:
        # only print a limited number in case of error (so doesn't print entire table)
        otus_allowed_toprint = otus_pulled[:ALLOWED_PRINT_AMT]
        otus_allowed_toprint.append(f'... (+ {len(otus_pulled) - ALLOWED_PRINT_AMT} more OTUs)')
        otus_pulled_fmt = '\n   '.join(otus_allowed_toprint)

    # print error message and interrupt script
    err_msg = (f'ERROR. {topn_otu_taxtab.shape[0]} OTUs were pulled from the OTU table instead of '
               f'the expected {top_n}:\n'
               f'   {otus_pulled_fmt}')

    raise KeyboardInterrupt(err_msg)

## MOVED BELOW BECAUSE I NEED TO ADD THE TAXONOMIC INFO TO THE OTU TABLE BEFORE WRITING TO FILE ########################
# ## WRITE OUT THE FILTERED OTU TABLE ##
#
# # use basename to create output OTU table file name
# topn_otutab_filename = f'{topn_otus_basename}{top_n}_{output_date_suffix}.csv'
# topn_otutab_pathout = illumina2505_troubleshoot_tax_pathout / topn_otutab_filename
# topn_otu_taxtab.to_csv(topn_otutab_pathout, index=False)
#
# # confirm file exists
# if topn_otutab_pathout.is_file():
#     print(f'SUCCESS. The top {top_n} OTU .csv output file, {topn_otutab_filename}, was created.\n'
#           f'   {topn_fasta_pathout}\n')
# # if not, raise error
# else:
#     err_msg = f'There was an error writing the top {top_n} OTU .csv output file, {topn_otutab_filename}, to a file.'
#     raise OSError(err_msg)

########################################################################################################################

## SUMMARIZE DATA TO DROP COLUMNS ######################################################################################

# create a copy of the original OTU table after filtering for top n OTUs
topn_otutab = topn_otu_taxtab.copy()

# rename the OTU column
topn_otutab.rename(
    columns={
        '#OTU ID': 'otu_id',
    },
    inplace=True,
)

# function to get the column names in the df matching the provided regex
def get_column_names(regex, df):
    return [ col for col in df.columns if re.search(regex, col, re.I) ]

def get_group_read_counts(regex, df):

    # get the name of the columns for this group
    matching_colnames = get_column_names(
        regex=regex,
        df=df
    )

    # calculate the read counts for each OTU in the input dataframe for this grouping of samples
    return [df[matching_colnames].iloc[i].sum() for i in np.arange(df.shape[0])]

column_regex = [
    r'(mock)|(_PT)',
    r'spore_\d{4}',
    r'soil_\d{4}',
    r'(?<=_)litter_\d{4}',
    r'(endo)|(MSC)|(root-)|(leaf-)',
]

# check that the regex I'm using captures all of the sample columns (all but otu_id column)
sample_cols_fnct = [ ]
for r in column_regex:
    colnames = get_column_names(r, df=topn_otutab)
    sample_cols_fnct = sample_cols_fnct + colnames

# check that all are unique
if len(sample_cols_fnct) == len(set(sample_cols_fnct)):
    pass
else:

    # figure out which are not unique (matching multiple groups)
    sample_cols_set = set()
    sample_cols_duplicates = []
    for col in sample_cols_fnct:
        if col in sample_cols_set:
            sample_cols_duplicates.append(col)
        else:
            sample_cols_set.add(col)

    # format non-unique in list
    sample_cols_duplicates.sort()
    sample_cols_duplicates_fmt = '\n   '.join(sample_cols_duplicates)

    raise KeyboardInterrupt(f'Some columns were captured by multiple groupings:\n   {sample_cols_duplicates_fmt}\n')

# check that all are included
if len(sample_cols_fnct) == (topn_otutab.shape[1] - 1):
    pass
else:

    # find which columns are missing, except for the one expected to be excluded, otu_id
    missing_cols = set(topn_otutab.columns).difference(set(sample_cols_fnct))
    missing_cols.remove('otu_id')
    missing_cols_fmt = '\n   '.join(list(missing_cols))

    # calculate percent missing
    percent_missing = (len(sample_cols_fnct)/(topn_otutab.shape[1] - 1))*100

    # compose error message
    err_msg = f'{len(sample_cols_fnct)} columns out of {topn_otutab.shape[1] - 1} '\
              f'({percent_missing:.1f}%) were captured by the grouping regex. The following '\
              f'were excluded:\n'\
              f'   {missing_cols_fmt}'

    if suppress_colerr:
        print(err_msg)
    else:
        raise KeyboardInterrupt(err_msg)


sample_group_read_counts = {
    'otu_id':topn_otutab['otu_id'].to_list(),
    'reads_mock':get_group_read_counts(r'(mock)|(_PT)',df=topn_otutab),
    'reads_spores':get_group_read_counts(r'spore_\d{4}',df=topn_otutab),
    'reads_soil':get_group_read_counts(r'soil',df=topn_otutab),
    'reads_litter':get_group_read_counts(r'(?<=_)litter_\d{4}',df=topn_otutab),
    'reads_endophytes':get_group_read_counts(r'(endo)|(MSC)',df=topn_otutab),
}

# create dataframe from dict
sample_group_read_counts_df = pd.DataFrame.from_dict(
    data=sample_group_read_counts,
    orient='columns',
)




## ADD TAX INFO TO FILTERED OTU TABLE ##################################################################################


# create a dictionary of new columns to add to the OTU table
topn_otutab_newcolumns = {
    'otu_id':[],            # need this to join correctly with existing OTU table
    'its1_sequence':[],     # ITS1 sequence of the OTU
    'total_read_count':[],        # read count of the OTU
    'taxonomy':[],          # taxonomy string assigned by amptk taxonomy
}

# go through the topn OTU records that have the taxonomy included and pull information on the OTUs
for otu_record in topn_tax_seqrecords:

    # OTU ID
    otu_read_id = re.search(OTU_ID_RE, otu_record.id).group(0)
    topn_otutab_newcolumns['otu_id'].append(otu_read_id)

    # ITS1 sequence
    topn_otutab_newcolumns['its1_sequence'].append(str(otu_record.seq))

    # read count
    otu_read_count = int(re.search(READ_COUNT_OG_RE, otu_record.id).group(0))
    topn_otutab_newcolumns['total_read_count'].append(otu_read_count)

    # taxonomy string
    otu_tax_str = re.search(TAX_STR_RE, otu_record.description).group(0)
    topn_otutab_newcolumns['taxonomy'].append(otu_tax_str)

# make sure each of the items in the new column dictionary have a value for every OTU
# then add to the table
for colname, values in topn_otutab_newcolumns.items():

    if len(values) == top_n:
        print(f'SUCCESS. {top_n} values for the {colname} new column were pulled from the OTU taxonomy .fasta file.\n')
    else:
        err_msg = (f'{len(values)} for the {colname} new column were pulled from the OTU taxonomy '
                   f'.fasta file when {top_n} were expected.\n')
        raise KeyboardInterrupt(err_msg)

# create a table of the new columns, to be joined with the original OTU table
newcol_table = pd.DataFrame.from_dict(
    data=topn_otutab_newcolumns,
    orient='columns',
)

########################################################################################################################

# join the grouped read counts and the new columns from the tax .fasta together
top_otutab_clean = newcol_table.merge(sample_group_read_counts_df, how='inner', on='otu_id')

# confirm number of OTUs is correct
if top_otutab_clean.shape[0]  == top_n:
    print(f'SUCCESS. Sample grouping read counts table and taxonomy information joined successfully.\n')
else:
    err_msg = 'error merging dataframes'
    raise KeyboardInterrupt(err_msg)


top_otutab_clean_tax = sort_taxonomy_info(
    input_df=top_otutab_clean,
    tax_column='taxonomy',
    drop_col=True,
)


## WRITE OUT THE FILTERED OTU TABLE ####################################################################################

# use basename to create output OTU table file name
topn_otutab_filename = f'{topn_otus_basename}{top_n}_{output_date_suffix}.csv'
topn_otutab_pathout = illumina2505_troubleshoot_tax_pathout / topn_otutab_filename
top_otutab_clean_tax.to_csv(topn_otutab_pathout, index=False)

# confirm file exists
if topn_otutab_pathout.is_file():
    print(f'SUCCESS. The top {top_n} OTU .csv output file, {topn_otutab_filename}, was created.\n'
          f'   {topn_fasta_pathout}\n')
# if not, raise error
else:
    err_msg = f'There was an error writing the top {top_n} OTU .csv output file, {topn_otutab_filename}, to a file.'
    raise OSError(err_msg)

########################################################################################################################