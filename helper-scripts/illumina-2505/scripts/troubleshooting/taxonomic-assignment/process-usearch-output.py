from pathlib import Path
import re, argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.style as stl
from climush import config
from climush.utilities import exit_process
from Bio import SeqIO

# set style for figures
stl.use('fivethirtyeight')

# one of the problematic OTUs that don't have a correct taxonomic assignment
PROBLEM_OTU = 'illumina_MSC0893_1101:21943:1227'




## FILE PATHS ##########################################################################################################

## FUNCTIONS ############################

# create an output file name based on the format of the input file name
def output_filepath(input_filename, output_parentdir, output_filetag, output_filefmt):

    # create a list of the input file name components
    input_filename_parts = input_filename.name.split('.')

    # create the file name string for the output file path
    output_filename = f'{input_filename_parts[0]}_{output_filetag}.{input_filename_parts[1]}.{output_filefmt}'

    # combine the file name string with the parent directory to create a full output file path
    return output_parentdir / output_filename

## INPUT ################################

## MAIN DIRECTORY ##

# main data directory for illumina-2505
illumina2505_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/')

## DATA DIRECTORIES ##

# path to the data directory for illumina-2505
illumina2505_data = illumina2505_main / 'data'

# path to any taxonomy-related data for illumina-2505
illumina2505_tax = illumina2505_data / 'illumina-2505_taxonomy'

## TROUBLESHOOTING DIRECTORIES ##

# output path for taxonomic assignment troubleshooting output
tax_troubleshoot_main = illumina2505_main / 'scripts' / 'troubleshooting' / 'taxonomic-assignment'
tax_troubleshoot_output = tax_troubleshoot_main / 'output'

## INPUT FILE(S) ##

# path to the output produced from running vsearch --usearch_global on the top12 OTUs from illumina-2505
top12_usearch_input = illumina2505_tax / 'illumina-2505_clust_top12_2025-06-26.usearch.txt'

# path to the top12 OTU sequence .fasta file
top12_seqs_input = tax_troubleshoot_output / 'illumina-2505_clust_top12_2025-06-26.fasta'


## OUTPUT ###############################

## TROUBLESHOOTING DIRECTORIES ##

# make a new subdirectory in the output directory for figures
tax_troubleshoot_output_figs = tax_troubleshoot_output / 'figures'
tax_troubleshoot_output_figs.mkdir(exist_ok=True)

## DATAFRAMES ##

# original dataframe with headers added
usearch_wheaders_out = output_filepath(
    input_filename=top12_usearch_input,
    output_parentdir=tax_troubleshoot_output,
    output_filetag='w-headers',
    output_filefmt='csv',
)

# reformatted dataframe, unfiltered
usearch_reformat_out = output_filepath(
    input_filename=top12_usearch_input,
    output_parentdir=tax_troubleshoot_output,
    output_filetag='reformatted',
    output_filefmt='csv',
)

## FIGURES ##

# distribution of reference squence hits with USEARCH per OTU
usearch_seqhit_hist = output_filepath(
    input_filename=top12_usearch_input,
    output_parentdir=tax_troubleshoot_output_figs,
    output_filetag='seqhit-histogram',
    output_filefmt='png',
)

########################################################################################################################


## IMPORT DATA #########################################################################################################

## FORMAT COLUMNS #######################

# usearch dataframe doesn't have column headers, so create some here
usearch_colnames = ['otu_id', 'taxonomy_str', 'confidence']
usearch_dtypes = [str, str, float]

# specify the datatype for each of the columns
usearch_dtypes_bycol = { name: type for name, type in zip(usearch_colnames, usearch_dtypes)}


## IMPORT DATAFRAME #####################

# import .txt file as pandas df
usearch_df_original = pd.read_table(
    top12_usearch_input,
    names=usearch_colnames,
    dtype=usearch_dtypes_bycol,
)

########################################################################################################################


## EXPORT INPUT DATA W/ ADDED HEADERS ##################################################################################

# make copy with column headers
usearch_df_original.to_csv(usearch_wheaders_out, index=False)

########################################################################################################################


## REFORMAT DATA #######################################################################################################


## MAKE DATAFRAME COPY ##################

# before making any changes to the dataframe, create a copy and make changes to copy
usearch_df = usearch_df_original.copy()


## REFORMATTING FUNCTIONS ###############

def subset_str(dataframe, column_label, regex, action, new_label=None, on_copy=False):
    '''
    Pull information from rows in a column of a dataframe using a
    regular expression. The input column may either be overwritten
    with the subsetting data or a copy of the dataframe with an
    updated column or new column will be returned.

    :param dataframe: dataframe containing the column to pull
    information from
    :param column_label: str; the name of the column in the input
    dataframe from which to pull information
    :param regex: a regular expression that generates a match to
    the subset of the string you want to pull from the larger
    string in a column's rows
    :param on_copy: True/False; subset the information on a copy of
    the input dataframe and return an updated version of the input
    dataframe without overwriting the original input version
    :return: if on_copy, returns a dataframe; if not on_copy, returns
    None and updates input dataframe in place
    '''

    ## CHECK FUNCTION INPUT #########################################

    ## ACTION ##

    # valid options for action
    valid_action = ['replace', 'append']

    # confirm that a valid option was chosen
    invalid_action_err = f'The input for action {action} is not valid option: {valid_action}'
    assert action in valid_action, invalid_action_err

    # if the action chosen is 'append', check for a new column label with new_label
    if (action == 'append') and (new_label is None):
        new_label = f'newcol_{column_label}'

    ## ON_COPY ##

    # create a copy if on_copy=True
    if on_copy:
        input_df = dataframe.copy()

    # otherwise, rename input variable, which will update input dataframe
    else:
        input_df = dataframe


    ## SUBSET COLUMN DATA ###########################################

    ## ACTION = 'REPLACE' ##

    if action == 'replace':

        # use .apply() on the input column, with a lambda function that will
        #   return matches to the provided regex; replace original input column
        input_df[column_label] = input_df[column_label].apply(
            lambda x: re.search(regex, x).group(0) if re.search(regex, x) else pd.NA,
        )

    ## ACTION = 'APPEND' ##

    elif action == 'append':

        # use .apply() on the input column, with a lambda function that will
        #   return matches to the provided regex; add new output column
        input_df[new_label] = input_df[column_label].apply(
            lambda x: re.search(regex, x).group(0) if re.search(regex, x) else pd.NA,
        )

    ## INVALID INPUT ##

    # this shouldn't be able to happen at this point but...
    else:
        return exit_process(invalid_action_err)


    ## RETURN DATAFRAME OR NONE #####################################

    if on_copy:
        return input_df
    else:
        return None


## REFORMAT ROW VALUES ##################

## OTU ID ##

# clean up the otu_id column to only include the OTU string
otu_regex= r'(?<=otu=).+?(?=;size)'
subset_str(
    dataframe=usearch_df,
    column_label='otu_id',
    regex=otu_regex,
    action='replace',
    new_label=None,
    on_copy=False,
)

## REFERENCE NUMBERS ##

# pull the first reference number, which may be a GenBank or UNITE (UDB) number
first_refnum_regex=r'^[A-Z]{1,3}\d+(?=|)'
subset_str(
    dataframe=usearch_df,
    column_label='taxonomy_str',
    regex=first_refnum_regex,
    action='append',
    new_label='refnum_01',
    on_copy=False,
)

# pull the second reference number, which is always a UNITE species hypothesis (SH)
second_refnum_regex = r'(?<=|)SH\d+\.\d{2}FU(?=;)'
subset_str(
    dataframe=usearch_df,
    column_label='taxonomy_str',
    regex=second_refnum_regex,
    action='append',
    new_label='refnum_02',
    on_copy=False,
)

## TAXONOMY LEVELS ##

# create a list of the taxonomic levels to use as new column labels
taxonomic_lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# for each taxonomic level, create a regex for pulling this info from the taxonomy_str column
#   join the taxonomic level string (new column labels) with their associated regex
taxonomic_regex = { tax_str: f'(?<={tax_str[0]}:).+?(?=,)' for tax_str in taxonomic_lvls }

# replace the regex for species with one that does not have the positive lookahead but instead terminator
taxonomic_regex.update({'species': f'(?<=s:).+?$'})

# iterate through each taxonomic level and create a column with values for each
for tax_lvl, tax_re in taxonomic_regex.items():

    subset_str(
        dataframe=usearch_df,
        column_label='taxonomy_str',
        regex=tax_re,
        action='append',
        new_label=tax_lvl,
        on_copy=False,
    )

## DROP COLUMNS #########################

# now that all relevent info is pulled from the taxonomy_str column, drop this column
usearch_df.drop(
    columns=['taxonomy_str'],
    axis=1,
    inplace=True,
)

########################################################################################################################


## EXPORT REFORMATTED TABLE ############################################################################################

# write the reformatted dataframe to a .csv file
usearch_df.to_csv(usearch_reformat_out, index=False)

########################################################################################################################


## MAKE A DICTIONARY OF OTU SEQUENCES ##################################################################################

# create an empty dictionary to add the OTU ID as the key and OTU sequence as the value
top12_otu_seqs = {}

# open the top12 OTU sequence .fasta file...
with open(top12_seqs_input, 'r') as fasta_in:

    # iterate through each sequence record in the file...
    for otu_record in SeqIO.parse(fasta_in, 'fasta'):

        # get the OTU ID, using the same regex as done with the USEARCH results table
        otu_id = re.search(otu_regex, otu_record.id).group(0)

        # get the OTU sequence
        otu_seq = str(otu_record.seq)

        # add these to dictionary in key/value pairs
        top12_otu_seqs.update({otu_id:otu_seq})

########################################################################################################################


## SUMMARIZE USEARCH HITS ##############################################################################################

## SUMMARIZE STATS FUNCTIONS ############

# create a function that will produce a dictionary of the standard summary stats
def summary_stats_dict(summary_column, dataframe, print_results):
    '''
    Produce a dictionary of standard summary statistics given a
    column of a dataframe. Summary statistics include: mean, median,
    standard deviation, minimum, and maximum.

    :param summary_column: the column of the dataframe containing the
    values to be summarized
    :param dataframe: the dataframe containing the column of values to
    summarize
    :param print_results: True/False; whether to print a summary of the
    statistics produced in addition to the output dictionary
    :return: a dictionary of summary statistics based on the values
    provided to the function
    '''

    # pull the values from the input dataframe from which to summarize data
    data_to_summarize = dataframe[summary_column]

    # create a dictionary where the key is the summary statistic type and
    #   the value is the statistic based on the data provided
    summary_stats = {
        'mean': data_to_summarize.mean(),
        'std': data_to_summarize.std(),
        'median': data_to_summarize.median(),
        'min': data_to_summarize.min(),
        'max': data_to_summarize.max(),
    }

    # print out a summary of the summary statistics calculated, if True
    if print_results:

        summary_stats_fmt = '\n   '.join(
            [ f'{stat_name} = {stat_val:.2f}' for stat_name, stat_val in summary_stats.items() ]
        )

        summary_msg = (f'The following are summary statistics on the {summary_column} of '
                       f'the input dataframe:\n'
                       f'   {summary_stats_fmt}')

        print(summary_msg)

    # otherwise, do not print anything
    else:
        pass

    # always return the summary stats dict
    return summary_stats

## HITS PER OTU #########################

# summarize the number of hits that an OTU has
otu_hit_counts = usearch_df.groupby('otu_id').count()

# what is the average number of hits per OTU?
count_proxy = 'confidence'
hits_per_otu = summary_stats_dict(
    summary_column=count_proxy,
    dataframe=otu_hit_counts,
    print_results=True,
)

# plot the distribution of the number of hits per OTU
fig, ax = plt.subplots(figsize=(6,6))

sns.histplot(
    data=otu_hit_counts,
    x=count_proxy,
    bins=10,
    ax=ax,
)

# change the x- and y-labels
ax.set_xlabel(
    'Number of Reference Sequence Hits with USEARCH',
    fontsize=14,
    labelpad=10,
)
ax.set_ylabel(
    'Number of OTUs',
    fontsize=14,
    labelpad=10,
)

# include a main title and subtitle
main_title = 'illumina-2505'
sub_title = 'Top 12 Most Abundant OTUs'

# main title
ax.text(x=0.5, y=1.1,
        s=main_title,
        fontsize=14,
        weight='bold',
        ha='center',
        va='bottom',
        transform=ax.transAxes,
        )

# subtitle
ax.text(x=0.5, y=1.05,
        s=sub_title,
        fontsize=12,
        ha='center',
        va='bottom',
        transform=ax.transAxes,
        )

fig.tight_layout()
plt.savefig(usearch_seqhit_hist)
plt.close()


## SEQ LENGTH VS OTU HIT COUNT ##########

# create a series from the sequence dictionary, where index is otu_id and values are sequence lengths
seqlength_df = pd.DataFrame(
    data={
        'otu_id': top12_otu_seqs.keys(),
        'otu_seqlength': [ len(seq) for seq in top12_otu_seqs.values() ],
        'otu_seq': top12_otu_seqs.values(),
    },
)


# add the sequence length to the summary table
otu_hit_counts = otu_hit_counts.reset_index('otu_id').merge(seqlength_df, on='otu_id')

# create two new columns for segments of the OTU sequences..

# one for the first 20 bp
otu_hit_counts['otu_seq_head'] = otu_hit_counts['otu_seq'].apply(lambda x: x[:20])

# one for the last 20 bp
otu_hit_counts['otu_seq_tail'] = otu_hit_counts['otu_seq'].apply(lambda x: x[-20:])


########################################################################################################################

# pulled from pattern recognized also in row 11,
its1catta_dete = 'GTAGGTGAACCTGCGGAAGGATCATTA'
its1catta_real =         'ACCWGCGGARGGATCATTA'

total_match_count = 0
for i in np.arange(otu_hit_counts.shape[0]):

    # search for the ITS1catta primer substring
    its1catta_hit = re.search(its1catta, otu_hit_counts['otu_seq'].iloc[i])

    # if a match is found...
    if its1catta_hit:

        # add to total match counter
        total_match_count += 1

        # get the OTU ID
        otu_id = otu_hit_counts['otu_id'].iloc[i]

        # get the location of match in sequence
        match_loc = its1catta_hit.span()

        # print the OTU ID
        print(f'ITS1catta match found: \n'
              f'   otu id      = {otu_id}\n'
              f'   match start = {match_loc[0]}\n'
              f'   match end   = {match_loc[1]}\n')


total_match_percent = (total_match_count / otu_hit_counts.shape[0]) * 100
print(f'Total ITS1catta match count: {total_match_count} / {otu_hit_counts.shape[0]} ({total_match_percent:.1f}%)')