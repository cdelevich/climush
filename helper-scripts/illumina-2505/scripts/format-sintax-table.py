from pathlib import Path
import re
import pandas as pd
import pandas.errors
from Bio import SeqIO
from climush.utilities import sort_taxonomy_info

# number of sequences reported to be classified by vsearch --sintax standard output
VSEARCH_CLASSIFIED_SEQS = 138919
VSEARCH_PROCESSED_SEQS = 139026


## FILE PATHS ##########################################################################################################

## INPUT #############################

# main data paths for the illumina-2505 sequence group taxonomy
illumina2505_data_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/data/')
illumina2505_data_tax = illumina2505_data_main / 'illumina-2505_taxonomy/'

# directory for the sintax-only taxonomic assignments
illumina2505_data_tax_sintax = illumina2505_data_tax / 'sintax-taxonomy_2025-07/'

# the tabularized SINTAX taxonomic assignment output dataframe
sintax_tab_in = illumina2505_data_tax_sintax / 'illumina-2505_98-clusters_nonsingletons_98-SINTAX.txt'

# the .fasta sequence file provided to the SINTAX taxonomic assignment algorithm
query_seqs_in = illumina2505_data_tax_sintax / 'illumina-2505_98-clusters_nonsingletons.fasta'

## OUTPUT ############################


########################################################################################################################


## IMPORT DATA #########################################################################################################

## TAXONOMY TABLE ####################

# sintax taxonomy tables don't include headers, so they need to be added
sintax_tax_tab_colnames = [
    'otu_id',
    'taxonomy_predicted',
    'match_orientation',
    'taxonomy_predicted_restricted',
]

def auto_handle_error_rows(df_path: Path, **kwargs):

    # if a keyword argument for error row is provided, get value, otherwise, set to None
    err_rows = kwargs.get('error_row', None)
    col_names = kwargs.get('column_names', None)

    try:

        output_table = pd.read_table(
            df_path,
            header=None,
            names=col_names,
            skiprows=err_rows,
        )

        return output_table

    # if any rows have parsing errors from pandas
    except pd.errors.ParserError as read_error:

        read_error_str = str(read_error)

        # if the error is a recognized ParserERror...
        if read_error_str.startswith('Error tokenizing data'):

            # get the index of the error row from the error message
            #  MUST SUBTRACT ONE BECAUSE DF DOESNT HAVE HEADER AND PANDAS WILL REPORT ROW VALUE AS ONE INDEX GREATER
            new_err_row = int(re.search(r'(?<=line )\d+(?=,)', read_error_str).group(0)) - 1

            # if no error rows have previously been detected...
            if err_rows is None:

                # create a list and add this value to it
                err_rows = [new_err_row]

            # if there is already a list of error rows...
            else:
                # append this new error row to the list
                err_rows.append(new_err_row)

            # recursively run this function providing the list of known error rows to skip over
            return auto_handle_error_rows(
                df_path=df_path,
                error_row=err_rows,
                column_names=col_names,
            )

        else:
            err_msg = (f'Unrecognized error when reading {sintax_tab_in.name}:\n '
                       f'{read_error}\n')
            print(err_msg)
            return None

sintax_df = auto_handle_error_rows(
    df_path=sintax_tab_in,
    column_names=sintax_tax_tab_colnames,
)

print(f'')


## NON-SINGLETON SEQS ################
with open(query_seqs_in, 'r') as seqs_in:
    seq_generator = SeqIO.parse(seqs_in, 'fasta')

########################################################################################################################


## REFORMAT OTU_ID #####################################################################################################

# clean up the otu_id column to only include the OTU string
# sintax_df['otu_id'] = sintax_df['otu_id'].apply(lambda x: re.search(r'(?<=otu=).+?(?=;size)', x).group(0))

otu_str_clean = []
for otu_val in sintax_df['otu_id']:
    otu_str_found = re.search(r'(?<=otu=).+?(?=;size)', otu_val)
    if otu_str_found:
        otu_str_clean.append(otu_str_found.group(0))
    else:
        otu_str_clean.append(otu_val)

sintax_df['otu_id'] = otu_str_clean

########################################################################################################################

## CHECK OTU ID OF ERROR OTU FROM AMPTK TAXONOMY #######################################################################

PROBLEM_OTU = 'illumina_MSC0893_1101:21943:1227'
sintax_df[sintax_df['otu_id'] == PROBLEM_OTU]['taxonomy_predicted_restricted'].iloc[0]

########################################################################################################################


## REFORMAT PREDICTED TAXONOMY W/O BOOTSTRAP CONFIDENCE ################################################################

sintax_df_reformat = sort_taxonomy_info(
    input_df=sintax_df,
    tax_column='taxonomy_predicted_restricted',
    drop_col=False,
)



## REFORMAT MATCH ORIENTATION ##########################################################################################

# reformat the reference_strand column so that symbols (+/-) are strings
sintax_df['reference_strand'] = sintax_df['reference_strand'].apply(lambda x: 'forward' if x == '+' else 'reverse')

PROBLEM_OTU = 'illumina_MSC0893_1101:21943:1227'
