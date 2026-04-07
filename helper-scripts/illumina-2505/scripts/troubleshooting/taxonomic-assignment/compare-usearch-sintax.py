from pathlib import Path
import pandas as pd
import numpy as np
import argparse, re
from climush import config
from climush.utilities import sort_taxonomy_info

PROBLEM_OTU = 'illumina_MSC0893_1101:21943:1227'


# run interactively or via terminal?
try:
    print(f'Running {Path(__file__).name}...\n')
    interactive_mode = False
except NameError:
    print(f'Running interactively...\n')
    interactive_mode = True

    # path to the default input file; must be str because parser will later convert to Path object
    amptk_table_in = '/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/data/illumina-2505_taxonomy/illumina-2505_98-clusters.taxonomy.txt'


## COMMAND LINE OPTIONS ################################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(
    prog='usearch-sintax-compare',
    description='Compare the taxonomy assigned by USEARCH and SINTAX taxonomic assignment algorithms '
                'as run by the amptk taxonomy.',
    epilog='This script is part of CliMush bioinformatics.',
)

## FILE PATHS ##

# REQUIRED; input file
parser.add_argument(
    '-i', '--input',
    required=True,
    type=Path,
    help='The path to the input taxonomy table produced by amptk taxonomy, which typically has the file name format '
         '<basename>.taxonomy.txt.',
)

## PARSE ARGUMENTS ##

# parse command line arguments into a dictionary

# if running interactively (within Python console in PyCharm)...
if interactive_mode:

    # feed default command line options directly to parser
    args = vars(
        parser.parse_args(
            ['--input', amptk_table_in]
        )
    )

# if running script from a Terminal...
else:
    # parse arguments from the command line input
    args = vars(parser.parse_args())

########################################################################################################################


## IMPORT TAXONOMY TABLE ###############################################################################################

# if running interactively in PyCharm...
if interactive_mode:

    # if the table is already loaded into the environment...
    if 'tax_results' in globals():

        # prompt whether you want to re-load it (time consuming)
        check_prompt = f'Are you sure you want to re-import the input taxonomy dataframe [True/False]?'
        reload_input = input(check_prompt)

        # if True, reload
        if reload_input == 'True':
            print(f'Reloading table...\r')
            tax_results = pd.read_table(args['input'])
        else:
            pass

    # if the table has not already been imported, then import
    else:
        tax_results = pd.read_table(args['input'])

# if running from Terminal, read table in without prompt
else:
    tax_results = pd.read_table(args['input'])

########################################################################################################################


## CLEAN DATAFRAME #####################################################################################################

# create a copy of the original input dataframe
tax_df = tax_results.copy()

# what columns are in the dataframe?
print(f'the following columns are in the input dataframe, {args["input"].name}')
print('\n'.join(tax_results.columns))

# print first few rows
tax_results.head()

## PULL CLUSTER/OTU INFO ##

# create a regex compiler for each of the info categories to pull from the #OTUID column
otu_regex_search = re.compile(r'(?<=otu=).+?(?=;size)')         # OTU ID
size_regex_search = re.compile(r'(?<=size=)\d+(?=;clusterid)')  # the size (read count) of the OTU
clustnum_regex_search = re.compile(r'(?<=clusterid=)\d+$')      # the cluster number assigned to the OTU by vsearch

otu_ids = []
read_sizes = []
cluster_nums = []

for otu_str in tax_df['#OTUID']:

    otu_ids.append(otu_regex_search.search(otu_str).group(0))
    read_sizes.append(size_regex_search.search(otu_str).group(0))
    cluster_nums.append(clustnum_regex_search.search(otu_str).group(0))


## CREATE NEW OTU INFO COLUMNS ##

# insert the OTU string as the first column in the dataframe
tax_df.insert(
    loc=0,
    column='otu_id',
    value=pd.Series(otu_ids),
    allow_duplicates=False,
)

# next add the cluster number
tax_df.insert(
    loc=1,
    column='cluster_number',
    value=pd.Series(cluster_nums),
    allow_duplicates=False,
)

# then add the read count of the cluster
tax_df.insert(
    loc=2,
    column='read_count',
    value=pd.Series(read_sizes),
)


## DROP COLUMNS ##

# now that info has been extracted from #OTUID, drop this column
tax_df.drop(
    '#OTUID',
    axis=1,
    inplace=True,
)

# this version of amptk doesn't run UTAX, so all will say No hit
tax_df.drop(
    'UTAX',
    axis=1,
    inplace=True,
)


## SORT TAXONOMY STRING ##

# create a list of each taxonomic level; will append algorithm-specific prefix for each species id method
taxonomic_lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# create a function to pull taxonomic information from the taxonomy column into organized columns
def pull_taxon_info(input_df, tax_info_order_dict, tax_str_sep, info_type, results_taxon_sep=';', input_df_tax_col='taxonomy'):

    for new_col, str_loc in tax_info_order_dict.items():

        accepted_info_types = {
            'amptk_method': 0,
            'taxonomic_levels': 1,
        }

        if info_type in accepted_info_types.keys():
            info_loc = accepted_info_types[info_type]
        else:
            err_msg = f'Input {info_type} not one of valid info_type values: {list(accepted_info_types.keys())}'
            raise ValueError(err_msg)

        # gather items into a list to later assign to new column name
        new_col_values = []

        # go through each taxonomy string in the dataframe column for taxonomy
        for tax_str in input_df[input_df_tax_col]:

            # split the full taxonomy string first into tax results info and tax taxonomy info
            # then split the taxonomy tax results part into its components, splitting by its separator
            tax_str_list = tax_str.split(results_taxon_sep)[info_loc].split(tax_str_sep)

            # check if this info is available for this particular row / OTU
            if len(tax_str_list) > str_loc:

                # if info_type is taxonomic_levels, strip the taxon level prefix before adding as new row value
                new_col_values.append(tax_str_list[str_loc].split(':')[-1])

            # if it isn't, append a pd.NA to list
            else:
                new_col_values.append(pd.NA)

        # create new column from the list of new column values
        tax_df[new_col] = new_col_values

# create a dictionary for the location of items related to taxonomic assignment results returned by amptk
tax_results_cols = ['amptk_method', 'amptk_confidence', 'amptk_genbank', 'amptk_unite']
tax_results_order = { result_col: order for result_col, order in zip(tax_results_cols, np.arange(0,len(tax_results_cols))) }
tax_results_sep = '|'

pull_taxon_info(
    input_df=tax_df,
    tax_info_order_dict=tax_results_order,
    tax_str_sep=tax_results_sep,
    info_type='amptk_method',
)

# create a dictionary for the location of items related to taxonomic levels (kingdom through species)
tax_taxon_cols = [ 'amptk_' + lvl for lvl in taxonomic_lvls ]
tax_taxon_order = { taxon_col: order for taxon_col, order in zip(tax_taxon_cols, np.arange(0,len(tax_taxon_cols))) }
tax_taxon_sep = ','

# pull taxonomic info for each taxonomic level to create new columns
pull_taxon_info(
    input_df=tax_df,
    tax_info_order_dict=tax_taxon_order,
    tax_str_sep=tax_taxon_sep,
    info_type='taxonomic_levels',
)

########################################################################################################################


sintax_tax = sort_taxonomy_info(
    input_df=tax_df,
    tax_column='SINTAX',
    drop_col=True,
)
