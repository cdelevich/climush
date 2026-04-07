from pathlib import Path
import re, argparse
import pandas as pd
from climush import config

PROBLEM_OTU = 'illumina_MSC0893_1101:21943:1227'


# sintax output file path
sintax_output_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/data/illumina-2505_taxonomy/illumina-2505_clust_top12_2025-06-26.sintax.txt')

# sintax output .txt file doesn't include headers, so create column names here
sintax_colnames = [
    'otu_id',
    'predicted_taxonomy',
    'reference_strand',
    'predicted_taxonomy_simple',
]

# read in sintax .txt output file as a pandas df
sintax_df = pd.read_table(
    sintax_output_path,
    names=sintax_colnames,
)

# reformat the reference_strand column so that symbols (+/-) are strings
sintax_df['reference_strand'] = sintax_df['reference_strand'].apply(lambda x: 'forward' if x == '+' else 'reverse')

# clean up the otu_id column to only include the OTU string
sintax_df['otu_id'] = sintax_df['otu_id'].apply(lambda x: re.search(r'(?<=otu=).+?(?=;size)', x).group(0))

# split the taxonomy in the simplified taxonomy column
# create a list of each taxonomic level
taxonomic_lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# show the species ID from SINTAX for the problem OTU

sintax_df[sintax_df['otu_id'] == PROBLEM_OTU]