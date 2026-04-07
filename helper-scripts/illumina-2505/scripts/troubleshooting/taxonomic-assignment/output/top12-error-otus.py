from pathlib import Path
import re
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 20)


## FILE PATHS ##########################################################################################################

## INPUT ##

# main directory path for the combined illumina sequence run, illumina-2505
illumina2505_main = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505')

# directory to the data files pulled from CliMush Sequences (Globus) for illumina-2505
illumina2505_data = illumina2505_main / 'data'

# path to the main troubleshooting directory within the illumina-2505 directory
illumina2505_troubleshoot = illumina2505_main / 'scripts' / 'troubleshooting'

# path to the troubleshooting directory for taxonomic assignment
illumina2505_troubleshoot_tax = illumina2505_troubleshoot / 'taxonomic-assignment'

# path to the output directory for this troubleshooting task
illumina2505_troubleshoot_tax_output = illumina2505_troubleshoot_tax / 'output'

# OTU table clean
otu_tab_path = illumina2505_troubleshoot_tax_output / 'illumina-2505_top12_2025-06-26.csv'


## OUTPUT ##

# date file suffix to include in output files
output_date_suffix = datetime.now().strftime('%Y-%m-%d')


########################################################################################################################

## READ IN FILE

otu_tab = pd.read_csv(otu_tab_path)

problem_otus = ['illumina_MSC0893_1101:21943:1227', 'illumina_MSC0953_1101:43340:31783']

problem_otu_tab = otu_tab[otu_tab['otu_id'].isin(problem_otus)]

problem_otu_its1 = problem_otu_tab['its1_sequence'].to_list()