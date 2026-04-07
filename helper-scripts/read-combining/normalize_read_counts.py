# functions used directly by this function
from climush.utilities import create_file_list


# functions for testing only
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# create dummy variables for the eventual input parameters
input_dir = Path('')        # can be a directory path or a list of files
output_dir = Path('')       # can be a file path or, if None, will replace the input files
reference_dir = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline')

## IMPORT BIOINFORMATICS CONFIGURATION SETTINGS ##

# read in the settings from the climush bioinformatics configuration file
settings = get_settings(reference_dir)
# get the run name from the configuration file
run_name = settings['run_details']['run_name']

## CREATE LIST OF INPUT FILES ##

# if input is already a list...
if isinstance(input_dir, list):
    # rename the variable to align with input_dir that's a directory
    input_files = input_dir
# if input is not a list
else:
    # create a list of sequence files in this directory
    input_files = create_file_list(input_dir)

