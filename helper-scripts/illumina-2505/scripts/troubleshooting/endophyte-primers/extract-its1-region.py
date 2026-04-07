## PACKAGE IMPORTS #####################################################################################################

# python standard library
import argparse, re
from pathlib import Path

# biopython modules + classes
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# climush package functions
from climush.bioinfo import dereplicate

########################################################################################################################


## COMMAND LINE OPTIONS ################################################################################################

## INSTANTIATE PARSER ##########################

parser = argparse.ArgumentParser(
    prog='extract-its1-region',
    description='Dereplicate sequences and extract the ITS1 region.',
    epilog='This script is part of the CliMush bioinformatics pipeline.'
)

## FILE PATHS ##################################

# input file / directory path
parser.add_argument(
    '-i', '--input',
    required=True,
    type=Path,
    help='Path to the input file or directory containing input files to process.',
)

# output file / directory path
parser.add_argument(
    '-o', '--output',
    required=True,
    type=Path,
    help='Path to the output directory.',
)

## PROCESSING OPTIONS ##########################

# the number of sequences to include from the input file
parser.add_argument(
    '-n','--n-samples',
    required=False,
    type=int,
    help='The number of samples to randomly draw from the input sequence file to run through itsx. '
         'If this option is not used, then the input sequences will not be subsetted and all sequences '
         'will be run through itsx after dereplication.',
)


## PARSE ARGS INTO DICTIONARY ##################

args = vars(parser.parse_args())

########################################################################################################################


## DEREPLICATE READS ###################################################################################################

derep_out = dereplicate(
    input_files=args['input'],
    output_dir=args['output'],
    reference_dir='',
    min_count=1,
    derep_step=1,
    keep_log=True,
)

########################################################################################################################