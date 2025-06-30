from pathlib import Path
from Bio import SeqIO
import argparse
from climush.utilities import run_subprocess, exit_process, convert_udb_format

## COMMAND LINE OPTIONS ################################################################################################


## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Display the first n sequence read headers from a sequence file.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

## FILE PATHS ##

# (required) input file path
parser.add_argument(
    '-i', '--input',
    type=Path,
    required=True,
    help='File path to a sequence file from which n sequence read headers will be displayed.'
)

## DISPLAY SETTINGS ##

# (optional) the number of sequence read headers to display
parser.add_argument(
    '-n','--n-headers',
    type=int,
    required=False,
    default=10,
    help='The number of sequence read headers to return/display from the input sequence file.'
)

## PARSE ARGUMENTS ##

# parse arguments into a dictionary
args = vars(parser.parse_args())


########################################################################################################################


## CHECK INPUT FILE FORMAT + CONVERT IF NEEDED #########################################################################

# use the function convert_udb_format() to check if a conversion from .udb needs to occur, which it will execute if needed
args['input'], db_converted = convert_udb_format(
    input_file = args['input'],
)

########################################################################################################################


## PRINT INPUT SEQUENCE FILE HEADERS ###################################################################################

# use a counter to determine when enough headers have been displayed
header_count = 0

# get the format of the converted or checked input file
input_fmt = args['input'].suffix
if input_fmt == '.fa':
    input_fmt = '.fasta'
else:
    pass

# read the input sequence file in using biopython's SeqIO parser
with open(args['input'], 'r') as seqfile_in:

    # iterate through the records in the input sequence file
    for record in SeqIO.parse(seqfile_in, input_fmt.replace('.','')):

        # update the counter to keep track of the number of headers displayed
        header_count += 1

        # print the header for this reference sequence in the input sequence file
        print(f'{record.id}\n')

        # check the counter to see if the iterator should continue (n_headers not yet reached)...
        if header_count < args['n_headers']:
            continue

        # or if the number of headers returned meets the n_headers argument, break and stop printing / iterating
        else:
            break

########################################################################################################################


## REMOVED CONVERTED DB ################################################################################################

# if the input file path needed to be converted...
if db_converted:

    # remove the output of the conversion (i.e., not the original)
    args['input'].unlink()

# if no conversion happened, end here
else:
    pass

########################################################################################################################