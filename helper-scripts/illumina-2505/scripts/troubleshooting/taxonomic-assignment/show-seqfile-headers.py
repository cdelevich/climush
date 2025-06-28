from pathlib import Path
from Bio import SeqIO
import argparse
from climush.utilities import run_subprocess, exit_process

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

# assign a boolean flag for whether a database (.udb) conversion took place
# this is so I can delete the converted copy after
db_converted = False

## DETERMINE INPUT FORMAT ##

# if the input format has multiple suffixes, this means it is likely a compressed .fastq.gz file
if len(args['input'].suffixes) > 1:

    # not expecting this right now, so I will not handle gzip
    err_msg = f'Input file is compressed and I have not written a solution to this.'
    exit_process(err_msg)

# if the input format has a single suffix, then continue on to check the format type
else:
    input_fmt = args['input'].suffix


## CHECK INPUT FORMAT ##

# check if the input format is a simple sequence file format
seqfile_fmts = ['.fastq', '.fasta', '.fa']

# if the input file format is in the list of expected file formats...
if input_fmt in seqfile_fmts:

    # do nothing else here but assign the sequence file to a variable
    input_seqfile = args['input']

# otherwise, check if something can be converted to be an expected file format...
else:

    # if the input format is .udb format...
    if input_fmt == '.udb':

        # first check if an .extracted.fa version of this database exists
        if args['input'].with_suffix('.extracted.fa').is_file():

            # if this exists, then use this version of the database as the input, no need to convert separately
            input_fmt = '.fasta'
            input_seqfile = args['input'].with_suffix('.extracted.fa')

        # if an .extracted.fa version of this database does not exist, convert to .fasta
        else:

            # switch the Boolean flag to True
            db_converted = True

            # assign a new input file format
            input_fmt = '.fasta'

            # create an output .fasta format file name
            fasta_out = args['input'].with_suffix(input_fmt)

            # compile the list of commands to run to convert udb to .fasta format
            udb_convert_cmd = [
                'vsearch',        # call vsearch
                '--udb2fasta',    # use its udb2fasta conversion function
                args['input'],    # input sequence file in the .udb file format
                '--output',       # designate the output file
                fasta_out,        # name of the output .fasta file format
            ]

            # use the run_subprocess() function from utilities.py to execute
            run_subprocess(
                cli_command_list=udb_convert_cmd,
                dest_dir=args['input'].parent,
                run_name='',
                program='vsearch-ubd2fasta',
                separate_sample_output=True,
            )

            # rename the .fasta path variable to input_seqfile, which is input for next section
            input_seqfile = fasta_out

    # if not one of the expected input sequence file formats or .udb to convert, print err and exit
    else:
        err_msg = f'Unrecognized / invalid input sequence file format: {input_fmt}'
        exit_process(err_msg)


########################################################################################################################


## PRINT INPUT SEQUENCE FILE HEADERS ###################################################################################

# use a counter to determine when enough headers have been displayed
header_count = 0

# read the input sequence file in using biopython's SeqIO parser
with open(input_seqfile, 'r') as seqfile_in:

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
    input_seqfile.unlink()

# if no conversion happened, end here
else:
    pass

########################################################################################################################