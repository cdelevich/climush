# import python packages
import argparse, re, pathlib
import numpy as np
from pathlib import Path
from Bio import SeqIO

# imports from climush package
from climush.constants import SEQ_FILE_GLOB, READ_COUNT_OG_RE
from climush.utilities import exit_process


## DEFINE + PARSE COMMAND LINE OPTIONS

# set up command line options
parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Select the top reads in a PacBio sequencing file.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# input directory containing the files with sequences to cluster; no default
parser.add_argument('-i', '--input',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory containing sequencing files to process.')

# output directory to write clustered sequence files to
parser.add_argument('-o', '--output',
                    default=None,
                    type=pathlib.PosixPath,
                    help='The path to a directory to write the processed sequencing files to.')

# the number of reads or percentage of total reads to keep per sample
parser.add_argument('--threshold',
                    default=None,
                    help='The cut off value for selecting the top sequences in a sample. If an integer is provided, '
                         'this will be used as the number of sequences to include in the output, ranked by the most '
                         'abundant sequences in terms of their full-length read count. If a float is provided, then '
                         'this will be used as the percentile of reads to include in the output, also ranked by the '
                         'most abundant sequences by full-length read count.')

# parse command line options and defaults into a dictionary
args = vars(parser.parse_args())

## REMOVE AFTER TESTING ########################################################################################
test_threshold = .98

args = {'input': Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pacbio_sporocarp-f_2024-06_Q40_dereplicated'),
        'output': None,
        'threshold': test_threshold}
################################################################################################################


## FORMAT THRESHOLD INPUT


# if its an integer...
if isinstance(args['threshold'], int):

    # make value even; pre-itsx, reads are combination of fwd and revcomp reads of same seq; will derep after
    if args['threshold'] % 2 != 0:
        args['threshold'] = args['threshold'] + 1
    else:
        pass

    # add a value to the args dictionary indicating threshold input type
    args.update({'threshold_type': 'top_number'})

# if its a float...
elif isinstance(args['threshold'], float):

    # add a value to the args dictionary indicating threshold input type
    args.update({'threshold_type': 'top_percentile'})

# if it is neither, this isn't a valid input
else:
    print(f'Invalid input type for filtering threshold. Value must either be an integer (number of top unique '
          f'reads to keep) or a float (percentage of top unique reads to keep).')
    # technically this would do better as a promopt and recursive function, but too lazy right now
    sys.exit()


## SET DEFAULT OUTPUT IF NONE PROVIDED

# create a default output path if one is not provided from command line
if args['output'] is None:

    # create output directory within the input directory, adding a folder name prefix of 'top-reads'
    args['output'] = args['input'] / (args['input'].name + '_top-reads')

# create the output directory
args['output'].mkdir(exist_ok=True)


## FILTER READS

# returns a numpy array of read counts from a sample
def get_read_counts(sequence_file):

    # set counters for total read count and unique read count
    sample_read_counts = []

    # get the format of the sequencing file; SeqIO does not want punctuation, so remove it
    file_fmt = sequence_file.suffix.replace('.', '')

    # ensure that the reads have been dereplicated; need read count from header
    for read in SeqIO.parse(sequence_file, file_fmt):

        # get the read header, which should contain the read count
        read_header = read.id

        # search the read header for the word 'size', which is the copy info added to the header by dereplication
        size_found = re.search(READ_COUNT_OG_RE, read_header, re.I)

        # get the number of copies of this unique read from the read header
        # if the read size value is detected in the header...
        if size_found:
            # get the read count value, cast as an integer, then add it to the total read count
            read_count = int(size_found.group(0))
            sample_read_counts.append(read_count)
        else:
            # print error and exit if the size variable isn't detected in the read header, meaning derep hasn't
            #  been carried out
            err_msg = (f'The \'size\' variable from dereplication was not detected in the read header for '
                       f'the sequencing file:\n'
                       f'   {input_fastx.name}\n'
                       f'Reads must be dereplicated before selecting the top reads for a sample.\n'
                       f'Exiting...')
            exit_process(err_msg)

    # after going through all reads in the sample, return the total and unqiue read counts
    return np.array(sample_read_counts)


# for each sequencing file in the input directory...
for input_fastx in args['input'].glob(SEQ_FILE_GLOB):

    ## REMOVE AFTER TESTING #######
    # input_fastx = next(args['input'].glob(SEQ_FILE_GLOB))
    ###############################

    # get the format of the sequencing file; SeqIO does not want punctuation, so remove it
    file_fmt = input_fastx.suffix.replace('.', '')

    # if the threshold value is a percentile, calculate the cut-off read count based on this percentile
    if args['threshold_type'] == 'top_percentile':

        # create an array of all read counts in this sample
        sample_read_counts = get_read_counts(input_fastx)

        # calculate the percentile based on the sample's read counts and the input percent value
        read_threshold = np.percentile(sample_read_counts, int(test_threshold*100))  # percent should be int

    # if the threshold value is the top x sequences, do nothing to it
    else:
        read_threshold = args['threshold']

    # counter for tracking how many reads have been cycled through
    reads_processed = 0

    # list of SeqIO records to write to output file
    top_read_records = []

    # create a variable that allows loop to go one value past the percentile threshold, if necessary, to get an even
    #   number of sequences in the output files (to increase chances of pairing fwd and revcomp seqs later)
    make_even = False

    # go through each read and determine which reads to keep in the filtered output fasta
    for read in SeqIO.parse(input_fastx, file_fmt):

        # get the read count of this sequence
        read_count = int(re.search(READ_COUNT_OG_RE, read.id, re.I).group(0))

        # check this sequence's read count against the threshold value

        # if selecting the top x (int) reads...
        if args['threshold_type'] == 'top_x_reads':

            # if the number of reads processed is still below the threshold...
            if reads_processed < args['threshold']:

                # add this sequence to the list of sequences to include in output fasta file
                top_read_records.append(read)

            # if the number of reads processed is equal to or above the threshold...
            else:

                # break here, and move on to writing output sequences to a file
                break

        # if selecting by the top % of reads...
        else:

            # compare the current read count against the percentile threshold
            if (read_count >= read_threshold) or (make_even == True):

                # add this sequence to the list of sequences to include in output fasta file
                top_read_records.append(read)

            # if the number of reads processed is equal to or above the threshold...
            else:

                # first check that there is an even number of sequences included (same orientation reason as above)
                if len(top_read_records) % 2 == 0:

                    # break here, and move on to writing output sequences to a file
                    break

                else:

                    # if there are an odd number of sequences, add one more sequence the output list to make it even
                    make_even = True


        # add to the reads processed counter
        reads_processed += 1

    # write out the top read records to a new fastx file (match the input file format)

    # print a warning if the threshold value used did not appreciably reduce the number of sequences in the sample
    # WARNING NOT COMPLETE; IS IT WORTH COUNTING TOTAL READS FOR ALL FILES TO DO THIS?

    # create the output file name based on the input file name; add top-reads suffix
    output_fastx = (args['output'] / (input_fastx.stem + '_top-reads')).with_suffix(('.' + file_fmt))

    # write out sequence file
    SeqIO.write(top_read_records, output_fastx, file_fmt)

print(f'Completed.')