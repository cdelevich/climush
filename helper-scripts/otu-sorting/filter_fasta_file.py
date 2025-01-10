import argparse, pathlib, re
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from climush.constants import SEQ_FILE_GLOB

TEST_MODE = True

if TEST_MODE:

    print(f'WARNING. Running in test mode, so command line options are not parsed.')

    args = {
        'input': Path('/Users/carolyndelevich/main/projects/climush/bioinfo-output/pacbio_sporocarp-f_2023-03/08_itsx/from-dropbox_2023-11-27/'),
        'output': None,
        'ftag': 'simplified',
        # 'regex': '(?<=2023a_)CS\d{1,3}(?=_)',  # match Heather's CS samples only
        'regex': '(?<=pacbio_2023a_).+?(?=_otu\d{1,5})',  # match all climush sample IDs
        'simplify': True,
    }

else:

    # set command line options
    parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                     description='Filter reads in a .fasta file, with options to simplify the read headers.',
                                     epilog='This script is part of the climush bioinformatics pipeline.')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=pathlib.PosixPath,
                        help='The path to the .fasta file to filter. If the filter is to be applied to multiple .fasta '
                             'files in a directory, a directory can be passed here and all .fasta files will be filtered.')

    # add a mutually exclusive group, so that you can either provide an output path or a suffix to use in the default output path
    output_group = parser.add_mutually_exclusive_group()

    output_group.add_argument('-o', '--output',
                              required=False,
                              type=pathlib.PosixPath,
                              help='The path to write output file to. If no path is provided, '
                                   'the file will be written out to the same location as the input file.')

    output_group.add_argument('-f', '--ftag',
                              required=False,
                              help='A file tag to use in the default output file name.')

    parser.add_argument('-r', '--regex',
                        required=True,
                        help='A regex-formatted search string to use to filter the .fasta file.')

    parser.add_argument('--simplify',
                        action='store_true',
                        help='Boolean flag; if you want to simply the read headers from their original form to only '
                             'the matching group of the provided regex, then use this flag. Otherwise, the original '
                             'sequence file headers will be preserved.')

    # store arguments in a dictionary
    args = vars(parser.parse_args())

# check if the input path is a single file or a directory; make a list of files either way
if args['input'].is_dir():
    file_list = list(args['input'].glob(SEQ_FILE_GLOB))
else:
    file_list = list(args['input'])

# unless an output path is provided, create an output directory to write the filtered sequence files to
if args['output'] is None:
    output_parent = file_list[0].parent / f'filtered-seqs_'

# create a dictionary to add new sequence file names and the filtered SeqRecords
filtered_seqs_dict = {}

# go through each of the sequence files or file
for seq_file in file_list:

    print(f'\nFiltering {seq_file.name}...')

    # create an output sequence file name based on the input file name

    # get the sequence region; drop the period, which will still be on the region
    seq_region = seq_file.suffixes[0].replace('.','')  # itsx output format will appear to pathlib like two file suffixes

    # get the input file basename; this will remove both the sequence region and the file extension
    input_file_basename = seq_file.stem.split('.')[0]

    # if neither an output path nor an output file tag is provided, create a generic output directory and name
    if (args['output'] is None) and (args['ftag'] is None):
        output_parent = seq_file.parent / f'filtered-seqs'
        output_parent.mkdir(exist_ok=True)

        # if you use .with_suffix(), it will overwrite any part of the file name that it thinks is a suffix already
        #   (i.e., the sequence region); so don't use .with_suffix() here
        seq_file_out = output_parent / f'{input_file_basename}_filtered.{seq_region}.fasta'

    # if an output path is provided, use this directly as the output path without changes
    elif args['ftag'] is None:
        seq_file_out = args['output']

    # otherwise, a file tag was provided, so use this as the file tag instead of the generic '_filtered' tag
    else:
        output_parent = seq_file.parent / f'filtered-seqs_{args["ftag"]}'
        output_parent.mkdir(exist_ok=True)
        seq_file_out = output_parent / f'{input_file_basename}_{args["ftag"]}.{seq_region}.fasta'

    # append this output file path to the sequence dictionary as a primary key, with an empty list as the value
    #   so that the seq records that pass the filter can be added to it
    filtered_seqs_dict.update({seq_file_out: []})

    # keep track of how many total reads are in this file and how many pass the filter
    total_read_count = 0
    filter_read_count = 0

    # go through each of the sequence records in the .fasta file
    with open(seq_file, 'rt') as fasta_in:
        for read in SeqIO.parse(fasta_in, 'fasta'):

            # add to the total read count
            total_read_count += 1

            # use the regex to decide whether this read passes the filter
            read_header_match = re.search(args['regex'], read.id, re.I)

            # if this read header matches the provided regex...
            if read_header_match:

                # add to the filtered read count
                filter_read_count += 1

                # check whether the read header should be simplified before adding this record to the dict
                if args['simplify']:

                    # update the read ID for this read with the regex match group
                    read.id = read_header_match.group(0)

                    # clear any of the other items that might be appended to the header by SeqIO
                    read.description = ''
                    read.name = ''

                # if the original read header is to be preserved...
                else:
                    # don't update the read header/id for this read
                    pass

                # append this read to the sequence record dictionary for this output file
                filtered_seqs_dict[seq_file_out].append(read)

            # if this isn't a matching read ID, move to the next read in the file
            else:
                continue

    # before continuing to the next sequence file, print a summary of the number of total reads checked in this file
    #   along with the number of read that pass the provided read ID filter
    print(f'   total reads    = {total_read_count}\n'
          f'   filtered reads = {filter_read_count}\n')


# go through each output file and the associated sequence records, writing out the reads passing the
#    filter to the output file
for output_path, record_list in filtered_seqs_dict.items():

    # if there aren't any sequences that passed the filter, don't write out an empty file
    if len(record_list) == 0:
        pass
    else:
        with open(output_path, 'wt') as fasta_out:
            SeqIO.write(record_list, fasta_out, 'fasta')