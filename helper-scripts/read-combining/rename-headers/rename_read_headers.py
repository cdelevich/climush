from climush.utilities import rename_read_header, get_sample_id, mkdir_exist_ok, add_prefix, get_settings
from pathlib import Path
from Bio import SeqRecord, SeqIO
from Bio.SeqRecord import SeqRecord
import re

# path to the test files
test_file_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/read-combining/rename-headers/test-files')

# create a dictionary of file paths, with keys being the seq collection dir and the value the list of files in that dir
file_path_dict = {dir:[] for dir in test_file_path.glob('*/')}
for seq_dir in file_path_dict:
    file_path_dict[seq_dir] = [ seq_file for seq_file in seq_dir.glob('*.fast*') ]

# pull single test file
# test_dir = list(file_path_dict.keys())[0]
# test_file = file_path_dict[test_dir][0]

# see if the rename_read_header() function from utilities.py works
# rename_read_header(
#     input_dir=test_dir,
#     run_name=test_dir.name,
#     file_format='.fasta',
#     unique_headers=True,
#     no_copy=False,
#     append_sample_str=True,
# )

# get files for one sequencing group
test_dir = list(file_path_dict.keys())[0]
test_file_list = file_path_dict[test_dir]


# use platform parameter to know what information will be in the read header
platform='illumina'
output_dir=Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/read-combining/rename-headers/test-files/test-output/')
# create a reference directory path for file-finding functions
ref_dir=Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline')

# import the settings for the bioinformatics configuration
settings = get_settings(ref_dir)
run_name = settings['run_details']['run_name']

# create a list to add new sequence records to (need new seq record to change headers)
combined_seq_records = []

## CREATE OUTPUT .FASTA FILE FOR COMBINED SEQUENCES ##

mkdir_exist_ok(output_dir)
combined_seqs_out = (output_dir / f'{run_name}_combined-reads').with_suffix('.fasta')

## GO THROUGH EACH SAMPLE'S SEQUENCE FILE ##

for seq_file_in in test_file_list:

    ## GET THIS SAMPLE'S ID ##

    sample_id = get_sample_id(file_path=seq_file_in).replace('_R1','')
    sample_header = f'sample={sample_id}'

    ## GO THROUGH EACH OF THIS SAMPLE'S READS ##

    with open(seq_file_in, 'rt') as seq_in:
        for record in SeqIO.parse(seq_in, format='fasta'):

            ## GET THIS READ'S COPY NUMBER ##

            # search for the read copy number in this read's header
            try:
                size_header = re.search(r'size=\d{1,}', record.description).group(0)
            except AttributeError:
                size_header = 'size=NA'

            ## CREATE THIS READ'S UNIQUE IDENTIFIER ##

            # create a unique read identifier using the sample Id and the sequencer's header identifier
            if platform=='illumina':

                # try to get a unique combination of values from the sequencer's read identifier
                try:
                    read_id = re.search(r'(?<=:)\d{1,}:\d{1,}:\d{1,}(?=;)', record.description).group(0)
                except AttributeError:
                    read_id = ''

                # tag the read ID to the end of this sample's ID and assign to the otu= identifier in the header
                read_id_header = f'otu={sample_id}_{read_id}'

            else:
                print(f'idk yet')

            ## COMBINE HEADER ELEMENTS INTO SINGLE HEADER ##

            # combine the sample ID, otu/read ID, and size into the new read header
            updated_read_header = ';'.join([sample_header, read_id_header, size_header])

            ## CREATE A NEW SEQUENCE RECORD WITH THE NEW READ HEADER ##

            # create a new SeqRecord with this updated read header
            new_seq_record = SeqRecord(
                id=updated_read_header,
                name=updated_read_header,
                description=updated_read_header,
                seq=record.seq
            )

            ## ADD UPDATED SEQUENCE RECORD TO LIST OF UPDATED RECORDS FOR THIS SAMPLE ##

            # add this sequence record with the updated read header to the list of records to update
            combined_seq_records.append(new_seq_record)

## WRITE ALL RENAMED READS TO SINGLE FILE FOR ALL SAMPLES IN THIS SEQUENCE GROUP ##

# once new sequence records have been created for all samples for this DNA region, write out .fasta
SeqIO.write(combined_seq_records, combined_seqs_out, 'fasta')




## CREATE FUNCTION FROM THIS ##
## copied into bioinfo.py on 2025-05-21 at 9:18AM
def combine_reads(input_dir, output_dir, platform, reference_dir):

    # create list of files from input dir if not already list of files
    if isinstance(input_dir, list):
        input_files=input_dir
    else:
        input_files=create_file_list(file_input=input_dir)

    # import the settings for the bioinformatics configuration
    settings = get_settings(reference_dir)
    run_name = settings['run_details']['run_name']

    # create a list to add new sequence records to (need new seq record to change headers)
    combined_seq_records = []

    ## CREATE OUTPUT .FASTA FILE FOR COMBINED SEQUENCES ##

    mkdir_exist_ok(output_dir)
    combined_seqs_out = (output_dir / f'{run_name}_combined-reads').with_suffix('.fasta')

    ## GO THROUGH EACH SAMPLE'S SEQUENCE FILE ##

    for seq_file_in in input_files:

        ## GET THIS SAMPLE'S ID ##

        sample_id = get_sample_id(file_path=seq_file_in).replace('_R1', '')
        sample_header = f'sample={sample_id}'

        ## GO THROUGH EACH OF THIS SAMPLE'S READS ##

        with open(seq_file_in, 'rt') as seq_in:
            for record in SeqIO.parse(seq_in, format='fasta'):

                ## GET THIS READ'S COPY NUMBER ##

                # search for the read copy number in this read's header
                try:
                    size_header = re.search(r'size=\d{1,}', record.description).group(0)
                except AttributeError:
                    size_header = 'size=NA'

                ## CREATE THIS READ'S UNIQUE IDENTIFIER ##

                # create a unique read identifier using the sample Id and the sequencer's header identifier
                if platform == 'illumina':

                    # try to get a unique combination of values from the sequencer's read identifier
                    try:
                        read_id = re.search(r'(?<=:)\d{1,}:\d{1,}:\d{1,}(?=;)', record.description).group(0)
                    except AttributeError:
                        read_id = ''

                    # tag the read ID to the end of this sample's ID and assign to the otu= identifier in the header
                    read_id_header = f'otu={sample_id}_{read_id}'

                else:
                    error_msg = f'This function is not yet formatted to work with any platform other than illumina'
                    return exit_process(err_msg)

                ## COMBINE HEADER ELEMENTS INTO SINGLE HEADER ##

                # combine the sample ID, otu/read ID, and size into the new read header
                updated_read_header = ';'.join([sample_header, read_id_header, size_header])

                ## CREATE A NEW SEQUENCE RECORD WITH THE NEW READ HEADER ##

                # create a new SeqRecord with this updated read header
                new_seq_record = SeqRecord(
                    id=updated_read_header,
                    name=updated_read_header,
                    description=updated_read_header,
                    seq=record.seq
                )

                ## ADD UPDATED SEQUENCE RECORD TO LIST OF UPDATED RECORDS FOR THIS SAMPLE ##

                # add this sequence record with the updated read header to the list of records to update
                combined_seq_records.append(new_seq_record)

    ## WRITE ALL RENAMED READS TO SINGLE FILE FOR ALL SAMPLES IN THIS SEQUENCE GROUP ##

    # once new sequence records have been created for all samples for this DNA region, write out .fasta
    SeqIO.write(combined_seq_records, combined_seqs_out, 'fasta')

    # confirm file was created
    if combined_seqs_out.is_file:

        print(f'{len(combined_seq_records)} sequences were combined across {len(input_files)} samples '
              f'in the {run_name} bioinformatics group.\n')


    # return the output path
    return combined_seqs_out
