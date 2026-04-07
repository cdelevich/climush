from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re, gzip

# create a function that will compare how many reads in a sample have identical first n bases
def identical_first_n(seq_file, n_bases):

    # create an empty list to add the first n bases to
    seq_list = []

    # attempt to open the file as if it is not compressed
    try:
        with open(seq_file, 'rt') as sfile_in:
            for record in SeqIO.parse(sfile_in, 'fastq'):
                seq_list.append(str(record.seq)[:n_bases])

    # if the file cannot be opened due to decoding error, its likely gzipped
    except UnicodeDecodeError:
        with gzip.open(seq_file, 'rt') as sfile_in:
            for record in SeqIO.parse(sfile_in, 'fastq'):
                seq_list.append(str(record.seq)[:n_bases])

    # calculate percent unique
    percent_unique_n_bases = (len(set(seq_list)) / len(seq_list))*100

    # print the summary of identical versus unique n bp across reads in sample
    print(f'for sample file: {seq_file.name}...\n'
          f'   total reads in sample: {len(seq_list)}\n'
          f'   reads w/ unique first {n_bases} bases: {len(set(seq_list))} ({percent_unique_n_bases:.2f}%)\n')

    return seq_list

## DEMUX IDENTICAL FIRST N BASES ##

# path to the demultiplexed files
demux_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/primer-detection/test-files/illumina_endophytes_2022/01_demultiplexed')

for demux_seq_file in demux_path.glob('*.fast*'):
    identical_first_n(seq_file=demux_seq_file,n_bases=16)


## CORRECTLY TRIMMED PRIMER IDENTICAL FIRST N BASES ##

# path to samples right prior to trimming
pretrim_right_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/primer-detection/test-files/illumina_soil-litter_2022-05/02_prefiltered/02B_no-ambig')

# path to sample reads from another sequencing run that had primers trimmed from them
trim_right_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/primer-detection/test-files/illumina_soil-litter_2022-05/03_primer-trimmed')

print(f'correctly trimmed:')
print(f'   pre-trimming\n')
for untrim_seq_file in pretrim_right_path.glob('*.fast*'):
    identical_first_n(seq_file=untrim_seq_file, n_bases=16)
print(f'   post-trimming\n')
for trim_seq_file in trimmed_right_path.glob('*.fast*'):
    identical_first_n(seq_file=trim_seq_file,n_bases=16)

