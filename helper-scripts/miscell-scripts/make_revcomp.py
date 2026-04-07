'''
THIS DOESNT WORK, UNICODEENCODE (OR DECODE?) ERROR WHEN READING FASTA WITH SEQIO PARSE
'''

from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

mock_in = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/miscell-scripts/MockCommunitySequences_clean.fasta')
mock_out = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/miscell-scripts/MockCommunitySequences_clean_w-revcomps.fasta')

err_records = []
revcomp_fwd_records = []
with open(mock_in, 'r') as seqs_in:
    for record in SeqIO.parse(seqs_in, 'fasta'):
            revcomp_seq = record.seq.reverse_complement()
            revcomp_id = record.id + '_revcomp'
            revcomp_record = SeqRecord(
                revcomp_seq,
                id=revcomp_id,
                name='',
                description='',
            )
            revcomp_fwd_records.append(record)
            revcomp_fwd_records.append(revcomp_record)



with open(mock_out, 'wt') as seqs_out:
    SeqIO.write(revcomp_fwd_records, seqs_out, 'fasta')
