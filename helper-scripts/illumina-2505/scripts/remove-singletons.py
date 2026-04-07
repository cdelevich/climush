from pathlib import Path
import pandas as pd
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re, zlib

main_data_dir = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/illumina-2505/data/')

# illumina-2505 98% clusters
clust98_main = main_data_dir / 'illumina-2505_clustered_2025-05-22'
otu_tab_path = clust98_main / 'illumina-2505_98-clusters_otu-table.txt'
otu_seq_path = clust98_main / 'illumina-2505_98-clusters.fasta'

# per-sequence-group .fasta files
seq_group_zip = main_data_dir / 'illumina_2505_per-sequence-group.zip'

otu_tab = pd.read_table(otu_tab_path, sep='\t')

# calculate the total read count for an OTU across all samples
otu_tab['otu_read_count'] = otu_tab[1:].sum(axis=1, numeric_only=True)

# calculate total number of OTUs including singletons
otu_count = otu_tab.shape[0]

# calculate the number and percent of singletons in the OTU table
n_reads=1
otu_tab_singletons = otu_tab[otu_tab['otu_read_count'] == n_reads]
num_singletons = otu_tab_singletons.shape[0]
perc_singletons = np.round((num_singletons / otu_count) * 100,2)
singleton_otu_id = otu_tab_singletons['#OTU ID'].to_list()

# calculate the number and percent of OTUs with less than 5 total reads
n_reads=5
num_lowabund5 = otu_tab[otu_tab['otu_read_count'] < n_reads].shape[0]
perc_lowabund5 = np.round((num_lowabund5 / otu_count) * 100,2)

# print summary
print(f'total number of 98% OTUs = {otu_count}\n'
      f'   percent singletons =               {perc_singletons}% ({num_singletons} OTUs)\n'
      f'   percent low abundance (<{n_reads} reads) = {perc_lowabund5}% ({num_lowabund5} OTUs)')

# create a table where singletons are removed
otu_tab_nosingletons = otu_tab[otu_tab['otu_read_count'] > 1]
otu_count_nosingletons = otu_tab_nosingletons.shape[0]
percent_remain_nosingletons = np.round((otu_count_nosingletons/otu_count)*100,2)

# print summary
print(f'After removing singletons, the OTU table has {otu_tab_nosingletons.shape[0]} OTUs '
      f'({percent_remain_nosingletons}% of original).')

# remove singletons from centroid .fasta file based on OTU ID
singleton_re = r';size=1;'
singleton_otu_records = []
nonsingleton_otu_records = []
with open(otu_seq_path, 'rt') as seqs_in:
      for record in SeqIO.parse(seqs_in, 'fasta'):
            if re.search(singleton_re, record.description):
                  singleton_otu_records.append(record)
            else:
                  nonsingleton_otu_records.append(record)

# confirm sorting was done correctly
assert len(nonsingleton_otu_records) == otu_count_nosingletons
assert len(singleton_otu_records) == num_singletons

singleton_fasta = clust98_main / 'singletons.fasta'
with open(singleton_fasta, 'w') as single_out:
      SeqIO.write(singleton_otu_records, single_out, format='fasta')

nonsingleton_fasta = clust98_main / 'nonsingletons.fasta'
with open(nonsingleton_fasta, 'w') as nonsingle_out:
      SeqIO.write(nonsingleton_otu_records, nonsingle_out, format='fasta')



## COUNT FINAL READS PER SEQ GROUP ##