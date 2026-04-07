from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from climush import config
from climush.constants import SAMPLE_COL_RE
import re

# path to the xlsx file that marcos sent me with the included collection numbers
coll_list_path = Path('/Users/carolyndelevich/main/projects/fungi-of-ordway/pacbio-samples/Pacbio_accession_numbers.xlsx')

# location of this script and 2024-04 file list
otu_sort_path = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/otu-sorting')
file_list_path = otu_sort_path / 'itsx_pacbio_sporocarp-f_2024-04_Q40_output-files.txt'

# import df, take first tab as the df
coll_df = pd.read_excel(coll_list_path, sheet_name=0)

# the column ORD seems to have the ORD### sample ID, along with what I think is the read number (e.g., ORD58_853169)
# create a dict with the sample ID as the key, and the read ID (with sample ID included) as the value
coll_dict = {read_id.split('_')[0]: read_id for read_id in coll_df['ORD']}

# downloaded a list of the files in the ITSx output for pacbio_sporocarp-f_2024-04_Q40 from Globus collection
with open(file_list_path, 'rt') as fin:

    # collect the name of the unique ORD sample IDs located among the pacbio_sporocarp-f_2024-04 itsx output
    ord_samples = set()

    # go through each file name in the list of itsx output files
    for filename in fin.readlines():

        # search for the ORD sample ID in the file name
        ord_sample_found = re.search(r'(?<=_)ORD\d{1,4}(?=\.)', filename, re.I)

        # if an ORD sample ID is located in the file name
        if ord_sample_found:

            ord_samples.add(ord_sample_found.group(0))

        else:

            continue


print(f'{len(ord_samples)} ORD samples located in the ITSx output on Globus for pacbio_sporocarp-f_2024-04_Q40.')


# see how many and which samples from 2024-04 are in the list of samples included in the manuscript
manuscript_samples = set(coll_dict.keys())
in_manuscript_globus = ord_samples.union(manuscript_samples)

if len(in_manuscript_globus) == len(ord_samples):
    print(f'All {len(in_manuscript_globus)} ORD samples in the manuscript list are in the 2024-04 collection.')
else:
    print(f'Only {len(in_manuscript_globus)} ORD samples are in both the manuscript list and the ITSx output from '
          f'2024-04 on Globus.')
    in_manuscript_only = manuscript_samples.difference(ord_samples)
    in_globus_only = ord_samples.difference(manuscript_samples)
    print(f'   only in manuscript         = {len(in_manuscript_globus)} \n'
          f'   only in 2024-04 collection = {len(ord_samples)}\n')


# this is confusing and I don't want to have to be searching like this every time going forward
# better to just create a table that has the sample ID and what sequencing run that sample was sequenced with
bc_mapping_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/config/barcode-mapping/')

# create an empty dictionary; the key will be the pacbio sequencing run name and the value is a list of the samples
#  that were included on that run
sample_mapping_dict = {}

for bc_map in bc_mapping_path.glob('pacbio*barcode-mapping*'):

    # get the pacbio sequencing run name belonging to this set of sequences
    sequencing_run = re.search(r'^pacbio_sporocarp-f_\d{4}-\d{2}(?=_barcode-mapping\.)', bc_map.name).group(0)

    # add the name of the sequencing run to the dict with an empty list as its value
    sample_mapping_dict.update({sequencing_run: []})

    # read in the mapping file to get the sample IDs included in the sequencing run
    if bc_map.suffix == '.xlsx':

        # most will be excel files with multiple file tabs; pools split across two tabs
        for tab_name, df in pd.read_excel(bc_map, sheet_name=None).items():

            if re.search('pool', tab_name, re.I):

                sample_id_col = 
