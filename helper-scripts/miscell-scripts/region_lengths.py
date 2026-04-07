from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import re
from climush.constants import POST_ITSX_SUFFIXES, SEQ_FILE_GLOB

csnm_path = Path('/Users/carolyndelevich/main/projects/climush/bioinfo-output/csnm_post-itsx/')

region = 'fullseq'
# sample_subset = r'(?<=_)CS\d{1,4}(?=_)'
sample_subset = None
for fasta in csnm_path.glob(SEQ_FILE_GLOB):

    if re.search(region, fasta.name, re.I):

        region_lens = {}

        with open(fasta, 'rt') as fasta_in:

            for read in SeqIO.parse(fasta_in, 'fasta'):

                if sample_subset is None:
                    region_lens.update({read.id: len(read.seq)})
                else:
                    csnm_found = re.search(sample_subset, read.description, re.I)

                    if csnm_found:
                        region_lens.update({csnm_found.group(0):len(read.seq)})
                    else:
                        continue

    else:
        continue

print(f'The mean length of the {region} for CSNM samples is:\n'
      f'   {np.mean(list(region_lens.values())):.2f} +/- {np.std(list(region_lens.values())):.2f}')