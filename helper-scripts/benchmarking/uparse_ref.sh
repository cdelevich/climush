#!/bin/bash

module load usearch

usearch -uparse_ref $1 -db ./MockCommunitySequences_clean.fasta -strand plus -fastaout $2.fasta --uparseout $2.up
