#!/bin/bash

alias sort_files='python3 pipeline/00_sort-files.py'
alias demux='python3 pipeline/01_demultiplex.py'
alias prefilter='python3 pipeline/02_prefilter.py'
alias primer_trim='python3 pipeline/03_remove-primers.py'
alias qual_filter='python3 pipeline/04_quality-filtering.py'
alias derep01='python3 pipeline/05_dereplicate-01.py'
alias itsx='python3 pipeline/06_separate-subregions.py'
alias chimera='python3 pipeline/07_chimera-check.py'
alias derep02='python3 pipeline/08_dereplicate-02.py'
alias cluster='python3 pipeline/09_otu-clustering.py'
alias taxonomy='python3 pipeline/10_assign-taxonomy.py'
alias otu_table='python3 pipeline/11_create-otu-table.py'
