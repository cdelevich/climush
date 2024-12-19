import re, json

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd
import biom

from climush.utilities import *
from climush.constants import STRICT_CLIMUSH_SAMPLE_ID_RE, LAX_CLIMUSH_SAMPLE_ID_RE

## SET FILE PATHS ######################################################################################################

# INPUT ######################

# most recent test of creating an OTU (a) with correct sample naming, (b) as a product of denoising
main_test_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/otu-tab_test/')

# size-based clustering output after proper sample labeling
clust_test_path = main_test_path / 'cluster-size'
clust_tab_path = clust_test_path / 'otu-table_renamed_test.txt'  # otu table
clust_tax_path = clust_test_path / 'otus_renamed_test.fasta'  # sequences with taxonomic assignments
clust_biom_path = clust_test_path / 'otus_renamed_test.biom'  # otu table in biom format - low memory

# denoise OTU table output after proper sample labeling
denoi_test_path = main_test_path / 'denoise'
denoi_tab_path = denoi_test_path / 'otu-table_renamed_denoise.txt'  # otu table
denoi_tax_path = denoi_test_path / 'denoise_renamed_test.fasta'  # sequences with ASV ID but *NO* taxonomic assignments
denoi_biom_path = denoi_test_path / 'otus_renamed_denoise.biom'  # otu table in biom format - low memory

# ASVs with taxonomy assigned by amptk; taxonomy not re-assigned after sample re-labeling
main_tax_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/otu-table_2024-09-28')
tax_seq_path = main_tax_path / 'otus_illumina_isl-2407.fasta.otus.taxonomy.fa'

# OTU table created by first filtering the taxonomy .fasta for certain tax group, then running size-based clustering
tax_subset_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/otu-tab_test/cluster-size/from-tax-subset')
tax_sub_otutab_path = tax_subset_path / 'tax-filt_otu-table.txt'

# OUTPUT #####################

# main directory to write output files to; make directory within the main input directory
main_output = main_test_path / 'test_summary'
main_output.mkdir(exist_ok=True)

# filtered ASV sequence file, only has ASVs that match a specified taxonomic group
# defined when used, since the taxonomic level is included in the file name
filtered_tax_prefix = 'filtered-tax'

# path for the amptk taxonomy assignment method key
tax_method_path = main_output / 'amptk-taxonomy-methods.json'


## DEFINE CONSTANTS ####################################################################################################

# delimiters
TAX_STR_DELIM = '|'
OTU_HEADER_DELIM = ';'
TAX_OTU_DELIM = ' '

# sample ID regex, by project
heather_sample_re = r'_Horseshoe_'
cams_sample_re = r'_J\d{1,3}_'
alaska_sample_re = r'_ak2022_.'

# amptk taxonomy assignment method codes
AMPTK_TAX_METHODS = ['GS', 'GSL', 'GDL', 'US', 'SS']

# build regex that will look for any of the possible amptk taxonomy assignment method codes (join w/ | 'or' operator)
AMPTK_TAX_METHOD_RE = '(' + '|'.join(AMPTK_TAX_METHODS) + ')'

# create a dictionary for the taxonomy method codes that describes what they mean
method_dict = {
    'GS': 'global alignment, bayesian and global alignment taxonomy is the same',
    'GSL': 'global alignment, bayesian and global alignment taxonomy is the same, last common ancestor algorithm used',
    'GDL': 'global alignment, bayesian and global alignment taxonomy discrepency, last common ancestor algorithm used',
    'US': 'UTAX classifier, bayesian and global alignment taxonomy is the same',
    'SS': 'SINTAX classifier, bayesian and global alignment taxonomy is the same'
}

# write taxonomy method dictionary out to .json file
with open(tax_method_path, 'wt') as json_out:
    json.dump(method_dict, json_out)


## FILTER ASVS BY TAXONOMY #############################################################################################

# keep only ASVs that contain this taxonomic designation in their amptk taxonomy string
wanted_tax = 'Thelephorales'

# create a list of ASV sequence records that match the wanted taxonomic group
wanted_tax_records = []
wanted_asv_ids = []

num_all_tax = 0
# go through post-amptk sequence fasta that has the taxonomy in the header
for record in SeqIO.parse(tax_seq_path, 'fasta'):
    num_all_tax += 1
    # add the ASV sequences record to the list if it is within the wanted taxonomic group
    if re.search(wanted_tax, record.description, re.I):
        wanted_asv_ids.append(record.id.split(OTU_HEADER_DELIM)[0])
        wanted_tax_records.append(record)
    else:
        continue

# how many ASVs matched the wanted taxonomic group?
if len(wanted_tax_records) > 0:
    print(f'{len(wanted_tax_records)} ASVs matched the taxonomic group {wanted_tax}.\n')
else:
    # stop the script from running (pause)
    raise KeyboardInterrupt(f'No records were located matching the taxonomic group {wanted_tax}.\n')

# write the matching ASV records to a new sequence file
filtered_tax_path = (main_output / f'{filtered_tax_prefix}_{wanted_tax}').with_suffix('.fasta')
with open(filtered_tax_path, 'w') as tax_out:
    SeqIO.write(wanted_tax_records, tax_out, 'fasta')

## READ IN OTU TABLES ##################################################################################################

# due to the size of these dataframes, use a custom function that will read in one row at a time, and return only what
#   the wanted values from the dataframe
def parse_otu_table(otu_table_path, rows_per_iter=1, asvs_to_keep='all', samples_to_keep='all', asv_col_name='#OTU ID', min_read_count=1):

    # READ IN DATAFRAME #

    # create a TextFileReader generator that will return a single row at a time, with headers, to reduce memory usage
    table_iter = pd.read_csv(otu_table_path, chunksize=rows_per_iter, delimiter='\t')


    # PROCESS DATAFRAME BY ROW #

    # create an output dictionary; primary key will be the sample ID, and its value is a dictionary of key/value pairs
    #    of each ASV that was present in this sample with its read count in this sample, as long as the read count
    #    is at or exceeds the minimum read count specified by min_read_count
    table_dict = {}

    # for each row (ASV) in the input dataframe...
    for read_counts in table_iter:

        # FILTER BY ASV #

        # get the ASV ID for the current row
        current_asv = read_counts[asv_col_name].iloc[0]

        # return all ASVs in table
        if asvs_to_keep == 'all':

            # if all ASVs are to be returned, do no filtering here and move to next section to filter by sample ID
            pass

        # return row of a single ASV
        elif isinstance(asvs_to_keep, str):

            # if this row matches the single ASV ID provided...
            if current_asv == asvs_to_keep:

                # move to the next section to filter by sample ID
                pass

            # if this row does not match the single ASV ID...
            else:

                # do nothing, continue to the next row in the table
                continue

        # return rows from a list of ASVs
        elif isinstance(asvs_to_keep, list):

            # if this current row ASV is in the list of ASVs to keep...
            if current_asv in asvs_to_keep:

                # move to the next section to filter by sample ID
                pass

            # if this ASV is not in the list of ASV IDs to keep...
            else:

                # do nothing, continue to next row in table
                continue

        # if the input value for asvs_to_keep cannot be interpreted...
        else:

            # throw an error
            error_msg = (f'ERROR. Input value for {parse_otu_table.__name__} argument asvs_to_keep must be one '
                         f'of three options: '
                         f'   (1) \'all\'\n'
                         f'   (2) a string matching a single sample ID\n'
                         f'   (3) a list of strings matching a sample ID\n')
            return exit_process(error_msg)


        # FILTER BY SAMPLE ID #

        # if all sample IDs are to be returned...
        if samples_to_keep == 'all':

            # do nothing, continue to sorting into dictionary
            pass

        # return ASV read counts for a single sample
        elif isinstance(samples_to_keep, str):

            # only climush samples; will only include samples, no mock community, etc.
            if samples_to_keep == 'climush':

                # tried on one line, but it gets angry?
                climush_ids = [col for col in read_counts.columns if re.search(STRICT_CLIMUSH_SAMPLE_ID_RE, col, re.I)]
                read_counts = read_counts[climush_ids]

            # only heather's oak savannah samples
            elif samples_to_keep == 'heather':

                heather_ids = [col for col in read_counts.columns if re.search(heather_sample_re, col, re.I)]
                read_counts = read_counts[heather_ids]

            # only jeremy's competition study samples
            elif samples_to_keep == 'cams':

                cams_ids = [col for col in read_counts.columns if re.search(cams_sample_re, col, re.I)]
                read_counts = read_counts[cams_ids]

            # only carolyn's alaska alder samples
            elif samples_to_keep == 'alaska_alders':

                alaska_ids = [col for col in read_counts if re.search(alaska_sample_re, col, re.I)]
                read_counts = read_counts[alaska_ids]

            # if the input string doesn't match a known project name...
            else:

                # assume it is a sample ID and filter to include only this sample
                read_counts = read_counts[samples_to_keep]

        # return ASV read counts for a list of sample IDs
        elif isinstance(samples_to_keep, list):

            # filter this row to include only the sample columns in the samples_to_keep list
            read_counts = read_counts[cols_to_keep]

        # if the input value for samples_to_keep cannot be interpreted...
        else:

            # throw an error
            error_msg = (f'ERROR. Input value for {parse_otu_table.__name__} argument samples_to_keep must be one '
                         f'of three options: '
                         f'   (1) \'all\'\n'
                         f'   (2) a string matching a single sample ID\n'
                         f'   (3) a list of strings matching a sample ID\n')
            return exit_process(error_msg)


        # DROP SAMPLES BELOW MIN READ COUNT #

        # create a list of sample IDs for those that meet or exceed minimum read count for this ASV
        mincount_samples = [s for s in read_counts.columns if (s != asv_col_name) and (read_counts[s].iloc[0] >= min_read_count)]  # idk why I couldn't think of better way

        # filter this row so only samples at or above minimum read count are retained
        read_counts = read_counts[mincount_samples]  # no longer has ASV ID column, not needed since we stored ASV ID in variable


        # APPEND FILTERED ROW TO OUTPUT DICT #

        # use the sample ID as the primary key, with its value being a dict of ASV IDs and the read counts for
        #    each of these ASVs in this sample (only if read count is at or above minimum value)

        # go through each sample ID, adding or updating to the output dictionary
        for samp in mincount_samples:

            # get the read count of this ASV in this sample
            asv_read_count = read_counts[samp].iloc[0]

            # if the sample ID is already in the dictionary...
            if samp in table_dict.keys():

                # no need to create a new primary key, continue to updating existing key
                pass

            # if the sample ID is not yet in the dictionary...
            else:

                # add the sample ID as a primary key, with an empty dict as a value
                table_dict[samp] = {}


            # update the sample ID with key/value pair of ASV ID and read count
            table_dict[samp].update({current_asv: asv_read_count})

            # continue to next row (ASV) in the input dataframe

    return table_dict

# CLUSTER ##

# this is the larger of the two, at 30 GB; currently not working even with this function
# clust_dict = parse_otu_table(otu_table_path=clust_tab_path, asvs_to_keep=wanted_asv_ids, samples_to_keep='climush')

# DENOISE ##
# this is also not loading
# denoi_dict = parse_otu_table(otu_table_path=denoi_tab_path, asvs_to_keep=wanted_asv_ids, samples_to_keep='climush')

# TAX-FILTERED OTU TABLE #

# I think I have to filter the full sequence fasta file before clustering then use this to cluster?
# create a list of sequences in the target taxonomic group
wanted_seq_list = []
with open(filtered_tax_path, 'r') as tax_in:
    for record in SeqIO.parse(tax_in, 'fasta'):
        wanted_seq_list.append(str(record.seq))

# use this sequence list to retain only those that match in this list in the full query fasta
query_fasta_path = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/otu-tab_test/cluster-size/from-tax-subset/query_isl-2407_renamed.fasta')
filtered_query_records = []
total_query_records = 0
with open(query_fasta_path, 'r') as query_in:
    for record in SeqIO.parse(query_in, 'fasta'):
        total_query_records += 1
        if str(record.seq) in wanted_seq_list:
            filtered_query_records.append(record)
        else:
            continue

print(f'The {len(wanted_seq_list)} {wanted_tax} ASVs were detected in {len(filtered_query_records)} samples in '
      f'this bioinformatics run.\n')

filtered_query_out = query_fasta_path.parent / 'filtered-query_isl-2407.fasta'
with open(filtered_query_out, 'w') as query_out:
    SeqIO.write(filtered_query_records, query_out, 'fasta')

# compare the number of input sequences for vsearch --cluster_size, in filtered-query above, matches its output
filtered_otu_out = query_fasta_path.parent / 'tax-filt_otus.fasta'
num_filtered_out = 0
for record in SeqIO.parse(filtered_otu_out, 'fasta'):
    num_filtered_out += 1
print(f'Number of input sequences: {len(filtered_query_records)}\n'
      f'Number of output sequences: {num_filtered_out}\n')



## MERGE OTU TABLE WITH FILTERED TAXONOMY ####


# how many clusters are there, based on the .uc output from vsearch --cluster-size
uc_path = tax_subset_path / 'tax-filt2.uc'
uc_df = pd.read_table(uc_path, delimiter='\t')
cluster_count = max(uc_df.iloc[:,1]) + 1
print(f'Number of OTUs in the .uc output = {cluster_count}')

# how many clusters are there, based on the OTU .fasta file output from vsearch --cluster-size
otus_path = tax_subset_path / 'tax-filt2_otus.fasta'
otus_count = 0
for record in SeqIO.parse(otus_path, 'fasta'):
    otus_count += 1
print(f'Number of OTUs in the OTU .fasta output = {otus_count}')


# check number of OTUs in second run of --cluster_size
tax2_path = tax_subset_path / 'tax-filt2_otu-table.txt'
tax2_otu_count = pd.read_table(tax2_path, delimiter='\t').shape[0]
print(f'Number of OTUs in the OTU table output = {tax2_otu_count}')



# read in produce of --cluster_size that used filtered-query_isl-2407.fasta as input
tax_filt_df = pd.read_table(tax2_path, delimiter='\t')

# rename the #OTU ID column to read_header, to match the query df and to better describe contents in column
tax_filt_df.rename(columns={'#OTU ID': 'read_header'}, inplace=True)

# add the sequence and taxonomy from the filtered fasta query file to the otu table

# create a dictionary of headers and sequences from the query fasta
query_seqs = {'read_header': [],
              'sequence': []}
for record in filtered_query_records:
    query_seqs['read_header'].append(record.id.split(OTU_HEADER_DELIM)[0])
    query_seqs['sequence'].append(str(record.seq))

# add sequences to the otu table by matching the record ids (headers), which are currently the OTU IDs

# convert dict from query fasta to a dataframe
query_seqs_df = pd.DataFrame.from_dict(query_seqs, orient='columns')

# check that two df are same length before merging
if query_seqs_df.shape[1] == tax_filt_df.shape[1]:
    # merge with tax filt dataframe, on #OTU ID
    merged_tax_df = query_seqs_df.merge(tax_filt_df, on='read_header')
else:
    # create a set of unique values in the read header cols
    query_unique = set(query_seqs_df['read_header'])
    tax_unique = set(tax_filt_df['read_header'])

    # IN QUERY, MISSING FROM FILTERED TAX
    diff_headers = query_unique.difference(tax_unique)
    percent_diff = (len(diff_headers) / len(query_unique)) * 100

    print(f'{len(diff_headers)} read header(s) ({percent_diff:.1f}%) are in the query dataframe but not the '
          f'filtered taxonomy dataframe. Would you like to see these read headers?\n')
    show_headers = input(f'y/n >')
    if show_headers == 'y':
        print('\n  '.join(diff_headers))
    else:
        pass

    # IN FILTERED TAX, MISSING FROM QUERY
    diff_headers = tax_unique.difference(query_unique)
    percent_diff = (len(diff_headers) / len(tax_unique)) * 100

    print(f'{len(diff_headers)} read header(s) ({percent_diff:.1f}%) are in the filtered taxonomy dataframe but not the '
          f'query dataframe. Would you like to see these read headers?\n')
    show_headers = input(f'y/n >')
    if show_headers == 'y':
        print('\n  '.join(diff_headers))
    else:
        pass

    raise KeyboardInterrupt('Read headers don\'t match between query and taxonomy dataframes. Fix before continuing.')


# drop non-climush samples; wait until after merging as it is easiest
climush_cols = [c for c in merged_tax_df.columns if re.search(STRICT_CLIMUSH_SAMPLE_ID_RE, c, re.I)]
climush_cols.insert(0, 'read_header')  # add in the non-sample columns too, which we want
climush_cols.insert(1, 'sequence')
merged_tax_df = merged_tax_df[climush_cols]


## EXPLORE BIOM OBJECTS ################################################################################################

# NOTE - i started this before switching to small chunksize values with pandas

# denoi_biom = biom.load_table(denoi_tab_path)
# clust_biom = biom.load_table(clust_tab_path)

# biom format Tables

# get the number of ASVs (rows) and number of samples (columns)
# num_asvs, num_samples = denoi_biom.shape


## GET TAXONOMY ########################################################################################################

# create a dict where the key is the ASV ID and the value is the full taxonomy string
tax_dict = {}
not_found_error = []

# go through each ASV in the taxonomy OTU .fasta file...
for tax_record in SeqIO.parse(tax_fasta_path, 'fasta'):

    # get the ASV ID from the header, without the read abundance included (size=)
    asv_id = tax_record.id.split(OTU_HEADER_DELIM)[0]

    # use regular expressions to pull information from the taxonomy string in the header

    # regex that will look for any of the amptk taxonomy assignment methods, and take everything after
    amptk_tax_all_re = AMPTK_TAX_METHOD_RE + '.+'

    # use regex to get amptk taxonomy information
    try:
        tax_str = re.search(amptk_tax_all_re, tax_record.description, re.I).group(0)

        # add asv ID and taxonomy info to dictionary
        tax_dict.update({asv_id: tax_str})

    # if regex doesn't return a match, then print a warning and append ASV ID to error list
    except AttributeError:
        print(f'Taxonomy string not found: {asv_id}')
        not_found_error.append(asv_id)




# check for any errors, pause if any detected
if len(not_found_error) > 0:
    raise KeyboardInterrupt(f'Some taxonomy strings could not be detected, check the error '
                            f'list to see which ASVs had issues.\n')
else:
    pass


## SORT TAXONOMY VALUES ################################################################################################

## GET TEST VALUES ##
def get_tax_w_sp(taxonomy_dict, tax_lvl):

    regex_dict = {'kingdom': r'(?<=k:)[a-zA-Z]+',
                  'phylum': r'(?<=p:)[a-zA-Z]+',
                  'class': r'(?<=c:)[a-zA-Z]+',
                  'order': r'(?<=o:)[a-zA-Z]+',
                  'family': r'(?<=f:)[a-zA-Z]+',
                  'genus': r'(?<=g:)[a-zA-Z]+',
                  'species': r'(?<=s:)\w.+\b'}

    for asv_id, tax in taxonomy_dict.items():
        if re.search(regex_dict[tax_lvl], tax, re.I):
            yield asv_id, tax
        else:
            continue
test_asv, test_tax_list = next(get_tax_w_sp(tax_dict, tax_lvl='species'))
#####################

# create a function that will run try / except for the given value
## THIS SHOULD HAVE **KWARGS
def try_to_find(search_re, search_str, add_to_dict, output_dict, dict_key, case_insensitive=True):

    # try to use regex to find match and return the match string
    try:
        if case_insensitive:
            match_str = re.search(search_re, search_str, re.I).group(0)
        else:
            match_str = re.search(search_re, search_str).group(0)

    # if the regex doesn't return a match, return None
    except AttributeError:
        return None

    # if you want to add output to dictionary, then following steps are executed
    if add_to_dict:
        output_dict[dict_key] = match_str
        return output_dict  # can do something with the updated dict, or it will just update in place too

    # if you don't want to add it to a dictionary, the matching string will be returned (if found)
    else:
        return match_str

def add_tax(tax_lvl, search_str, output_dict):

    regex_dict = {'kingdom': r'(?<=k:)[a-zA-Z]+',
                  'phylum': r'(?<=p:)[a-zA-Z]+',
                  'class': r'(?<=c:)[a-zA-Z]+',
                  'order': r'(?<=o:)[a-zA-Z]+',
                  'family': r'(?<=f:)[a-zA-Z]+',
                  'genus': r'(?<=g:)[a-zA-Z]+',
                  'species': r'(?<=s:)\w.+\b'}

    # search for tax level in string
    tax_result = re.search(regex_dict[tax_lvl], search_str, re.I)

    # if a match is found, add it to the output dictionary
    if tax_result:
        output_dict[tax_lvl] = tax_result.group(0)

    # if no match is found, add 'NA' for this tax level
    else:
        output_dict[tax_lvl] = 'NA'

    return output_dict


# create a list of the new taxonomy columns that will be added to the OTU table
new_tax_columns = ['taxonomy_method', 'accession_number', 'percent_match',
                   'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# create a list of just the taxonomy lvl column names, used in loop below
tax_only_cols = [c for c in new_tax_columns if not '_' in list(c)]

# create a list for any ASVs that have issues with pulling information from the taxonomy string
tax_str_error = []

# create a dictionary to add sorted values to, with key as ASV ID and values the sorted taxonomic information
new_tax_dict = {k:{} for k in tax_dict.keys()}

# go through each of the ASVs and their taxonomy strings
for asv_id, taxonomy in tax_dict.items():

    # create a dictionary for each of the new taxonomy columns, will append values to this
    tax_column_values = {k:'' for k in new_tax_columns}

    ## TOP MATCH DETAILS #################################################################

    # TAXONOMY METHOD

    # taxonomy method regex, must be at start of taxonomy string
    tax_method_re = '^' + AMPTK_TAX_METHOD_RE

    # find the taxonomy method for this ASV
    try_to_find(search_re=tax_method_re, search_str=taxonomy, add_to_dict=True,
                output_dict=tax_column_values, dict_key='taxonomy_method')

    # ACCESSION NUMBER

    # no good way to use regex here, variable num of accession nums (1-2) both divided by | and ends in num
    access_num_str = taxonomy.split(';')[0].split(TAX_STR_DELIM)[2:]

    # add the list of accession numbers to the dict
    tax_column_values['accession_number'] = access_num_str

    # PERCENT MATCH

    # percent match regex
    perc_match_re = r'(?<=|)\d{1,3}\.\d{1,4}(?=|)'

    # search for percent match, add to dictionary
    try_to_find(search_re=perc_match_re, search_str=taxonomy, add_to_dict=True,
                    output_dict=tax_column_values, dict_key='percent_match')

    ## TAXONOMIC ASSIGNMENTS ##################################################################

    # taxonomic assignments and top match details are separate by a semicolon
    tax_assign = taxonomy.split(';')[-1]

    for tax_lvl in tax_only_cols:
        add_tax(tax_lvl=tax_lvl, search_str=tax_assign, output_dict=tax_column_values)


    ## UPDATE NEW COLUMN DICT #################################################################

    new_tax_dict.update({asv_id: tax_column_values})


## CREATE TAXONOMY TABLE ###############################################################################################

# create a pandas dataframe from the taxonomy dictionary, format to match the OTU table
tax_df = pd.DataFrame(new_tax_dict).transpose().reset_index().rename(columns={'index':'#OTU ID'})


## MERGE ASV AND TAXONOMY TABLES #######################################################################################

# read in the OTU table as a pandas dataframe
asv_df = pd.read_table(otu_table_path, delimiter='\t')

# merge the OTU table and the taxonomy table, arranging so information from the taxonomy dataframe is first
tax_asv_df = pd.merge(left=tax_df, right=asv_df, how='inner', on='#OTU ID')


## SUBSET NON-CLIMUSH SAMPLES ##########################################################################################

# Heather's samples

# subset heather's horseshoe table, then merge with tax_df
heather_samp_cols = [c for c in asv_df.columns if c.startswith('illumina_Horse')]
heather_keep_cols = ['#OTU ID'] + heather_samp_cols
heather_otus = pd.merge(left=tax_df, right=asv_df[heather_keep_cols], how='inner', on='#OTU ID')

# sum read count of each ASV
heather_otus['total_read_count'] = heather_otus[heather_samp_cols].sum(axis=1)

# filter any ASVs with zero read count; i.e., remove ASVs not in Heather's samples
heather_filtered = heather_otus[heather_otus['total_read_count'] > 0]

# filter the ASV .fasta to only include the ASVs in Heather's OTU table
# get ASV IDs in Heather's samples
heather_asv_ids = heather_filtered['#OTU ID'].to_list()

# generator to get only certain ASVs from
def filter_otu_fasta(full_otu_fasta, subset_list):
    for record in SeqIO.parse(full_otu_fasta, 'fasta'):
        if record.id.split(OTU_HEADER_DELIM)[0] in subset_list:
            yield record
        else:
            continue

if heather_filtered.shape[0] == len(filter_otu_fasta(tax_fasta_path, heather_asv_ids)):
    # write out filtered table to a .csv OTU table
    heather_filtered.to_csv(heather_tab_out, index=False)
    # write out filtered sequence records to OTU .fasta
    SeqIO.write(filter_otu_fasta(tax_fasta_path, heather_asv_ids), heather_otu_out, 'fasta')
elif heather_filtered.shape[0] < len(filter_otu_fasta(tax_fasta_path, heather_asv_ids)):
    print(f'ERROR. ')

## CREATE CLIMUSH THELEPHORALES FASTA ##################################################################################

# create a table with only climush samples

# all sample columns start with illumina; only climush samples have litter or soil as the second string
# climush_sample_cols = [c for c in tax_asv_df.columns if re.search(r'^illumina_soil|^illumina_litter', c, re.I)]
# taxonomy_cols = [c for c in tax_asv_df.columns if not re.search(r'^illumina', c, re.I)]
#
# # create a list of all columns to keep from original dataframe (taxonomy and climush samples only)
# cols_to_keep = taxonomy_cols + climush_sample_cols
#
# # subset the original dataframe to create a new dataframe with wanted columns only
# climush_df = tax_asv_df[cols_to_keep]
#
# # calculate the total read count for each ASV, add as new column at end of table
# climush_df['total_read_count'] = climush_df[climush_sample_cols].sum(axis=1)
#
# # check the min and max read counts, print out the taxonomy columns for these
#
# # max
# max_asv = climush_df[climush_df['total_read_count'] == climush_df['total_read_count'].max()]
# max_read_count = max_asv['total_read_count'].iloc[0]
# max_tax = max_asv[taxonomy_cols]
#
# # print(f'The top ASV is:'
# #       f'   {}'
# #       f'with a total read count of:'
# #       f'   {max_read_count}\n')
#
# # min
# min_read_count = climush_df['total_read_count'].min()
#
#
# #######
#
# # print sample IDs in fasta file
# t = Path('/Users/carolyndelevich/main/github_repos/climush/bioinformatics-pipeline/pipeline-output/otu-table/query_sm.fasta')
#
# s = []
# with open(t, 'rt') as fin:
#     for line in fin.readlines():
#         if line.startswith('>'):
#             s.append(line.split(';')[0].replace('>', ''))
#
# print('\n'.join(s))
#
#
