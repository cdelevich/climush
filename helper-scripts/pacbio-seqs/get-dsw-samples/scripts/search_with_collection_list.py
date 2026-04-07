import argparse
import re
from pathlib import Path
from warnings import warn

## REMOVE AFTER TESTING ##

args = {
    'coll_list': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/pacbio-seqs/get-dsw-samples/data/output/climush_sporocarp-collection-list_DSW_without-sequences_2025-09-26.txt'),
    'start': Path.cwd(),
}


## COMMAND LINE OPTIONS ################################################################################################

## INSTANTIATE PARSER ###########################################################

parser = argparse.ArgumentParser(
    prog=Path(__file__).stem,
    description='Use a collection list .txt file to search for files matching'
                'the collection list names.',
    epilog='This script is part of the climush bioinformatics pipeline.',
)


## ADD ARGUMENTS ################################################################

## I/O ##

parser.add_argument(
    '--coll-list',
    type=Path,
    required=True,
    help='Path to a collection list .txt file, in which there is one collection name per line that includes '
         'the sequence run date and a unique collection number as the last value in the collection file basename.',
)

parser.add_argument(
    '--path-list',
    nargs='',
    const=,
    help='',
)

## PATH WALK SETTINGS ##

parser.add_argument(
    '--start',
    type=Path,
    required=False,
    default=Path.cwd(),
    help='The location in which to start searching for files that will match the collections in the input collection '
         'list. If none provided, the search will start in the current working directory. Starting at a location '
         'as close to the files as possible will increase the speed at which this script will produce results.',
)



## PARSE ARGUMENTS ##############################################################

# parse arguments received from the command line into a dictionary
args = vars(parser.parse_args())

########################################################################################################################


## IMPORT COLLECTION LIST ##############################################################################################

# import each sample ID from the collection list as a key in a dictionary, with an empty value
with open(args['coll_list'], 'r') as colls_in:
    collection_list = { line.strip(): '' for line in colls_in.readlines() }

########################################################################################################################


## SET PATHWALK LIMITS #################################################################################################

## GET UNIQUE SEQUENCE RUN DATES ################################################

## COMPILE SEQRUN REGEX ##

# compile a regex that matches the sequence run date in the pacbio sporocarp sample ID
seqrun_regex = re.compile(r'\d{4}-\d{2}')

## CREATE OUTPUT LIST/SET ##

# add all matches to the sequence run date regex to a set of unique run dates from the input collection list
seqrun_dates = set()
with_seqrun_date = 0

# if any samples in the collection do not match the sequence run date regex, add to a list
no_seqrun_dates = []

## SEARCH EACH SAMPLE ID ##

# go through each sample ID in the input collection list...
for sample_id in collection_list.keys():

    # use the compiled sequence run date regex to search for the sequence run date in the sample ID
    seqrun_found = seqrun_regex.search(sample_id)

    # if the regex returned a match to the sequence run date...
    if seqrun_found:

        # add the matching substring (i.e., the sequence run date) to the set of unique sequence run dates
        seqrun_dates.add(seqrun_found.group(0))

        # add to counter for number of sample IDs with sequence run date in name, since
        #   set of sequence run dates won't tell you how many samples
        with_seqrun_date += 1

    # if no match to the sequence run date regex was made...
    else:

        # add this sample ID to the list of sample IDs without a sequence run date
        no_seqrun_dates.append(sample_id)

## CHECK FOR SAMPLES W/O SEQRUN DATE ##

# summarize matches vs no matches
total_input_samples = len(collection_list.keys())
without_seqrun_date = len(no_seqrun_dates)

# warn that no sequence run date in some samples will require search of all directories within args['start']
if without_seqrun_date > 0:
    warn_msg = (f'{without_seqrun_date} out of {total_input_samples} sample IDs in the input collection list,\n'
                f'{args["coll_list"].name}, do not have a sequence run date in the sample ID '
                f'(e.g., {no_seqrun_dates[0]}).\n'
                f'Therefore, all directories will be searched, starting with the search start directory:\n'
                f'   {args["start"]}\n')
    warn(warn_msg)

# if all sample IDs in the input collection list have sequence run dates...
else:
    # do nothing here, will print a summary message further down
    pass


## COMPILE DIRECTORIES TO SEARCH ################################################

## CREATE SEQRUN DIRECTORY LIST ##

# create a list of directories to search in based on the sequence run dates from the sample IDs in the
#   input collection list
dirs_to_search = [ f'pacbio_sporocarp-f_{seqrun_date}' for seqrun_date in seqrun_dates ]

########################################################################################################################


## 2023-12 samp count
samps_2023_11 = [f for f in collection_list if f.startswith('pacbio_sporocarp-f_2023-11')]
samps_2023_12 = [f for f in collection_list if f.startswith('pacbio_sporocarp-f_2023-12')]
samps_2025_02 = [f for f in collection_list if f.startswith('pacbio_sporocarp-f_2025-02')]

sampcount_2023_11 = len(samps_2023_11)
sampcount_2023_12 = len(samps_2023_12)
sampcount_2025_02 = len(samps_2025_02)

print(f'2023-11: {sampcount_2023_11}\n'
      f'2023-12: {sampcount_2023_12}\n'
      f'2025-02: {sampcount_2025_02}\n')

total = sampcount_2023_11 + sampcount_2023_12 + sampcount_2025_02



## PATHWALK FILE SEARCH ################################################################################################

# do a top-down path walk starting at the command-line input for start (or current working directory)
for root, dirs, files in args['start'].walk():

    # if the currently searched directory is among the directories to search...
    if root.name in dirs_to_search:
        ...

    # if the cu


########################################################################################################################