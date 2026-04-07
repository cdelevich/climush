from pathlib import Path
from Bio import SeqIO
import re, argparse
from datetime import datetime
from climush.utilities import convert_udb_format, exit_process


TEST_MODE = False
if TEST_MODE:
    # less stable but works interactively (i.e., testing in Python console in PyCharm)
    prog_name = 'pull-reference-seqs'
    print(f'WARNING - {prog_name} RUNNING IN TEST MODE')

    # test arguments
    test_input = [
        '--db', 'parent/directory/of/ITS.udb',
        '--udb', 'SH1174631.09FU',
        # '--gb','K12345'
        '--tax', 'Amanita_exitialis',
    ]

else:
    # more stable but does not work interactively
    prog_name = Path(__file__).stem

## COMMAND LINE OPTIONS ################################################################################################

## INSTANTIATE PARSER ##

parser = argparse.ArgumentParser(
    prog=prog_name,
    description='Pull reference sequences from a sequence reference database.',
    epilog='This script is part of the CliMush bioinformatics pipeline.'
)


## INPUT / OUTPUT ##

# REQUIRED; path to the sequence reference database / fasta file to search for reference sequences in
parser.add_argument(
    '-d', '--db',
    type=Path,
    required=True,
    help='File path to a reference database or .fasta file from which to pull reference sequences from '
         'matching the input search parameters'
)

# OPTIONAL; path to the output .fasta file that contains the reference sequences pulled from -d / --db
parser.add_argument(
    '-o', '--output',
    type=Path,
    required=False,
    default=None,
    help='File path to the output .fasta file that contains the subsetted reference sequences from the '
         'input reference database provided to the -d / --db argument. If nothing is provided, the '
         'search results will be dated and numbered by the order in which the search occurred for that '
         'day.'
)


## INPUT FILE HEADER FORMAT ##

# OPTIONAL; the delimiter in the read headers that separate the accession numbers and taxonomy info
parser.add_argument(
    '--delim',
    type=str,
    nargs='?',
    required=False,
    default=';',
    help='The delimiter used in the reference sequence headers to separate the accession numbers from the taxonomy. '
         'In the UNITE reference database format, a semi-colon delimits these two fields of information (default).',
)


## SEARCH PARAMETER SUBARGUMENTS ##

# set a group for the following search arguments
search_args = parser.add_argument_group(
    title='search parameters',
    description='the search parameters used to pull matching reference sequences from input reference sequence file',
)

# search by taxonomy string
search_args.add_argument(
    '-t', '--tax',
    type=str,
    nargs='*',
    required=False,
    help='A taxonomic string (e.g., species binomial, class) to pull reference sequences for. Any sequences '
         'with taxonomy matching this string will be returned. If searching for a match to a species binomial, '
         'include an underscore (_) between the genus and species.',
)

# search by UNITE database reference number
search_args.add_argument(
    '-u', '--udb',
    type=str,
    nargs='*',
    required=False,
    help='The UNITE database reference number of the reference sequence(s) to pull. This should start with '
         'either SH for a species hypothesis accession number or UDB a INSDC / UNITE accession number of '
         'the representative / reference sequence of a species hypothesis.',
)

# search by GenBank database reference number
search_args.add_argument(
    '-g', '--gb',
    type=str,
    nargs='*',
    required=False,
    help='The GenBank accession number of the reference sequence(s) to pull.',
)

## PARSE ARGUMENTS INTO DICTIONARY ##

if TEST_MODE:
    args = vars(parser.parse_args(test_input))
else:
    args = vars(parser.parse_args())

########################################################################################################################


## CHECK COMMAND LINE INPUT ############################################################################################

## INPUT / OUTPUT ##

# if an output file path isn't provided, create one
if args['output'] is None:

    # create a date suffix with today's date to add to the output file name
    date_suffix = datetime.today().strftime('%Y-%m-%d')

    # get the name of the searched database to use in the output file name
    search_db_name = args['db'].stem

    # check if there was already a search today
    output_files_today = [ file for file in args['db'].parent.glob(f'*{date_suffix}*')]
    # if there was already a search done today...
    if len(output_files_today) > 0:

        # regex to find the search number in an output file
        search_num_re = r'(?<=search)\d{2}(?=_)'

        # iterate through todays auto-named output files...
        search_nums = []
        for today_file in output_files_today:

            # search for the search number of each file...
            search_num_found = re.search(search_num_re, today_file.name)

            # once located...
            if search_num_found:
                # add the number of the search to the list of search numbers already used today (as an integer)
                search_nums.append(int(search_num_found.group(0)))

        # use the next number in the series of search numbers for today for this file
        next_search_num = max(search_nums)

        # if the number is only a single digit, add a leading zero
        if next_search_num < 10:
            days_search_num = '0' + str(next_search_num)

        # if double-digits, don't add leading zero
        else:
            days_search_num = str(next_search_num)

    # if this was today's first search with an auto-generated output file path...
    else:
        # use 01 as the search of today number
        days_search_num = '01'

    # assemble output file's basename into a single string
    output_basename = '_'.join([search_db_name, 'search' + days_search_num, date_suffix])

    # replace None value with created output file path
    args['output'] = (args['db'].parent / output_basename).with_suffix('.fasta')

# if an output file path is provided, then use this (do nothing here)
else:
    pass


## SEARCH PARAMETER SUBARGUMENTS ##

# check that at least one of the arguments in the search_args argument group are provided

# if all are None, raise Error
if (args['tax'] is None) and (args['udb'] is None) and (args['gb'] is None):

    # create error message
    err_msg = (f'At least one search option needs to be provided.')

    raise parser.error(err_msg)

# otherwise, check each format to make sure the input is valid
else:

    if args['udb'] is not None:

        # must start with UDB or SH
        for udb_ref in args['udb']:
            if udb_ref.startswith('UDB') or udb_ref.startswith('SH'):
                pass

            # otherwise, raise error
            else:
                err_msg = f'Input to the argument -u / --udb must start with UDB or SH'
                parser.error(err_msg)

    else:
        pass

    if args['gb'] is not None:

        for gb_ref in args['gb']:
            # must start with letter and have at least 5 digits
            valid_gb_match = re.search(r'^[A-Z].+?\d{5,}', gb_ref)
            if valid_gb_match:
                pass
            else:
                err_msg = (f'A valid GenBank accession number, as provided to -g / --gb, should '
                           f'start with a letter and contain at least 5 numbers.')
                parser.error(err_msg)

    else:
        pass



########################################################################################################################


## FORMAT INPUT REFERENCE SEQ FILE #####################################################################################

# check the input file format, and if .udb, return a .fasta-formatted version
input_refseqs, udb_converted = convert_udb_format(
    input_file=args['db'],
)

# if convert_udb_format locates an extracted.fa format, it will use this, but SeqIO doesn't like .fa over .fasta
# check the file format, and if .fa, create a variable for the file format set to .fasta for use with SeqIO
input_refseqs_fmt = input_refseqs.suffix
if input_refseqs_fmt == '.fa':
    input_refseqs_fmt = '.fasta'
else:
    pass

########################################################################################################################


## FORMAT SEARCH PARAMETER STRINGS #####################################################################################

# function will any OR regex searches
def compile_search_regex(search_list, check_for_binomial):

    if search_list is None:
        return False, re.compile(pattern='')

    else:

        if check_for_binomial:

            # create an empty list to add formatted search strings to
            search_list_fmt = []

            # iterate through the search_list looking for a species binomial, which should have an underscore (_)
            for search_str in search_list:

                # if an underscore is located within the string, format regex for a binomial (two-word match)
                if re.search('_', search_str):

                    # first replace the underscore with a space
                    # then add additional regex for positive look-behind and terminal string match
                    search_str_fmt = '(?<=s:)' + search_str.replace('_', ' ') + '$'

                    # add this species binomial formatted search regex to the list of search strings formatted for regex
                    search_list_fmt.append(search_str_fmt)

                # if there is no underscore within the string, add just start and end formatting
                else:

                    # add to list without additional regex formatting
                    search_list_fmt.append(f'^{search_str}$')

            # once all strings in the search list have been checked for binomials and reformatted when needed...

            # return the list as a single regex string, searching for matches to any of the components
            return True, re.compile('|'.join(search_list_fmt))


        # if used on a accession number search list, no need to check for species binomials to reformat
        else:
            return True, re.compile('|'.join(search_list))



## CREATE REGEX COMPILER DICTIONARY ##

# sort the argparse input related to search strings into a dictionary
search_dict = {
    'taxonomy': {
        'requires_search': compile_search_regex(search_list=args['tax'], check_for_binomial=True)[0],
        'compiler': compile_search_regex(search_list=args['tax'], check_for_binomial=True)[1],
    },
    'accession_num': {
        'unite': {
            'requires_search': compile_search_regex(search_list=args['udb'], check_for_binomial=False)[0],
            'compiler': compile_search_regex(search_list=args['udb'], check_for_binomial=False)[1],
        },
        'genbank': {
            'requires_search': compile_search_regex(search_list=args['gb'], check_for_binomial=False)[0],
            'compiler': compile_search_regex(search_list=args['gb'], check_for_binomial=False)[1],
        }
    }
}

########################################################################################################################


## SEARCH REFERENCE SEQUENCES ##########################################################################################

# create a function that will check a search group
def search_seqrecord(seqrecord, search_category, search_info_dict):

    ## CHECK FOR VALID SEARCH CATEGORY ##

    # get a list of the primary keys in the input dictionary, which are the search categories
    valid_primary_categories = list(search_dict.keys())

    # can provide a dictionary to search_category if it is a multilevel category (e.g., accession_num > genbank)
    if isinstance(search_category, dict):

        # if the input provided to search_category is a dictionary, get the secondary key
        valid_secondary_categories = list(search_category.values())

        # pull the primary (key) and secondary (value) categories from the input list
        primary_category = list(search_category.keys())[0]
        secondary_category = list(search_category.values())[0]

        # confirm that the search category provided is a primary key in the input dictionary
        if primary_category in valid_primary_categories:
            if secondary_category in valid_secondary_categories:
                pass
            else:
                err_msg = (f'The secondary category, {secondary_category}, is not a valid secondary search '
                           f'category: {valid_secondary_categories}')
                return exit_process(err_msg)
        else:
            if secondary_category in valid_secondary_categories:
                err_msg = (f'The primary category, {primary_category}, is not a valid primary search category: '
                           f'{valid_primary_categories}.')
                return exit_process(err_msg)
            else:
                err_msg = (f'Neither the primary category, {primary_category}, nor the secondary search category, '
                           f'{secondary_category}, are valid search categories:\n'
                           f'   valid primary search categories: {valid_primary_categories}\n'
                           f'   valid secondary search categories: {valid_secondary_categories}\n')
                return exit_process(err_msg)


        # if all tests pass for dictionary version of search_category input, subset dictionary
        search_category = search_info_dict[primary_category][secondary_category]

    elif isinstance(search_category, str):

        primary_category = search_category

        # confirm that the search category provided is a primary key in the input dictionary
        if search_category in valid_primary_categories:

            # subset dictionary
            search_category = search_info_dict[search_category]

        # if it is not a valid primary key, return an error and exit
        else:
            err_msg = f'{search_category} is not a valid search category: {valid_search_categories}\n'
            return exit_process(err_msg)

    else:
        raise ValueError('The search category needs to be either a dictionary or string.')


    ## CREATE LIST OF ITEMS IN REFSEQ TO SEARCH ##

    if primary_category == 'accession_num':
        refseq_search_list = seqrecord.description.split(args['delim'])[0].split('|')
    elif primary_category == 'taxonomy':
        refseq_search_list = [tax.split(':')[-1] for tax in seqrecord.description.split(args['delim'])[1].split(',')]
    else:
        err_msg = 'issue detecting input errors'
        return exit_process(err_msg)

    ## SEARCH REFSEQ LIST FOR MATCHING ITEMS ##

    # check if this search category was provided a string / regex to search with
    if search_category['requires_search']:

        # if it does contain a string to use in search, use its compiler to perform the search
        match_list = list(filter(search_category['compiler'].match, refseq_search_list))

        # if a match is found...
        if len(match_list) > 0:

            return True, seqrecord

        # if a match is not found, return None
        else:
            pass

    # if this search category does not have a string / regex to use to search with, return None
    else:
        pass

    return False, ''

# create a list of matching sequence records
matching_seqrecords = []

# iterate through each record in the .fasta-formatted reference sequence file
with open(input_refseqs, 'r') as refseqs_in:

    # for each reference sequence in the input reference sequence database...
    for refseq in SeqIO.parse(refseqs_in, input_refseqs_fmt.replace('.', '')):

        ## SEARCH FOR ACCESSION NUMBER ##

        # use search_seqrecord() to search the seqrecord for match to argparse input search strings, will add seqrecord
        #   to list if match is found

        for category in ['taxonomy', {'accession_num': 'unite'}, {'accession_num': 'genbank'}]:

            match_found, seqrecord_match = search_seqrecord(
                seqrecord=refseq,
                search_category=category,
                search_info_dict=search_dict,
            )

            # as soon as a match to any search string is found, stop searching for others, move onto next refseq
            if match_found:
                # add the matching seqrecord to the list of matching seqrecords
                matching_seqrecords.append(seqrecord_match)

                # stop searching this refseq and move on to next refseq
                break

            # if no match found, try next search regex (if any)
            else:
                continue

########################################################################################################################


## WRITE MATCHES OUT TO FILE ###########################################################################################

# only write out file if matches were found:
if len(matching_seqrecords) > 0:
    # write the matching sequence records out to the output file path
    with open(args['output'], 'w') as matchseqs_out:
        SeqIO.write(matching_seqrecords, matchseqs_out, format='fasta')

    # confirm file was created and print number of matching records:
    if args['output'].is_file():
        num_matches = len(matching_seqrecords)
        print(f'{num_matches} sequence records in the reference dataset {input_refseqs.name} were matched to the following '
              f'search parameters:\n'
              f'   ...\n'
              f'and written to the following output file:\n'
              f'   {args["output"]}\n')

# if no matching records found, then print message
else:
    print(f'No sequence records in the reference dataset {input_refseqs.name} were matched to the following '
          f'search parameters:\n'
          f'   ...\n')

########################################################################################################################