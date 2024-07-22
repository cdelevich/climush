import pathlib, argparse, re
from pathlib import Path
from climush.utilities import get_sample_id

# DEFINE COMMAND LINE OPTIONS #######################################################################################

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Export a sorted list of sample names in the input directory.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

parser.add_argument('-i', '--input',
                    type=pathlib.PosixPath,
                    required=True,
                    help='The path containing directories or files that have the sample ID within the directory or '
                         'path names.')

parser.add_argument('-o', '--output',
                    type=pathlib.PosixPath,
                    required=False,
                    help='The name to use for the output file containing the sorted sample IDs in the input '
                         'directory')

parser.add_argument('--search',
                    required=False,
                    nargs='*',  # gathers all arguments after into a list, regardless of how many
                    help='The substring to search for in the input directory. For example, if you want a list of '
                         'only the samples from Florida, you can use the substring \'ORD\' for PacBio reads')

# PARSE COMMAND LINE OPTIONS ######################################################################################

# create a dictionary of arg names and input values
args = vars(parser.parse_args())

# if an output file name is not provided, create one based off the input directory name
if args['output'] is None:
    output_parent = args['input']
    output_name = args['input'].name + '_sample-list'
    args.update({'output': (output_parent / output_name).with_suffix('.txt')})

# create a regex search string from the search parameter
if args['search'] is None:
    search_for = '.'  # if this arg is not provided, then search for everything in dir
else:
    search_for = '|'.join(args['search'])

# SEARCH INPUT DIRECTORY AND CREATE LIST OF SAMPLE IDS #############################################################

# create a list of sample names in the input directory that match the search string
sample_list = []
for file in list(args['input'].glob('*')):
    if re.search(search_for, file.name, re.I):  # used regex because much easier to join multiple search strings with or ('|')
        # intended to use the get_sample_id() function but it is currently shit for pacbio
        sample_id = file.stem.split('_')[-1]  # use .stem so that the file extension is not included
        sample_list.append(sample_id)


# define a function that is used to determine how the sample IDs are sorted
def sort_sampleid(filename):
    num = int(filename.replace('ORD',''))  # numbers are strings, so will not sort correctly between 2 and 3 digit nums
    return num

# sort the sample list based on this sorting function
sample_list.sort(key=sort_sampleid)

# EXPORT THE SORTED LIST INTO THE OUTPUT FILE PATH ##################################################################

with open(args['output'], 'wt') as fout:
    fout.write(f'directory: {args["input"]}\n\n')  # create line at top of file containing the input directory path
    fout.write(f'search str: {args["search"]}\n\n')  # create line underneath this containing the search string(s) used
    fout.write(f'sample names:\n')
    for s in sample_list:
        fout.write(f'  {s}\n')  # write the indented list of sorted sample IDs, one per line