from pathlib import Path
import argparse

parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Locate samples in the CliMush Sequences Globus directory.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

# required positional argument; provide a sample ID or partial ID that will convert to glob or regex
parser.add_argument('sample_id',
                    required=True,
                    help='The sample ID that you want to locate in the files; can be a full sample ID or a partial '
                         'match if a glob or regex format is provided as the sample ID string.')

# option to return only exact matches to the sample ID input string
parser.add_argument('--exact',
                    action='store_true',
                    help='Use flag if you want only an exact match to the provided sample ID string. By default, exact '
                         'matches only is turned off.')

# option to export the results as a .tsv file
parser.add_argument('--export',
                    action='store_true',
                    help='Use flag if you want the results of the search to be exported to a .tsv file. By default, '
                         'results will only be displayed through the terminal and not written to a file.')

