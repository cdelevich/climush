from pathlib import Path
import pathlib, re, warnings, argparse

TEST_MODE = False

# if using test mode
if TEST_MODE:

    # warn that script is running in test mode
    warn_msg = f'Test mode running for change-file-extension.py, so command line options are not being processed.\n'
    warnings.warn(warn_msg)

    # set default test args
    args = {
        'input': Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/endophytes/test-files_2023-endos_prefix-replaced'),
        'rep': '.fastq.gz',
    }

# if not using test mode
else:

    # instantiate the parser
    parser = argparse.ArgumentParser(
        prog=Path(__file__).stem,
        description='Replace input file\'s file extensions with a different file extension.',
        epilog='This script is part of the CliMush bioinformatics pipeline.',
    )

    # add required arg for input file path
    parser.add_argument(
        '-i', '--input',
        type=pathlib.PosixPath,
        required=True,
        help='Path to a directory containing the files that require the file extension to be replaced.',
    )

    # add required arg for replacement file extension to use
    parser.add_argument(
        '--rep',
        type=str,
        help='The file extension to use when replacing the input file\'s file extension.',
    )

    # parser arguments into a dictionary
    args = vars(parser.parse_args())


## FORMAT THE --REP STRING TO FILE EXTENSION ##

# check it is formatted correctly, need leading '.'
if args['rep'].startswith('.'):
    pass

# if it doesn't have a leading '.', add one here
else:
    args['rep'] = '.' + args['rep']


## ITERATE THROUGH INPUT FILES TO CREATE REPLACEMENT FILE NAME (WITH --REP FILE EXTENSION)

# store input / output file names in dictionary
rename_dict = {}

for input_seqfile in args['input'].glob('*.fast*'):

    # get the input file's basename without its file extension
    # use the suffixes and not just suffix to ensure all file extensions removed if there are multiple
    input_basename = input_seqfile.name.replace(''.join(input_seqfile.suffixes), '')

    # add the --rep file extension to the input_basename to create the output basename
    output_seqfile = (input_seqfile.parent / input_basename).with_suffix(args['rep'])

    # update the dictionary with the input / output file paths
    rename_dict.update({input_seqfile: output_seqfile})


## RENAME FILES TO REPLACE FILE EXTENSION
for input_seqfile, output_seqfile in rename_dict.items():

    # rename the input files with the output file extension
    input_seqfile.replace(output_seqfile)

