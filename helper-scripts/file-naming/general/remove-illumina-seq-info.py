from pathlib import Path
import pathlib, argparse, re
import pandas as pd

TEST_MODE=False

if TEST_MODE:

    args = {
        'input': Path('/helper-scripts/file-naming/endophytes/test-files_2023-endos_prefix-replaced'),
        'log': True,
    }

else:

    # instantiate parser
    parser = argparse.ArgumentParser(
        prog=Path(__file__).stem,
        description='Remove illumina sequencer information from file names',
        epilog='This script is part of the CliMush bioinformatics pipeline.',
    )

    # add required input path parameter
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=pathlib.PosixPath,
        help='path to sequence files that need illumina sequencer info removed.',
    )

    # optional flag to keep a original to new file name conversion table
    parser.add_argument(
        '--log',
        action='store_true',
        help='flag that, when used, will write out a table with the original file names and '
             'their associated updated file name.',
    )

    # parse arguments into a dictionary
    args = vars(parser.parse_args())


# create a dictionary of the input sequence files (old file names) to add updated file names to
rename_dict = {
    original_filepath:'' for original_filepath in args['input'].glob('*.fast*')
}

for original_filepath in rename_dict:

    # get just the file name from the file path
    original_filename = original_filepath.name

    # search for the S## string in the original file name
    updated_filename01 = re.sub(r'_S\d{1,3}(?=_)','',original_filename)

    # search for the trailing 001 in the original file name
    updated_filename02 = re.sub(r'_001(?=\.fast)','', updated_filename01)

    # create the filepath for the replacement file name
    updated_filepath = original_filepath.parent / updated_filename02

    # add the updated file path to the dictionary, paired with the original file path
    rename_dict.update({original_filepath:updated_filepath})

# go through the renaming dictionary and rename the original filenames iwth the updated filenames
for original_filepath, updated_filepath in rename_dict.items():

    # replace original with updated
    original_filepath.replace(updated_filepath)

# if --log flag used, write out conversion to .csv file
if args['log']:

    # create output file log name
    output_log_path = args['input'].parent / f'remove-illumina-seq-info_{args["input"].name}.csv'

    # get only the file names for the output table
    original_filenames = [file_og.name for file_og in rename_dict.keys()]
    updated_filenames = [file_up.name for file_up in rename_dict.values()]

    # create dataframe from dictionary
    output_dict_df = pd.DataFrame(
        data=zip(original_filenames, updated_filenames),
        columns=['original_filenames', 'updated_filenames']
    )

    # output to .csv file
    output_dict_df.to_csv(output_log_path, index=False)

# otherwise, do not write out this file
else:
    pass