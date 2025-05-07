from climush.utilities import compress_data
from climush.constants import GZIP_GLOB
from pathlib import Path
import argparse, pathlib, warnings, re


parser = argparse.ArgumentParser(prog=Path(__file__).stem,
                                 description='Locate files that did not get compressed and compress those files.',
                                 epilog='This script is part of the CliMush bioinformatics pipeline.')

parser.add_argument(
    '-i', '--input',
    type=pathlib.PosixPath,
    required=True,
    help='Path to a directory that contains a mix of compressed and uncompressed sequence files that need to all be '
         'converted to compressed format.'
)

# parser.add_argument(
#     '--decompress',
#     action='store_true',
#     help='Flag that, when used, will search for the compressed files and decompress them, which is opposite the '
#          'normal behavior of this script in which the uncompressed files are replaced by their compressed format.'
# )

args = vars(parser.parse_args())


# keep track of the number of compressed files, uncompressed files, and final compressed-only files
uncompressed_filecount_initial = 0

# go through each of the input files that are uncompressed
for uncompressed_file in args['input'].glob('*.fastq'):

    # create a Boolean flag that will switch to True if there's a compressed counterpart file
    compressed_found = False

    # add to the starting uncompressed file counter
    uncompressed_filecount_initial += 1

    # initiate the compressed file counter here, it loops below so I want it to restart each time, keep last
    compressed_filecount_initial = 0

    # get the file's basename (i.e., with suffixes removed)
    uncompressed_file_basename = uncompressed_file.name.replace(''.join(uncompressed_file.suffixes),'')

    # search for a corresponding compressed file with this same basename
    for compressed_file in args['input'].glob(GZIP_GLOB):

        # add to compressed file count initiall counter
        compressed_filecount_initial += 1

        # search for basename in the compressed file names
        compressed_file_match = re.search(uncompressed_file_basename, compressed_file.name, re.I)

        # if there's a compressed file that also has the uncompressed file's basename...
        if compressed_file_match:

            # switch the Boolean to True
            compressed_found = True

            # don't keep looping, break here
            break

        # if this compressed file isn't the one matchign the uncompressed file, keep searching
        else:
            continue


    # if the compressed file was located...
    if compressed_found:

        # remove the uncompressed file
        uncompressed_file.unlink()

    # if the compressed file was not located...
    else:

        # use the compress_data() function on this file to compress it, removing the uncompressed version
        compress_data(
            input_path=uncompressed_file,
            output_path=None,
            compress_fmt='gzip',
            keep_input=False,
        )


# count the final number of compressed files, compare to initial count
compressed_filecount_final = len([ comp_f for comp_f in args['input'].glob(GZIP_GLOB) ])
uncompressed_filecount_final = len([ uncomp_f for uncomp_f in args['input'].glob('*.fastq') ])

print(f'There are {compressed_filecount_final} gzip-compressed files in the input directory, compared to '
      f'an initial {compressed_filecount_initial} gzip-compressed files.\n')
print(f'There are {uncompressed_filecount_final} uncompressed files in the input directory, compared to '
      f'an initial {uncompressed_filecount_initial} uncompressed files.\n')