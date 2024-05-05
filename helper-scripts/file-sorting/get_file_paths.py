import pathlib, argparse
from pathlib import Path

# define command line options
parser = argparse.ArgumentParser(prog='GetFilePaths',
                                 description='Return a file of all directories in the input directory, including '
                                             'all subdirectories and file contents.')

# required; specify the path to the directory to be sorted
parser.add_argument('input_path',
                    help='Path to the main directory for the sequencing files that require sorting.')

# optional; otherwise, will write to the input file path
parser.add_argument('-o', '--output',
                    default=None,
                    help='Optional path to the output file containing the file paths. If none is provided, will '
                         'create an output file in the same location as the input path.')

args = parser.parse_args()

search_path = Path(args.input_path)

if args.output is None:
    output_path = (search_path / f'{search_path.name}_filepaths').with_suffix('.txt')
else:
    output_path = Path(args.output)

# list of items returned for each directory from Path.walk()
subitem_list = ['directory path', 'subdirectories', 'files']

# write out to output_path
with open(output_path, 'w+') as fout:
    for root in search_path.walk():  # go through each directory found by Path.walk()
        for i,items in enumerate(root):  # then access the root path, subdirectories, and filenames
            if i == 0:  # don't indent if it is the first, which is the root path
                fout.write('---------------------------------------------\n')  # for each new root, insert a separator
                fout.write(f'{subitem_list[i]}: {items}\n')
                fout.write('---------------------------------------------\n')  # for each new root, insert a separator
            else:
                if isinstance(items, list):  # format with one print out per line if its a file path
                    fout.write(f'\n  {subitem_list[i]}')  # print the item description
                    if len(items) > 0:
                        items.insert(0, '\t')  # add tab so there's one before first item
                        formatted_items = '\n\t'.join(items)  # format the strings in the item to print one line at a time
                        fout.write(f'{formatted_items}\n')  # then the formatted items
                    else:
                        fout.write(f'\n\t(none)\n')  # write out 'none' for empty list (i.e., no items)
                else:  # if the item is a string and not a list, just write the string
                    fout.write(f'\n  {subitem_list[i]}')
                    fout.write(f'\t{items}\n')
        fout.write('\n\n\n')  # write extra space between directories
