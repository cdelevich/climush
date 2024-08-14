import pathlib, argparse, subprocess, re, sys, json, shutil
from pathlib import Path
from climush.utilities import check_dir_exists, is_pathclass

# define command line options
parser = argparse.ArgumentParser(prog='SortPacbioFiles',
                                 description='Sorts the PacBio sequencing file output from the University of '
                                             'Oregon\'s sequencing core (GC3F), so that input files are more '
                                             'easily distinguished from extraneous information from the output.')

# required; specify the path to the directory to be sorted
parser.add_argument('input_path',
                    help='Path to the main directory for the sequencing files that require sorting.')

# include if you want to allow previous directories to be overwritten (will not clobber, pretty sure...)
parser.add_argument('--overwrite', action='store_true',
                    help='If flag is used, the new directories created will overwrite previously existing '
                         'directories with the same name.')

# prevent a parent directory from being created, if non-existent parent included in a new directory
parser.add_argument('--no-parent', action='store_false',
                    help='Turn off the -p flag, which has the same meaning as in mkdir from the command line. '
                         'If a new directory path to be created has a parent that does not yet exist, this flag '
                         'prevents the parent (and therefore child) from being created.')

# for testing, will create a copy of everything in the input directory, but with empty files
# it will then perform the file sorting on this empty version to ensure everything works without risk
parser.add_argument('--dry-run', action='store_true',
                    help='Perform a dry run of the file sorting, where an exact copy of the original directory '
                         'tree')

# for testing, will keep the empty directory created when --dry-run is used, which is normally removed after the
# dry run successfully completed
parser.add_argument('--keep-run', action='store_true',
                    help='Keep the empty directories and files that were produced by --dry-run.')

# if files have been sorted using this script, will restore the files/directories to their original structure
parser.add_argument('--restore', action='store_true',
                    help='Restore the files and directories to their original file structure before sorting.')

args = parser.parse_args()


#############################################
# GET ALL DIRECTORIES + FILES IN RUN DIR ####
#############################################

# create a Path object from the input path string from args; confirm that this Path exists
input_path = check_dir_exists(dir_path=args.input_path)

# get the name of the sequencing run from the input_path directory name
run_id = input_path.name

# if the dry run flag is used, create a dummy copy of the input directory to use as the input
if args.dry_run:

    # run copy w/ the --attributes-only flag to create empty copy of all files/subdirs (with -r command)
    dummydir_path = input_path.parent / f'{run_id}_dry-run' # put in same location as the input path's parent
    dummydir_cmd = ['cp', '-r', '--attributes-only', input_path, dummydir_path]
    subprocess.run(dummydir_cmd)

    # then replace the input_path with the path to the dummy copy of the input_path
    input_path = dummydir_path

#####################################################
# RESTORE ORIGINAL FILE STRUCTURE (IF REQUESTED) ####
#####################################################

# define path to the summary directory, which is used later but required now for a restore
summary_main = input_path / f'{run_id}_summary'

# define the path to the .json file sorting output here; saved at end but need path now for restore
file_rename_path = (summary_main / f'{run_id}_file-sorting').with_suffix('.json')

if args.restore:

    if file_rename_path.is_file():

        # find the file sorting .json file that is output when this script is run
        with open(file_rename_path) as json_in:
            restore_dict = json.load(json_in)

        # json package does not allow PosixPath as dict key, so this dict only contains strings
        for src, trg in restore_dict.items():
            Path(trg).rename(Path(src))  # need to reset string paths to PosixPaths

        print(f'The original file structure was successfully restored.\n')

    else:

        print(f'This directory:\n'
              f'  {input_path}\n'
              f'does not appear to have been renamed using the script {Path(__file__).name}. It attempted to locate '
              f'the file:\n'
              f'  {file_rename_path}\n'
              f'to restore the original file structure but could not locate it.\n')

    sys.exit()

#############################################
# SAVE COPY OF ORIGINAL FILE TREE ###########
#############################################

# do not do this for a dry run
if args.dry_run:
    pass
else:
    # original file tree path
    original_tree_path = (summary_main / f'{run_id}_original-file-tree').with_suffix('.txt')

    def save_file_tree(parent_dir, out_path, save=True):
        '''
        Create a text file containing the file tree of the input parent directory.

        :param parent_dir: top level directory within the tree
        :param out_path: file path to save the tree to
        :return: None; saves file in background
        '''

        # assemble Unix command to get tree
        tree_cmd = ['tree', parent_dir]

        # execute command
        tree_out = subprocess.run(tree_cmd, capture_output=True, text=True)

        if save:
            # write stdout to file
            with open(out_path, 'w') as fout:
                fout.write(tree_out.stdout)

            return None
        else:
            return tree_out

    # do not save now, since output path does not yet exist
    # cannot just use function later, because don't want new directory (summary_main) in original tree
    og_tree = save_file_tree(parent_dir=input_path,
                             out_path=original_tree_path,
                             save=False)

#############################################
# CREATE NEW DIRECTORIES FOR SORTING ########
#############################################

# QUALITY SCORE FOLDERS

# quality score folder dictionary
qscore_paths = {k:{'main': '',  # root dir of each qscore folder
                   'ccs': '',   # circular consensus reads for each qscore
                   'pool-demux': {'main': '',
                                  'lima-info': ''}}
                for k in ['Q20', 'Q30', 'Q40']}

# create the subdirectories for each quality score directory
for q in qscore_paths:

    # create path to main quality score dir
    q_path_main = input_path / f'{run_id}_{q}'

    # move through the subdirectories to create these directories
    for p in qscore_paths[q]:

        # if it's the main/root val for this path, add to the dict
        if p == 'main':
            qscore_paths[q].update({p: q_path_main})

        # this is the only subdir with its own subdirs
        elif p == 'pool-demux':
            # create main directory for this sub, 'main' tag
            qsubsub_path = q_path_main / f'{run_id}_{q}_{p}'

            # for through each subdir
            for s in qscore_paths[q][p]:

                # if its the 'main' tag, just assign the main dir for this sub
                if s == 'main':
                    qscore_paths[q][p].update({s: qsubsub_path})

                # otherwise, create path based on this sub's main, and mkdir (needs parents=True)
                else:
                    qss_path = qsubsub_path / f'{qsubsub_path.name}_{s}'
                    qscore_paths[q][p].update({s: qss_path})

                    qss_path.mkdir(mode=0o777, parents=args.no_parent, exist_ok=args.overwrite)

        # otherwise, create a path, then directory under this main parent dir
        else:

            # create a path for this directory using the key value
            qsub_path = q_path_main / f'{q_path_main.name}_{p}'

            # use this qscore root path to create the subdirectories
            qsub_path.mkdir(mode=0o777, parents=args.no_parent, exist_ok=args.overwrite)

            # confirm the path was created successfully
            if qsub_path.is_dir():  # check if the directory exists (i.e., was properly created)
                qscore_paths[q].update({p: qsub_path})  # if it was, then add to path dict
            else:  # if it was not, then print error and exit
                print(f'ERROR. The quality score folder was not created:\n'
                      f'\t{qsub_path}\n'
                      f'Exiting {__file__.name}...')
                sys.exit()


# SUBREAD INFORMATION

# subread files, prior to create consensus reads (I believe)
subread_path = input_path / f'{run_id}_subreads'
subread_path.mkdir(mode=0o777, parents=args.no_parent, exist_ok=args.overwrite)

# READ INFORMATION

# I think the reads are after creating consensus reads, but before separating into quality scores
read_path = input_path / f'{run_id}_reads'
read_path.mkdir(mode=0o777, parents=args.no_parent, exist_ok=args.overwrite)


# SEQUENCE RUN SUMMARY INFORMATION

summary_paths = {'main': '',
                 'file-transfer': ''}

for p in summary_paths:

    # add the main directory path to the dictionary
    if p == 'main':
        summary_paths[p] = summary_main

    # if it isn't main, then mkdir and add to dict for this dir
    else:
        subdir = summary_main / f'{run_id}_{p}'
        subdir.mkdir(mode=0o777, parents=args.no_parent, exist_ok=args.overwrite)
        summary_paths[p] = subdir

if not args.dry_run:
    with open(original_tree_path, 'w') as fout:
        fout.write(og_tree.stdout)

#########################################################
# WALK UNSORTED PATH AND RECORD FILE LOCATIONS ##########
#########################################################

# walk through all directories in the input_path
# do not sort files during the walk, or it will make a mess, I think
# instead, create a dictionary, where the key is the source (old path) and the value is the target (new path)
# this is also a good way to keep track of what changes were made, as this dict can be written out as .json log

# dictionary to store the old file paths (key) and new file paths (values)
src_trg_dict = {}
# this dictionary will be used to rename (i.e., move) the files after the walk completes

# keep master list of all files, then make sure they are accounted for in srg_trg_dict
full_file_list = []

# for each directory and subdirectory in the input path, inspect the directory's path, subdirectories, and files
for root, subdir, files in input_path.walk():

    # sort the files in this main directory
    for file in files:

        # the source path will be the same set-up for all files
        src = root / file

        # add all files into full_file_list, to ensure all are transferred to src_trg_dict
        full_file_list.append(src)

        # if this is the main root directory...
        if root == input_path:

            # any files I made to summary directory
            if re.search(f'^{run_id}', file, re.I):
                trg = summary_paths['main'] / file

            # if any of the files contain a subread tag, move to subreads directory
            # unsure what .xml file is, looks like some read info (pre-consensus)
            # not totally sure what the .log file is either, so put here
            # barcodes.fa is for the core's demux process
            elif (re.search(r'.+?subread.+?', file, re.I)
                  or src.suffix == '.xml'
                  or src.suffix == '.log'
                  or src.suffix == '.fa'):
                trg = subread_path / file

            # if the word 'transfer' is in the file, put in file-transfer subdir
            elif re.search(r'.+?transfer.+?|^tmp', file, re.I):
                trg = summary_paths['file-transfer'] / file

            # if this is the file the describes the directory contents according to GC3F
            elif re.search(r'^info.+?\.txt', file, re.I):
                trg = summary_paths['main'] / f'pacbio_gc3f-file-descriptions.txt'

            else:
                trg = summary_paths['main'] / file

        # if this path contains ccs.Q20, ccs.Q30, or ccs.Q40 anywhere in the path...
        elif re.search(r'ccs.Q\d{2}', str(root), re.I):

            qscore = re.search(r'(?<=ccs\.)Q\d{2}', str(root), re.I).group(0)

            # if the path name contains ccs.Q20, ccs.Q30, or ccs.Q40...
            if re.search(r'ccs.Q\d{2}', root.name, re.I):

                # update the filename so that it includes the quality score
                file_parts_all = file.split('.')
                file_name = f'{file_parts_all[0]}_{qscore}'
                file_ext = '.'.join(file_parts_all[1:])
                new_file = f'{file_name}.{file_ext}'
                trg = qscore_paths[qscore]['ccs'] / file

            # if the directory is a Q## subdirectory, demultiplexed...
            # this is GC3F's demultiplexing, not mine (pool demux, not sample demux)
            elif re.search(r'demultiplex.+?', root.name, re.I):

                # if it is a reads file separated into a pool...
                if re.search(r'(?<=^\d{4}\.)P[ol]{0,3}\d--P[ol]{0,3}\d(?!\.consensusreadset\.)',
                             file, re.I):

                    # update the filename so that it includes the quality score
                    file_parts = file.split('.')
                    file_parts.insert(2, qscore)
                    new_file = '.'.join(file_parts)
                    trg = qscore_paths[qscore]['pool-demux']['main'] / new_file

                # if not one of these demux sequence files, then put into lima-info subdir
                else:
                    trg = qscore_paths[qscore]['pool-demux']['lima-info'] / file

            # I'm not sure what files this would encompass, but I'll send to main Q## directory
            else:
                trg = qscore_paths[qscore]['main'] / file

        # if this is the /reads directory within main
        elif re.search(r'read.?', str(root), re.I):  # flexible to read or reads as directory name
            trg = read_path / file

        # the INFO directory contains info related to the file transfer from GC3F to Globus shared dir
        elif re.search(r'info', str(root), re.I):
            trg = summary_paths['file-transfer'] / file

        else:
            trg = summary_paths['main'] / file

        # add the old and new paths to the dictionary
        src_trg_dict[src] = trg


##################################################
# MOVE THE FILES TO THEIR UPDATE LOCATION ########
##################################################


# CONFIRM ALL FILES ARE ACCOUNTED FOR IN DICTIONARY

# create a set of the dict and list, in order to find difference between the two
dict_set = set(src_trg_dict.keys())
list_set = set(full_file_list)
missing_files = list(list_set.difference(dict_set))

# check that the difference is 0 before continuing
# if the difference is not zero, print the missing files and exit process before anything is moved
if len(missing_files) > 0:
    missing_files.insert(0, '\n\t')
    missing_files.sort()
    formatted_missing = '\n\t'.join(missing_files)
    print(f'The following {len(missing_files)} files are missing from the file shuffling dictionary. Please review '
          f'this list, and make any changes to the {Path(__file__).name} script, then rerun.')
    print(formatted_missing)
    print('\nNo file sorting has occurred. Exiting...')
    sys.exit()


# USE DICTIONARY TO MOVE FILES

# go through the dict, and sort the files

# json package does not allow PosixPath objects as keys, so create new dict to add paths to as strings
src_trg_str_dict = {}

for src, trg in src_trg_dict.items():

    # rename the original file with the updated file path (i.e., 'move')
    src.rename(trg)

    # add source/target paths as strings to .json output dict
    src_trg_str_dict[str(src)] = str(trg)  # cannot be PosixPath for .json output

def find_and_delete_empty_dirs(parent_dir):
    '''
    Locate and remove empty directories within a parent directory.

    :param parent_dir: the top-level directory in which to search and remove empty
    directories
    :return: none
    '''

    if args.dry_run:

        shutil.rmtree(parent_dir)
        print(f'The dry run of file sorting was successful.')
        return sys.exit()

    else:
        empty_dirs = []

        for root, subdir, _ in parent_dir.walk(top_down=False):

            for sd in subdir:
                full_path = root / sd
                try:
                    full_path.rmdir()  # should only remove if the directory is empty
                    empty_dirs.append(('  ' + str(full_path) + '\n'))  # format for printing before adding to list
                except OSError:  # should raise OSError if directory if not empty
                    continue

        print(f'The following {len(empty_dirs)} directories were empty and therefore removed:\n')
        print(*empty_dirs)

        return None

# REMOVE ALL EMPTY DIRECTORIES
# if this is a dry run, it will remove all of the directories it produced
find_and_delete_empty_dirs(parent_dir=input_path)

# SAVE THIS DICTIONARY FOR RECORD OF OLD FILE LOCATIONS AND NEW FILE LOCATIONS

# export this dictionary, as record of old and new file paths
with open(file_rename_path, 'w+') as jout:
    json.dump(src_trg_str_dict, jout)

#############################################
# SAVE COPY OF NEW FILE TREE ################
#############################################

# original file tree path
new_tree_path = (summary_main / f'{run_id}_updated-file-tree').with_suffix('.txt')

save_file_tree(parent_dir=input_path,
               out_path=new_tree_path)
