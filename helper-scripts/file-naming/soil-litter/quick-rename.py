import argparse, tomlkit, itertools, os, shutil, re, sys
from pathlib import Path
import pandas as pd
from datetime import datetime
from climush.utilities import prompt_yes_no_quit, is_pathclass, import_filepath, flag_multiple_files, get_settings, mkdir_exist_ok
from climush.constants import CONFIG_FILETYPE, SEQ_FILE_RE, GZIP_REGEX, GZIP_GLOB, MOCK_COMM_RE, NEG_CTRL_RE, UNDET_RE

# guide names
rename_guide = pd.read_csv('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/illumina_soil-litter_2023-05_file-rename-conversion.csv').to_dict(orient='dict')

# files to rename
rename_dir = Path('/Users/carolyndelevich/main/github_repos/climush/helper-scripts/file-naming/illumina_soil-litter_2023-10_raw-reads_empty')
rename_file_paths = list(rename_dir.glob('*.fastq.gz'))
rename_file_names = [f.stem.replace('.fastq', '') for f in rename_file_paths]


# get new name by aligning guide with old names
for file in rename_file_names:
    sample_id = '_'.join(file.split('_')[:5])
    for guide in guide_old:
        if guide.startswith(sample_id):