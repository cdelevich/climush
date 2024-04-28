from mapping import filepath_map as fpm

import argparse, sys, subprocess
from pathlib import Path
from climush.constants import *
# from climush.bioinfo import create_blast_db, assign_taxonomy
from climush.utilities import *

settings = get_settings(fpm)