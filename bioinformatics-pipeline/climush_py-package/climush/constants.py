
NEEDS_ACTION_REGEX = r'needs_(?!demux)'

# PACBIO_DIR = 'pacbio'
# ILLUMINA_DIR = 'illumina'
# SANGER_DIR = 'sanger'

# PLATFORM LIST
SEQ_PLATFORM_OPTS = ['pacbio', 'illumina', 'sanger']  # no longer used in get_seq_platform()

# CONFIGURATION FILE HANDLES
CONFIG_FILETYPE = '.toml'
RENAME_CONFIG_HANDLE = 'file-rename'
PIPELINE_CONFIG_HANDLE = 'pipeline-settings'
NAME_CONFIG_HANDLE = 'filename-components'

# CONFIGURATION FILE SECTIONS
PIPELINE_CONFIG_MPID = 'pacbio-multiplex-ids'

# SEQUENCE FILE NAMES
MOCK_COMM_RE = r'^mock'
NEG_CTRL_RE = r'^ntc'
UNDET_RE = r'^undet'

# MISCELL REGEX
AFFIRM_REGEX = r'^c|(continue)|^y'

YES_RE = r'^yes$|^y$'
NO_RE = r'^no$|^n$'
QUIT_RE = r'^quit.*?|^exit.*?'  # avoids confusion about whether to use quit or quit() or exit or exit()

HIDDEN_FILE_REGEX = r'^\.'
GZIP_REGEX = r'(\.fast.\.gz)$'
GZIP_GLOB = '*.fastq.gz'
ANY_PLATFORM_REGEX = r'^illumina|^pacbio|^sanger'
PLATFORM_ANYWHERE_RE = r'illumina|pacbio|sanger'
SEQ_FILE_GLOB = '*.fast*'
SEQ_FILE_RE = r'\.fast.$'


ORIENT_RE = r'(?<=_0[1-9]_)R[1,2](?=\.)'

# DETECT COLUMN NAMES IN BARCODE MAPPING TABLES #######
# used in demultiplex() in bioinfo.py to detect columns
POOL_NUM_RE = r'(?<=pool).?(\d)'  # detect multiplexed pool number in barcode table tabs
FWD_COL_RE = r'(fwd)|(forward)'  # find columns with forward barcodes
REV_COL_RE = r'(rev)|(reverse)'  # find columns with reverse barcodes
SAMPLE_COL_RE = r'^sample'  # find columns that contain the sample ID


# to use as default in get_sample_id() if platform not found
SAMPLE_ID_RE = r'(?<=\w_)(pacbio|sanger|illumina).+?(?=\.)'


#######################################################################################################################
# FOR FASTA HEADERS THAT **HAVE NOT** BEEN RENAMED BY RENAMED_READ_HEADER() FROM UTILITIES.PY #########################
# these regex only apply if the read headers have not been renamedj
# currently, these regex also only apply if the reads have been run through ITSx and have the ITSx headers

# used to locate a simplified read ID from the long read ID output in PacBio read headers
# this read ID is used as a unique identifier for that read, and used to create a more detailed sample ID within a fasta
# example of detailed sample ID = pacbio_sporocarp-f_2022-04_154_748596 (where 748596 is the read ID)
READ_ID_OG_RE = r'[0-9]{1,}(?=\/ccs)'  # used in rename_read_header() in utilities.py, check_itsx_output() in bioinfo.py

# get the gene region of a read from the ITSx-formatted read header
READ_REGION_OG_RE = r'(?<=\|\w\|).+(?=\sExtracted)'  # used in rename_read_header() in utilities.py

# get the read length of a read from the ITSx-formatted read header
READ_LEN_OG_RE = r'[0-9]{1,}(?=\sbp)'  # used in rename_read_header() in utilities.py

# get the read count of a read from the ITSx-formatted read header
READ_COUNT_OG_RE = r'(?<=size=)[0-9]{1,}'  # used in rename_read_header() in utilities.py


#######################################################################################################################
# FOR FASTA HEADERS THAT **HAVE** BEEN RENAMED BY RENAMED_READ_HEADER() FROM UTILITIES.PY #############################
# these regex can only be used after the function rename_read_header() from utilities.py has renamed the read header

# for concatenating reads after ITSx
PREFIX_RE = r'^\w.+?(?=_pacbio|_sanger|_illumina)'  # used in search_path_with_regex() in concat_regions()

# used to locate the updated sample ID within a read header
# may still return a substring even if the header has not been renamed, but will not be the intended substring
READ_ID_RENAMED_RE = r'^.+?(?=;)'  # used in concat_regions() in bioinfo.py

# get the gene region of a read after the header has been reformatted by rename_read_header() from utilities.py
READ_REGION_RENAMED_RE = r'(?<=[A-Z]\|)\w.+?(?= Extracted)'  # used in concat_regions() in bioinfo.py

# get the read length of a read after the header has been reformatted by rename_read_header() from utilities.py
READ_LEN_RENAMED_RE = r'(?<=\()[0-9]{1,}(?= bp\))'  # used in concat_regions() in bioinfo.py

# get the read count of a read after the header has been reformatted by rename_read_header() from utilities.py
READ_COUNT_RENAMED_RE = r'(?<=;size=)[0-9]{1,}(?=\|)'  # used in concat_regions() in bioinfo.py


# PIPELINE FILE PREFIXES
DEMUX_PREFIX = 'demux'
NOPHIX_PREFIX = 'no-phix'
NOAMBIG_PREFIX = 'no-ambig'
TRIMMED_PREFIX = 'trim'
DEREP_PREFIX = 'derep'
MERGED_PREFIX = 'merged'
QUALFILT_PREFIX = 'qualfilt'
ITSX_PREFIX = 'itsx'
NOCHIM_PREFIX = 'no-chim'
CLUSTER_PREFIX = 'clust'
REPREAD_PREFIX = 'rep-read'

# PIPELINE FILE SUFFIXES
LOG_SUFFIX = '.log'