from pathlib import Path

NEEDS_ACTION_REGEX = r'needs_(?!demux)'

# PACBIO_DIR = 'pacbio'
# ILLUMINA_DIR = 'illumina'
# SANGER_DIR = 'sanger'

# PLATFORM LIST
SEQ_PLATFORM_OPTS = ['pacbio', 'illumina', 'sanger']

# CONFIGURATION FILE HANDLES
CONFIG_FILETYPE = '.toml'
RENAME_CONFIG_HANDLE = 'file-rename'
PIPELINE_CONFIG_HANDLE = 'pipeline-settings'
NAME_CONFIG_HANDLE = 'filename-components'

# CONFIGURATION FILE SECTIONS
PIPELINE_CONFIG_MPID = 'pacbio-multiplex-ids'

# MISCELL REGEX
AFFIRM_REGEX = '^c|(continue)|^y'

YES_RE = '^yes$|^y$'
NO_RE = '^no$|^n$'
QUIT_RE = '^quit.*?|^exit.*?'  # avoids confusion about whether to use quit or quit() or exit or exit()

HIDDEN_FILE_REGEX = '^\.'
GZIP_REGEX = '(\.fast.\.gz)$'
ANY_PLATFORM_REGEX = '^illumina|^pacbio|^sanger'
SEQ_FILE_GLOB = '*.fast*'

# PIPELINE FILE PREFIXES
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