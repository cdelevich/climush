from pathlib import Path
# BE CAUTIOUS OF USING RESOLVE/ABSOLUTE -- IT OFTEN WILL CHANGE THINGS SO THAT IT IS DEPENDENT
# ON WHERE YOU RUN THE SCRIPT FROM!!!
# TRY RELATIVE_TO?

OUTPUT_DIRS = {'demultiplex':'demultiplexed',
               'prefilter':'prefiltered',
               'remove-primers':'primers-removed',
               'dereplicate-01':'dereplicated-01',
               'quality-filter':'quality-filtered',
               'separate-subregions':'subregions-separated',
               'chimera-check':'chimera-checked',
               'dereplicate-02':'dereplicated-02',
               'cluster-otus':'otus-clustered',
               'create-otu-table':'otu-table',
               'assign-taxonomy':'taxonomy-assigned'}

RELATIVE_PATHS = {'pipeline-output':OUTPUT_DIRS,
                  'climush_py-package':{},
                  'pipeline':{},
                  'sequences':{'needs_rename':'',
                               'needs_demux':'',
                               'needs_clarification':'',
                               'pacbio':'',
                               'illumina':'',
                               'sanger':''},
                  'config':{'barcode_mapping':{}}}

# PIPELINE OUTPUT PATHS
# PIPE_OUT_MAIN = Path('pipeline-output').resolve()
# remove these, create within script by pulling name of script
# PIPE_OUT_DEMUX = Path(PIPE_OUT_MAIN / 'demultiplexing')
# PIPE_OUT_PREFILTER = Path(PIPE_OUT_MAIN / 'prefilter')
# PIPE_OUT_PRIMTRIM = Path(PIPE_OUT_MAIN / 'remove-primers')

# DOCKER CONTAINER PATHS RELATIVE TO PIPELINE SCRIPTS
# PIPELINE_PATH = Path('pipeline')
# SEQUENCES_PATH = Path('sequences')
# CONFIG_PATH = Path('config')
# MAPPING_PATH = Path(CONFIG_PATH / 'barcode-mapping')

# MOCK COMMUNITY PATH TO FILE
# MOCK_PATH = Path(CONFIG_PATH / 'illumina_mock-community')

# FILE SORTING NEW DIRECTORIES
# DEMUX_DIR = 'needs_demux'
# RENAME_DIR = 'needs_rename'
# UNCLEAR_DIR = 'needs_clarification'

# NO FILE PATHS BELOW


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