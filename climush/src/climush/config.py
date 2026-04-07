import warnings
import pandas as pd
import matplotlib.pyplot as plt

# set user preferences for warnings, figure quality, etc.
pd.options.mode.chained_assignment = None  # default='warn', blocks copy on df warning
warnings.simplefilter(action='ignore', category=FutureWarning)  # supress future warnings
pd.set_option('display.max_columns', 3000)
plt.rcParams['savefig.dpi'] = 300