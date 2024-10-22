# list of required packages
required_packages <- c('tidyverse', 'BiocManager')

# if any of the required packages are not installed, install them
for (pkg in required_packages){
  if (!pkg %in% c(installed.packages()[,1])){
  install.packages(pkg, dependencies = TRUE)
  }
}

# then load required packages into library
library(required_packages)

# then have BiocManage install 'biomformat', if not already installed
if (!'biomformat' %in% c(installed.packages()[,1])){
  BiocManager::install("biomformat")
}
