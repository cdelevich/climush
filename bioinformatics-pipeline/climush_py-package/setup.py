from setuptools import setup

setup(name='climush',
      version='0.1',
      description='Custom Python functions for bioinformatics for the NSF Macrosystems CliMush project',
      author='Carolyn Delevich',
      author_email='cdelevic@uoregon.edu',
      packages=['climush'],
      install_requires=['pandas', 'datetime', 'scipy', 'numpy', 'matplotlib',
                        'seaborn', 'argparse', 'Bio'])  # Python modules excluded (e.g., os, sys)