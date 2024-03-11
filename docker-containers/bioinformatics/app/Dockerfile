# syntax=docker/dockerfile:1

# Copyright (c) Carolyn Delevich
# Distributed under the terms of the Modified BSD License.

# latest version of Python to match Conda install version
FROM python:3.12-slim-bookworm

LABEL maintainer="Carolyn Delevich <cdelevic@uoregon.edu>"

# The base systems is a Debian release configured with Python 3.12.  Add
# a Java runtime environment and the bioinformatics and Python apps.

# import variables defined in the Docker compose file
ARG VERSION_NUM

# install base packages then clear the apt cache
RUN apt update --yes && \
    apt install --yes --no-install-recommends \
    default-jre \
    make \
    nano  \
    micro \
    unzip \
    wget && \
    rm -rf /var/lib/apt/lists/*

# install the climush Python package, which contains the functions required for the pipeline
#RUN pip install climush

# install conda and give permission to write to the conda folder (for installing other packages?)
ENV CONDA_DIR /opt/conda
RUN wget --quiet -O ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ~/miniconda.sh -b -u -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH

# We don't want the container to run as root, so make a user account
# (the -m option makes the user's home directory as well).  
ENV USER=mycophile
RUN useradd -m $USER

COPY mycoformatics.yml /tmp/mycoformatics.yml
RUN conda env create -f /tmp/mycoformatics.yml

RUN echo "source activate mycoformatics" > ~/.bashrc
ENV PATH $CONDA_DIR/envs/mycoformatics/bin:$PATH

# Add the home folder to PATH:
ENV PATH=/home/$USER:$PATH

# mount volumes for pipeline-specific directories, copy content over
WORKDIR /home/$USER
VOLUME ["/home/$USER/climush_py-package", "/home/$USER/config", "/home/$USER/pipeline", "/home/$USER/sequences"]
COPY --chown=$USER:users climush_py-package /home/$USER/climush_py-package
COPY --chown=$USER:users config /home/$USER/config
COPY --chown=$USER:users pipeline /home/$USER/pipeline
COPY --chown=$USER:users bash_aliases /home/$USER/.bashrc

# Define an environment var to indicate script is running in a container
ENV CONTAINERIZED=''
USER $USER

# set bash as entrypoint, start by sourcing the .bashrc file
CMD ["source", "/home/$USER/.bashrc"]
ENTRYPOINT ["/bin/bash"]