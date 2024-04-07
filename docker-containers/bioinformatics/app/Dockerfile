# syntax=docker/dockerfile:1

# Copyright (c) Carolyn Delevich
# Distributed under the terms of the Open Source Initiative MIT License.

# latest version of Python to match Conda install version
FROM python:3.12-slim-bookworm

LABEL maintainer="Carolyn Delevich <cdelevic@uoregon.edu>"

# The base systems is a Debian release configured with Python 3.12.

# version number of this image from the Docker compose file
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

# install the climush Python package, which contains custom functions required for the pipeline
RUN pip install climush

# REMOVE AFTER TESTING
RUN echo "just a test" > test.txt

# install conda and give permission to write to the conda folder (for installing other packages?)
ENV CONDA_DIR /opt/conda
RUN wget --quiet -O ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ~/miniconda.sh -b -u -p $CONDA_DIR && \
    rm ~/miniconda.sh

# add conda executables to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# create an environment from the conda mycoformatics environment
COPY mycoformatics.yml /tmp/mycoformatics.yml
RUN conda env create -f /tmp/mycoformatics.yml

# export a list of packages and versions in this environment
ENV PKG_EXP=climush-packages_docker$VERSION_NUM.txt
RUN conda list --export > $PKG_EXP

# We don't want the container to run as root, so make a user account
# (the -m option makes the user's home directory as well).
ARG USER
RUN useradd -m $USER

# add the home folder to PATH:
ENV PATH=/home/$USER:$PATH

# before adding activation of conda environment to .bashrc, copy over the version containing bash aliases
COPY --chown=$USER:users .bashrc /home/$USER/.bashrc

# copy conda env activation to container's bash script so it runs upon entry
RUN echo "source activate mycoformatics" >> /home/$USER/.bashrc
ENV PATH $CONDA_DIR/envs/mycoformatics/bin:$PATH

# mount volumes for pipeline-specific directories, copy content over
WORKDIR /home/$USER
#VOLUME ["/home/$USER/config", "/home/$USER/pipeline"]
#COPY --chown=$USER:users config /home/$USER/config
#COPY --chown=$USER:users pipeline /home/$USER/pipeline
#COPY --chown=$USER:users $PKG_EXP /home/$USER/config/
COPY --chown=$USER:users test.txt /home/$USER/testdir/

# Define an environment var to indicate script is running in a container
ENV CONTAINERIZED=''
USER $USER

# set bash as entrypoint, want the .bashrc file to be first thing to run
ENTRYPOINT ["/bin/bash", "-l", "-c"]
