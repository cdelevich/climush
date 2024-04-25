#!/usr/bin/env bash

# get the current working directory
echo "PWD=$(pwd)" >> .env

# pull the compose.yaml file from GitHub to set up Docker image
if which wget >/dev/null ; then
  echo "Downloading via wget."
  wget -O compose.yaml https://github.com/cdelevich/climush/docker-containers/bioinformatics/deploy/compose.yaml
elif which curl >/dev/null ; then
  echo "Downloading via curl."
  curl -O https://github.com/cdelevich/climush/docker-containers/bioinformatics/deploy/compose.yaml
else
  echo "Cannot download, neither wget nor curl is available."
fi

# create an empty sequence directory, with full permissions
mkdir -m 777 sequences

# create a config directory with full permissions
mkdir -m 777 config

# download the configuration file to this directory
cd config
if which wget >/dev/null ; then
  echo "Downloading via wget."
  wget -O https://github.com/cdelevich/climush/bioinformatics-pipeline/config/climush-bioinfo_user-settings.toml
elif which curl >/dev/null ; then
  echo "Downloading via curl."
  curl -O https://github.com/cdelevich/climush/bioinformatics-pipeline/config/climush-bioinfo_user-settings.toml
else
  echo "Cannot download, neither wget nor curl is available."
fi
# had to do this annoying way where I change directory, would not work otherwise
cd ../

# run Docker container from Docker image via compose.yaml
docker compose run climush-bioinformatics

# if on HPC cluster
#module load apptainer
#apptainer run docker://cdelevich/climush-bioinfo:latest