#!/usr/bin/env bash

# get the current working directory
echo "PWD=$(pwd)" >> .env

# pull the compose.yaml file from GitHub to set up Docker image
#wget -O compose.yaml https://github.com/cdelevich/climush.git#main:/docker-containers/bioinformatics/deploy/compose.yaml
curl -o . https://github.com/cdelevich/climush.git#main:/docker-containers/bioinformatics/deploy/compose.yaml

# create an empty sequence directory
mkdir sequences

# create a config directory and add configuration file into it
mkdir config
#wget -O /config/climush-bioinfo_pipeline-settings.toml https://github.com/cdelevich/climush.git#main:/bioinformatics-pipeline/config/climush-bioinfo_pipeline-settings.toml
curl -o /config https://github.com/cdelevich/climush.git#main:/bioinformatics-pipeline/config/climush-bioinfo_pipeline-settings.toml

# run Docker container from Docker image via compose.yaml
docker compose run climush-bioinformatics