#!/bin/bash

module load singularity

singularity build climush-test.sif docker-archive://climush-test.tar 

singularity instance start climush_test.sif climushtest 

singularity shell instance://climushtest
