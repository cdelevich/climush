#!/bin/bash

# convert input value for the log frequency from minutes to seconds
#LOGFREQ=`expr $2 \* 60`
LOGFREQ=$2
DEFAULT_LOGFREQ=`expr 30 \* 60`  # couldn't combine minutes with 'm' string, so converting to seconds

# set a maximum number of times that this script will log progress (20 = 10 hrs)
MAXITER=20

# get the absolute path from the input path
LOGDIR=$(realpath $1)  # this works, confirmed with echo

# create a variable for the name of the log file
FOUT="${LOGDIR}/job_progress.log"

# change directory to the input path, positional variable 1
# if this failed, then print error message and exit
#cd ${LOGDIR} || printf "An input path is required for this function to work." & exit
cd ${LOGDIR} || exit

# print directory; will also clear any previous log file versions by writing (>) and not appending (>>)
echo "Logging file size and number of files every ${LOGFREQ:=$DEFAULT_LOGFREQ} seconds in:
   ${LOGDIR}
" > $FOUT

# define a function that will log the number of files and total directory size
function log_progress () {
  # get the log time, number of files, and the size of the main directory
  local TIME=$(date +"%Y-%m-%d %T")
  local NFILE=$(find . -maxdepth 1 -type f | wc -l | tr -d ' ')
  local FSIZE=$(du -sh | cut -f1)

  # format these values and write out to the log file
  echo "${TIME}
  number of files = ${NFILE}
  directory size =  ${FSIZE}
  " >> $FOUT
}

# run the log function every x minutes; defaults to every 30 min if no time provided from command line
#while true; do log_progress & sleep ${LOGFREQ:=30m}; done
#while true; do log_progress & sleep 30m; done
for i in `seq 1 $MAXITER`; do log_progress & sleep ${LOGFREQ:=$DEFAULT_LOGFREQ}; done