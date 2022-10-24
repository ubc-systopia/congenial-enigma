#!/bin/bash

#################################################################################################
#
# Params:
#  -r <REPODIR> mount the congenial repo directory in docker image
#  -d <DATADIR> mount the data directory to store graphs, results
#
# Example usage:
#  to run using a pre-generated configuration file
#
#  $ bash build_and_run_docker.sh -r results-v1 -d 
#
#  to generate the configuration file for running it:
#
# 
################################################################################################


while getopts ":r:d:" opt; do
  case ${opt} in
    r )
      REPO_DIR=$(readlink -f $OPTARG)
      ;;
    d )
      DATA_DIR=$(readlink -f $OPTARG)
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
echo ${REPO_DIR} ${DATA_DIR}