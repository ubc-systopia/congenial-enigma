#!/bin/bash
#!/bin/bash

#################################################################################################
#
# Params:
#  -r <REQSPATH> absolute path to local requirements.txt 
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


# while getopts ":r:d:" opt; do
#   case ${opt} in
#     r )
#     #   REQS_PATH=$(readlink -f $OPTARG)
#       REQS_PATH=${OPTARG}
#       ;;
#     \? )
#       echo "Invalid option: $OPTARG" 1>&2
#       ;;
#     : )
#       echo "Invalid option: $OPTARG requires an argument" 1>&2
#       ;;
#   esac
# done
# echo ${REQS_PATH} 

docker build -t install-deps .
    # --build-arg reqs_path=${REQS_PATH} \
    # .