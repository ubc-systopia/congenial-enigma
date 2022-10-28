#!/bin/bash
# set the variables
CFG_FILE=$1     # the configuration csv file
CFG_INDEX=$2    # the index into the csv file
MODE=${3}       # execution mode -- one of {download, preprocess, reorder, plot, pr_expt}

cd /congenial-enigma

echo "python3.10 -m konect_scraper.cluster.slurm.${MODE} ${CFG_FILE} ${CFG_INDEX} /data"
python3.10 -m \
    konect_scraper.cluster.slurm.${MODE} \
        ${CFG_FILE} \
        ${CFG_INDEX} \
        /data
