#!/bin/bash
# set the variables
CFG_FILE=$1     # the configuration csv file
CFG_INDEX=$2    # the index into the csv file

cd /congenial-enigma

echo "python3.10 -m konect_scraper.cluster.slurm.download ${CFG_FILE} ${CFG_INDEX}"
python3.10 -m konect_scraper.cluster.slurm.download ${CFG_FILE} ${CFG_INDEX}
