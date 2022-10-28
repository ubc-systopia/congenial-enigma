#!/bin/bash

#SBATCH --time=00-00:05:00     	# DD-HH:MM:SS
#SBATCH --job-name=congenial-engima-download
#SBATCH --output=%x-%j.out
module load singularity/3.7

# set the variables
if [[ $# -eq 4 ]]; then
    CFG_FILE=$1
    IMAGE=$2
    REPO_HOME=$3
    DATA_DIR=$4
    ARRAY_OFFSET=0
else
    CFG_FILE=$1
    ARRAY_OFFSET=$2
    IMAGE=$3
    REPO_HOME=$4
    DATA_DIR=$5
fi

if [[ -v SLURM_ARRAY_TASK_ID ]]; then
    CONFIG_ID=$((ARRAY_OFFSET + SLURM_ARRAY_TASK_ID))
    echo "config id from slurm array id and offset: ${CONFIG_ID} = ${ARRAY_OFFSET} + ${SLURM_ARRAY_TASK_ID}"
else
    CONFIG_ID=$ARRAY_OFFSET
    echo "config id from supplied offset: ${CONFIG_ID} = ${ARRAY_OFFSET}"
fi

echo ${CFG_FILE}
echo ${IMAGE}
echo ${CONFIG_ID}
echo ${REPO_HOME}
echo ${DATA_DIR}

singularity exec --bind ${DATA_DIR},${REPO_HOME}:/congenial-enigma ${IMAGE} /congenial-enigma/konect_scraper/cluster/scripts/download/singularity-exec.sh ${CFG_FILE} ${CONFIG_ID}
#SBATCH --constraint=broadwell 	# Request Broadwell processor
#SBATCH --mem=1G       	# Memory proportional to GPUs: 32000 Cedar, 47000 Béluga, 64000 Graham.
#SBATCH --time=0-01:00     	# DD-HH:MM:SS
#SBATCH --nodes=1-1
