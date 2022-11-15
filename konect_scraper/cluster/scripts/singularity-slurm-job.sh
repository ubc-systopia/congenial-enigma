#!/bin/bash

#SBATCH --output=%x-%j.out
#SBATCH --nodes=1-1

module load singularity/3.7

# set the variables
if [[ $# -eq 6 ]]; then
    CFG_FILE=$1
    IMAGE=$2
    REPO_HOME=$3
    DATA_DIR=$4
    MODE=$5
    LOCAL_CFG=$6
    ARRAY_OFFSET=0
else
    CFG_FILE=$1
    ARRAY_OFFSET=$2
    IMAGE=$3
    REPO_HOME=$4
    DATA_DIR=$5
    MODE=$6
    LOCAL_CFG=$7
fi

if [[ -v SLURM_ARRAY_TASK_ID ]]; then
    CONFIG_ID=$((ARRAY_OFFSET + SLURM_ARRAY_TASK_ID))
    echo "config id from slurm array id and offset: ${CONFIG_ID} = ${ARRAY_OFFSET} + ${SLURM_ARRAY_TASK_ID}"
else
    CONFIG_ID=$ARRAY_OFFSET
    echo "config id from supplied offset: ${CONFIG_ID} = ${ARRAY_OFFSET}"
fi

offset_line_number=$((CONFIG_ID + 2))
line=$(sed "${offset_line_number}q;d" ${LOCAL_CFG})
echo ${line}
SCRIPTS_DIR=/congenial-enigma/konect_scraper/cluster/scripts/

# echo ${CFG_FILE}
# echo ${IMAGE}
# echo ${CONFIG_ID}
# echo ${REPO_HOME}
# echo ${DATA_DIR}
# echo "${SCRIPTS_DIR}singularity-exec-${MODE}.sh"

singularity exec --bind ${DATA_DIR}:/data,${REPO_HOME}:/congenial-enigma \
    ${IMAGE} \
    ${SCRIPTS_DIR}singularity-exec.sh ${CFG_FILE} ${CONFIG_ID} ${MODE}
