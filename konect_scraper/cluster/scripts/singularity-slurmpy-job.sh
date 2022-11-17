#!/bin/bash

module load singularity/3.7

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

echo "singularity exec --bind ${DATA_DIR}:/data,${REPO_HOME}:/congenial-enigma \
        ${IMAGE} \
        ${SCRIPTS_DIR}singularity-exec.sh ${CFG_FILE} ${CONFIG_ID} ${MODE}"

# singularity exec --bind ${DATA_DIR}:/data,${REPO_HOME}:/congenial-enigma \
#     ${IMAGE} \
#     ${SCRIPTS_DIR}singularity-exec.sh ${CFG_FILE} ${CONFIG_ID} ${MODE}
