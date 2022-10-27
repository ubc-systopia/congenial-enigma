#!/bin/bash

CONFIG_FILE=$1
LOGDIR=$2
IMAGE=$3
# get the number of configuration lines
NUM_LINES=$(wc -l ${CONFIG_FILE} | cut -f 1 -d ' ')
NUM_CONFIGS=$((NUM_LINES - 2)) # include csv header
CFG_ID_END=$((NUM_CONFIGS - 1))

OUT_FILE=${LOGDIR}/download-%A-%a.out
JOB_FILE=$(pwd)/singularity-start.sh
echo "slurm output file: ${OUT_FILE}"
echo "slurm job file: ${JOB_FILE}"

echo "---------------------------------------"
echo "enqueue ${NUM_CONFIGS} jobs..."

MAX_ARRAY_JOBS=$(scontrol show config | grep MaxArraySize | cut -d= -f2)
ARRAY_START=0

ARRAY_END=$((MAX_ARRAY_JOBS - 1))
while [[ $NUM_CONFIGS -gt $MAX_ARRAY_JOBS ]]; do
  # calculate the end of the array

  # enqueue the batch for the current values
  echo "sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${CONFIG_FILE} ${ARRAY_START} ${IMAGE}"
  sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${CONFIG_FILE} ${ARRAY_START} ${IMAGE}

  # calculate the new array start
  ARRAY_START=$((ARRAY_START + MAX_ARRAY_JOBS))
  # calculate the remaining configurations
  NUM_CONFIGS=$((NUM_CONFIGS - MAX_ARRAY_JOBS))

  # wait a little bit, 2 seconds
  sleep 2
done

ARRAY_END=$((NUM_CONFIGS - 1))

# enqueue the remainder
echo "sbatch --array=0-${ARRAY_END} -o ${OUT_FILE}  ${JOB_FILE} ${CONFIG_FILE} ${ARRAY_START} ${IMAGE}"
sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${CONFIG_FILE} ${ARRAY_START} ${IMAGE}