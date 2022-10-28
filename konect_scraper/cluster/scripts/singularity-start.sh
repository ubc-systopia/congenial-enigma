#!/bin/bash

LOCAL_CONFIG=$1
MOUNTED_CONFIG=$2
LOGDIR=$3
IMAGE=$4
REPO_HOME=$5
DATA_DIR=$6
MODE=$7
# get the number of configuration lines
NUM_LINES=$(wc -l ${LOCAL_CONFIG} | cut -f 1 -d ' ')
NUM_CONFIGS=$((NUM_LINES - 2)) # include csv header
CFG_ID_END=$((NUM_CONFIGS - 1))

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

OUT_FILE=${LOGDIR}/%A-%a.out
JOB_FILE=${SCRIPT_DIR}/singularity-slurm-job.sh
echo "log directory: " ${LOGDIR}
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
  # todo depending on MODE, pass additional arguments to sbatch to ask for appropriate --cpus-per-task, --mem
  echo "sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${MOUNTED_CONFIG} ${ARRAY_START} ${IMAGE} ${REPO_HOME} ${DATA_DIR} ${MODE}"
  sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${MOUNTED_CONFIG} ${ARRAY_START} ${IMAGE} ${REPO_HOME} ${DATA_DIR} ${MODE}

  # calculate the new array start
  ARRAY_START=$((ARRAY_START + MAX_ARRAY_JOBS))
  # calculate the remaining configurations
  NUM_CONFIGS=$((NUM_CONFIGS - MAX_ARRAY_JOBS))

  # wait a little bit, 2 seconds
  sleep 2
done

ARRAY_END=$((NUM_CONFIGS - 1))

# enqueue the remainder
echo "sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${MOUNTED_CONFIG} ${ARRAY_START} ${IMAGE} ${REPO_HOME} ${DATA_DIR} ${MODE}"
sbatch --array=0-${ARRAY_END}  -o ${OUT_FILE}  ${JOB_FILE} ${MOUNTED_CONFIG} ${ARRAY_START} ${IMAGE} ${REPO_HOME} ${DATA_DIR} ${MODE}