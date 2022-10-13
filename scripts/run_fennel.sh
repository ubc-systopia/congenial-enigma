#!/bin/bash 

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

LOG () {
		# add timestamps to log; add the name of script to log
		# extract the name of the script from $0 
		SCRIPTNAME=$(basename $0)

		echo -e "${GREEN}$(date +"%Y-%m-%d %H:%M:%S") $SCRIPTNAME:${NC} $1"
}

# TODO: add fnl_edges if fnl_edges is implemented
ALL_FENNEL_BIN=("fnl_nodes")

for FENNEL_BIN in "${ALL_FENNEL_BIN[@]}"; do 
		# check if the binary exists 
		if [ -z "$(which $FENNEL_BIN)" ]; then 
				LOG "${RED}ERROR: $FENNEL_BIN not found${NC}"
				exit 1
		fi
done

# 1st arg: path to the input file
# 2nd arg: path to the output file 
# 3rd arg: specifies the number of partitions 
# 4th arg: number of iterations per partition 
# 5th arg: -y to skip the user confirmation

if [ $# -lt 4 ]; then 
		LOG "${RED}ERROR: not enough arguments provided(given $#; expected 4)${NC}"

		echo "Usage: $(basename $0) <input_grpah_file> <output_file> <number_of_partitions> <number_of_iterations_per_partition> [-y]"
		echo ""
		echo -e "\tExample: $(basename $0) ./INPUT_FILE.txt ./FENNEL 10 100"
		echo ""
		echo -e "\tThe example above will run Fennel on the ./INPUT_FILE.txt and produces two files: 1- Last partitioning results in ./fnl_nodes.FENNEL.part.X where X is the number partitions and 2- ./fnl_nodes.FENNEL.part.X.output which includes all fennel's output data and perf information. In this example, setting X is 10 means that the number of partitions starts from 2 and goes to 2^10, i.e. {2, 2^2, 2^3, ..., 2^10}. For each partition, it repeats the experiment 100 times. All repetitions are done by perf tool as 'perf stat -d -d -d --table -r 100'"
		exit 1 
fi

# check if perf is installed
PERF=$(which perf)
if [ -z "$PERF" ]; then
		LOG "perf not found; check if it is installed (run 'which perf')"
		exit 1
fi

# check if perf stat is installed
perf stat -d -d -d --table -r 3 ls > /dev/null 2>&1
if [ $? -ne 0 ]; then
		LOG "perf stat not found; check if it is installed or has the right permissionsrun 'perf stat -d -d -d --table -r 3 ls')"
		exit 1
fi
if [ ! -f $1 ]; then 
		echo "Input file $1 does not exist"
		exit 1 
fi
INPUT_FILE=$1

OUTPUT_FILE=$2

if [ $3 -lt 2 ]; then 
		echo "Number of partitions must be at least 2"
		exit 1 
fi
NUM_PARTITIONS=$3

if [ $4 -lt 1 ]; then echo "Number of iterations per partition must be at least 1"
		exit 1 
fi
NUM_ITERATIONS=$4

# print all input args 
LOG "Input file: $INPUT_FILE"
LOG "Output file: $OUTPUT_FILE"
LOG "Number of partitions: {2^1, ..., 2^$NUM_PARTITIONS}"
LOG "Number of iterations per partition: $NUM_ITERATIONS"
LOG "Fennel binary files: 'fnl_nodes'"

# confirm user wants to continue if -y is not specified
if [ $# -eq 5 ]; then 
		LOG "Auto confirmed(-y option specified)"
else 
		echo -e "${RED}Are you sure you want to continue?${NC}"
		read -p "Enter y to continue: " -n 1 -r
		echo    # (optional) move to a new line
		if [[ ! $REPLY =~ ^[Yy]$ ]]
		then
				exit 1
		fi
fi

# create powers of 2 for NUM_PARTITIONS
for i in $(seq 1 $NUM_PARTITIONS); do 
		num_partitions_array[$i]=$((2**$i))
done

LOG "Number of partitions: ${num_partitions_array[*]}"

# TODO: read the load balance factor from arguments
load_balance_factor=0.03

# iterate over the number of partitions
for FNL_BIN in "${ALL_FENNEL_BIN[@]}"; do 
		LOG "Running $FNL_BIN"

		count=1
		for num_partitions in "${num_partitions_array[@]}"
		do 
				# ./output.fnl_edges.part.X
				OUTPUT_NAME=$OUTPUT_FILE.$FNL_BIN.part.$num_partitions.output
				# log the command to run
				LOG "Running: perf stat -d -d -d --table -r $NUM_ITERATIONS $FNL_BIN $INPUT_FILE $num_partitions $load_balance_factor &> $OUTPUT_NAME"
				perf stat -d -d -d --table -r $NUM_ITERATIONS $FNL_BIN $INPUT_FILE $num_partitions $load_balance_factor &> $OUTPUT_NAME

		# check if fennel ran successfully
		if [ $? -ne 0 ]; then 
				LOG "Fennel failed: error logged to ./error.log"
				# log the file name in an error file 
				echo $OUTPUT_NAME >> ./error.log
		else 
				LOG "Fennel succeeded"
		fi
		# done
		count=$((count+1))
done
done

LOG "If there is any error, it is logged to ./error.log"
LOG "DONE"
