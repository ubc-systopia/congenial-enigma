#!/bin/bash 

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

LOG () {
		# add timestamps to log; add the name of script to log
		echo -e "${GREEN}$(date +"%Y-%m-%d %H:%M:%S") $0:${NC} $1"
}

# metis binary file 
METIS_PARTITIONER=$(which gpmetis)

if [ -z "$METIS_PARTITIONER" ]; then
		LOG "gpmetis not found; check if it is installed (run 'which gpmetis')"
		exit 1
fi

# 1st arg: path to the input file
# 2nd arg: path to the output file 
# 3rd arg: specifies the number of partitions 
# 4th arg: number of iterations per partition 
# 5th arg: -y to skip the user confirmation

if [ $# -lt 4 ]; then 
		echo "Usage: $0 <input_grpah_file> <output_file> <number_of_partitions> <number_of_iterations_per_partition> [-y]"
		echo ""
		echo -e "\tExample: $0 ./input.txt ./output.txt 10 100"
		echo ""
		echo -e "\tThe example above will run the algorithm on the ./input.txt for Y times and produces two files: 1- last partitioning results in ./output.txt.part.X where X is the number partitions and 2- ./output.txt.part.X.output which includes all metis output data and perf information. In this example, setting X to 10 means that the number of partitions starts from 2 and goes to 2^10, i.e. {2, 2^2, 2^3, ..., 2^10}, and Y starts from 1 and goes to 100. All repetitions are done by perf tool as 'perf stat -d -d -d --table -r 3'"
		exit 1 
fi

if [ ! -f $1 ]; then 
		echo "Input file $1 does not exist"
		exit 1 
fi
INPUT_FILE=$1

# if [ -f $2 ]; then 
# 		echo "Output file $2 already exists"
# 		exit 1 
# fi
OUTPUT_FILE=$2

if [ $3 -lt 2 ]; then 
		echo "Number of partitions must be at least 2"
		exit 1 
fi
NUM_PARTITIONS=$3

if [ $4 -lt 1 ]; then 
		echo "Number of iterations per partition must be at least 1"
		exit 1 
fi
NUM_ITERATIONS=$4


# print all input args 
LOG "Input file: $INPUT_FILE"
LOG "Output file: $OUTPUT_FILE"
LOG "Number of partitions: {2^1, ..., 2^$NUM_PARTITIONS}"
LOG "Number of iterations per partition: $NUM_ITERATIONS"
LOG "Metis binary file: $METIS_PARTITIONER" 

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

# iterate over the number of partitions
count=1
for num_partitions in "${num_partitions_array[@]}"
do 
		LOG "Number of partitions: $num_partitions ($count/$NUM_PARTITIONS)"
		# iterate over the number of iterations per partition
		for ((i=1; i<$NUM_ITERATIONS+1; i++))
		do
				LOG "Iteration: $i"
				OUTPUT_NAME=$OUTPUT_FILE.part.$num_partitions.iter.$i.output
				LOG "Output file: $OUTPUT_NAME"

		# run metis 
		LOG "Running metis"
		$METIS_PARTITIONER -objtype vol -seed 0 $INPUT_FILE $num_partitions > $OUTPUT_NAME 

		# check if metis ran successfully
		if [ $? -ne 0 ]; then 
				LOG "Metis failed: error logged to ./error.log"
				# log the file name in an error file 
				echo $OUTPUT_NAME >> ./error.log
		else 
				LOG "Metis succeeded"
		fi
done
count=$((count+1))
done

LOG "If there is any error, it is logged to ./error.log"
LOG "DONE"
