#!/bin/bash 

NC='\033[0m'
BLUE='\033[1;34m'

LOG () {
		# date, name of the script, and the message
		echo -e "${BLUE}$(date +"%Y-%m-%d %H:%M:%S") ${0}:${NC} ${1}"
}
# 1st arg is the number of partitions
# 2nd arg is the number of iterations per partition (default is 1)

if [ $# -eq 0 ]; then
		echo "Usage: $0 <num. partitions: length of powers of 2 list> [<number of iterations per partition>: default is 1]"
		echo "Example: $0 4 3 : will partition all .metis graphs in the current directory in {2^1, 2^2, 2^3, 2^4} partitions, repeating 3 times per number of partition"
		exit 1
fi

NUM_PARTITIONS=$1
if [ $# -eq 2 ]; then
		NUM_ITERATIONS=$2
else
		NUM_ITERATIONS=1
fi

# confirm the input arg by the user 
LOG "Input arguments: $NUM_PARTITIONS $NUM_ITERATIONS"
echo "Number of partitions: $NUM_PARTITIONS"
echo "Number of iterations per partition: $NUM_ITERATIONS"
read -p "Is this correct? (y/n) " -n 1 -r
echo    # (optional) move to a new line

if [[ $REPLY =~ ^[Yy]$ ]]
then
		LOG "confirmed"
else
		LOG "Exiting..."
		exit 1
fi

# confirm that all .output and .part.* files will be deleted recursivelly in this directory
to_be_deleted_extensions=("*.output" "*.part.*")

read -p "All $to_be_deleted_extensions files will be removed. Are you sure you want to delete all $to_be_deleted_extensions files in this directory? (y/n) " -n 1 -r
echo    # (optional) move to a new line

if [[ $REPLY =~ ^[Yy]$ ]]
then
		LOG "Deleting..."
else
		LOG "Exiting..."
		exit 1
fi

for extension in "${to_be_deleted_extensions[@]}"
do
		LOG "Deleting files with extension $extension"
		find . -name $extension -delete
done

LOG "Deleted all .output and .part.* files"
LOG "Starting script"

# find all .metis files in the current directory
METIS_FILES=$(find . -name "*.metis")
# sort the files by its size (in bytes)
METIS_FILES=$(find . -name "*.metis" -type f -exec ls -lhSr {} + | tr -s '[ \t]' ' ' | cut -d ' ' -f 9)

# print the files;
LOG "Found the following metis files (sorted by their size asc):"
echo "${METIS_FILES}"

# run run_metis.sh for each file in the list;
total_metis_files=$(echo "$METIS_FILES" | wc -l)
count=1
for file in $METIS_FILES; do
		LOG "Running metis on ${file} ($count/$total_metis_files)"
		OUTPUT_FILE="${file}"
		./scripts/run_metis.sh $file $OUTPUT_FILE $NUM_PARTITIONS $NUM_ITERATIONS -y
		count=$((count+1))
done
