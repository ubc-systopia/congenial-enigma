#!/bin/bash 

NC='\033[0m'
BLUE='\033[1;34m'

LOG () {
		# date, name of the script, and the message
		SHORTNAME=$(basename $0)
		echo -e "${BLUE}$(date +"%Y-%m-%d %H:%M:%S") $SHORTNAME:${NC} ${1}"
}

RUN_FENNEL_SH=$(which run_fennel.sh)
if [ ! -f $RUN_FENNEL_SH ]; then
		LOG "run_fennel.sh not found; place it in your PATH"
		exit 1
fi

# 1st arg is the number of partitions
# 2nd arg is the number of iterations per partition (default is 1)

if [ $# -eq 0 ]; then
		echo "Usage: $(basename $0) <num. partitions: length of powers of 2 list> [<number of iterations per partition>: default is 1]"
		echo "Example: $(basename $0) 4 3 : will partition all .metis graphs in the current directory in {2^1, 2^2, 2^3, 2^4} partitions, repeating 3 times per number of partition"
		exit 1
fi

NUM_PARTITIONS=$1
# make sure it is a number 
if ! [[ $NUM_PARTITIONS =~ ^[0-9]+$ ]]; then
		echo "Error: $NUM_PARTITIONS is not a number (given $NUM_PARTITIONS)"
		exit 1
fi

if [ $# -eq 2 ]; then
		NUM_ITERATIONS=$2
else
		NUM_ITERATIONS=1
fi
# make sure it is a number
if ! [[ $NUM_ITERATIONS =~ ^[0-9]+$ ]]; then
		echo "Error: $NUM_ITERATIONS is not a number (given $NUM_ITERATIONS)"
		exit 1
fi

# confirm the input arg by the user 
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
to_be_deleted_extensions=(*fennel*part*)

# show the first 10 files to be deleted if ther is any 
for ext in "${to_be_deleted_extensions[@]}"
do
		# find all files with the extension
		files_to_be_deleted=$(find $PWD -name "$ext" -print )
		if [ -n "$files_to_be_deleted" ]; then
				LOG "The following files will be deleted (first 10):"
				echo "$files_to_be_deleted" | head -n 10
				break
		fi
done

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
		find $PWD -name $extension -delete
done

LOG "Deleted all *.output and .part.* files"
LOG "Starting script"

# find all .metis files in the current directory
METIS_FILES=$(find . -name "*.metis")
# sort the files by its size (in bytes)
METIS_FILES=$(find . -name "*.metis" -type f -exec ls -lhSr {} + | tr -s '[ \t]' ' ' | cut -d ' ' -f 9)

# print the files;
LOG "Found the following metis files (sorted by their size asc):"
echo "${METIS_FILES}"

# run run_fennel.sh for each file in the list;
total_metis_files=$(echo "$METIS_FILES" | wc -l)
count=1
for file in $METIS_FILES; do
		LOG "Running Fennel on ${file} ($count/$total_metis_files metis files)"

		$RUN_FENNEL_SH $file $NUM_PARTITIONS $NUM_ITERATIONS -y

		count=$((count+1))
done
