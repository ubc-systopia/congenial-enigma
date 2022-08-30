#!/bin/bash 

RED='\033[0;31m'
GREEN='\033[0;32m'
GREY='\033[0;37m'
YELLOW='\033[1;33m'
NC='\033[0m'

LOG () {
		COLOR=$2
		TEXT=$1

		if [ -z $COLOR ]; then
			COLOR=$GREY
		fi

		# date time format + message 
		echo -e "${COLOR}$(date +"%Y-%m-%d %H:%M:%S") $0: $TEXT ${NC}"
}

process_file () {
		file=$(cat $1 | tr -s '[ \t]*' ' ')

		name=$(echo $file | grep -Po "Name: \K[^,]*")

		num_vertices=$(echo $file | grep -Po "#Vertices: \K[^,]*")
		num_edges=$(echo $file | grep -Po "#Edges: \K[^,]*")
		num_parts=$(echo $file | grep -Po "#Parts: \K[^,\n ]*")

		# timing information 
		part_time=$(echo $file | grep -Po "Partitioning:[ ]*\K[0-9]*.[0-9]*")
		io_time=$(echo $file | grep -Po "I/O:[ ]*\K[0-9]*.[0-9]*")

		# memory information 
		max_mem_used=$(echo $file | grep -Po "Max memory used:[ ]*\K[0-9]*.[0-9]*")
		rusage_ru_maxrss=$(echo $file | grep -Po "rusage\.ru_maxrss:[ ]*\K[0-9]*.[0-9]*")

		# if any of these are empty, then the file is not valid
		if [ -z "$name" ] || [ -z "$num_vertices" ] || [ -z "$num_edges" ] || [ -z "$num_parts" ] || [ -z "$part_time" ] || [ -z "$io_time" ] || [ -z "$max_mem_used" ] || [ -z "$rusage_ru_maxrss" ]; then
				LOG "File $1 is not valid (missing fields, empty files, truncated, ... ) - Skipped" ${YELLOW}
				return 1
		fi


		# if the file is valid, then we can print the results
		write_to_file $name $num_parts $num_vertices $num_edges $part_time $io_time $max_mem_used $rusage_ru_maxrss $output_file 
}

write_to_file () {
		# write to a csv file 
		echo "$1,$2,$3,$4,$5,$6,$7,$8" >> $9
		# if there is an error, print it
		if [ $? -ne 0 ]; then
				echo "Error writing to file $9"
				echo "$1,$2,$3,$4,$5,$6,$7,$8" 
		fi
}

# 1st arg is the output file name 
# last arg is -y to overwrite
if [ $# -lt 1 ]; then
		LOG "Given args: $@"
		LOG "Usage: $0 <output file name> [-y]"
		exit 1
fi

output_file=$1

if [ -f $output_file ]; then
		LOG "Output file $output_file exists. Overwrite? [y/n]"
		read -n 1 -s answer
		
		if [ $answer != "y" ]; then
				LOG "Exiting"
				exit 1
		fi

		true > $output_file
else
		LOG "Creating output file $output_file"
		touch $output_file
fi

write_to_file "name" "num_parts" "num_vertices" "num_edges" "part_time" "io_time" "max_mem_used" "rusage_ru_maxrss" $output_file

LOG "Total number of files to process: $(find . -name "*.output" | wc -l)"
LOG "Processing files in current directory"

find . -name "*.output" | while read file; do process_file "$file"; done
LOG "Successfully processed files" ${GREEN}
