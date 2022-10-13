#!/bin/bash

help() {
		echo "Usage: $(basename $0) <metis log file path>"
		echo "Reads a metis log output and extracts data out of its lines. The output is a csv file at <metis log file path>.time.csv"
		echo "\tExample: $(basename $0) /home/username/metis.log"
}

# first arg is the name of the file 
input_log_file=$1

# check if the file exists
if [ ! -f $input_log_file ]; then
    echo "File $input_log_file does not exist"
    exit 1
fi

output_csv_file=${input_log_file}.time.csv
# check if the file exists; remove it if it does and ask for confirmation
if [ -f $output_csv_file ]; then
    echo "File $output_csv_file already exists (will be overwritten)"
	echo "the first 10 lines of the file are:"
	head -n 10 $output_csv_file
	echo
    read -p "Do you want to overwrite it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
	rm $output_csv_file
    else
	exit 1
    fi
fi

echo "Creating file $output_csv_file"
grep "Running: perf stat" $input_log_file | cut -d ' ' -f 1,2,3,5,6,18,19 |  sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | tr -s ':' ' ' | tr -s '[ \t]+' ',' > $output_csv_file

echo "adding the header"
sed -i '1s/^/date,h,m,s,script_name,perf_cmd,stat_cmd,dataset_name,num_parts\n/' $output_csv_file 

echo "Done: results in $output_csv_file"
