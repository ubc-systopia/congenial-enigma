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

		BASENAME=$(basename $0)
		# date time format + message 
		echo -e "${COLOR}$(date +"%Y-%m-%d %H:%M:%S") $BASENAME: $TEXT ${NC}"
}

process_file () {
		file=$(cat $1 | tr -s '[ \t]*' ' ')
		file_name=$(basename $1)
		echo $file_name
		# if file is empty, skip it
		if [ -z "$file" ]; then
				LOG "File $1 is empty, skipping" $YELLOW
				return
		fi
		# exctract vol or cut from file name; it is the the second element from the last part of the file name
		vol_or_cut=$(awk -F. '{print $(NF-1)}' <<< $file_name)
		# echo $vol_or_cut

		names=( $(echo $file | grep -Po "Name: \K[^, ]*") )
		# echo ${names[@]} "length:" ${#names[@]}
		# echo "-----------------"
		num_vertices=( $(echo $file | grep -Po "#Vertices: \K[^, .]*") )
		num_edges=( $(echo $file | grep -Po "#Edges: \K[^, .]*") )

		num_parts=( $(echo $file | grep -Po "#Parts: \K[^, .]*") )

		# communication volume and edge-cut for each partition 
		comm_vols=( $(echo $file | grep -Po "communication volume: \K[^, .]*") )
		edge_cuts=( $(echo $file | grep -Po "Edgecut: \K[^, .]*") )

		# timing information 
		part_times=( $(echo $file | grep -Po "Partitioning: \K[^, ]*") )
		io_times=( $(echo $file | grep -Po "I/O: \K[^, ]*") )
		# TODO: add total time
		reporting_times=( $(echo $file | grep -Po "Reporting: \K[^, ]*") )

		# memory information 
		max_mem_used_list=( $(echo $file | grep -Po "Max memory used: \K[^, ]*") )
		rusage_ru_maxrss_list=( $(echo $file | grep -Po "rusage.ru_maxrss: \K[^, ]*") )
		

		# if any of these are empty list then something went wrong 
		if [ -z "$names" ] || [ -z "$num_vertices" ] || [ -z "$num_edges" ] || [ -z "$num_parts" ] || [ -z "$comm_vols" ] || [ -z "$edge_cuts" ] || [ -z "$part_times" ] || [ -z "$io_times" ] || [ -z "$max_mem_used_list" ] || [ -z "$rusage_ru_maxrss_list" ]; then
				LOG "File $1 is missing some information, exiting" $RED
				exit 1
		fi

		expr_length=${#names[@]}
		expr_length=$((expr_length-1))

		# get the perf output from the last 10 lines 
		perf_output=$(tail -n 7 $file_name)
		all_perf_times=( $(echo $perf_output | tr -s '[ \t]*' ' ' | grep -Po "\K([+-]?\d+.\d+) ") )

		# last index of all_perf_times is the total time
		echo "Total time: ${all_perf_times[@]}"
		perf_total_time=${all_perf_times[-2]}
		perf_std=${all_perf_times[-1]}
		
		echo "Total time: $perf_total_time" "std: $perf_std"
		# iterate from 0 to the lengthg of names 
		for i in $(seq 0 ${expr_length}); do 
				echo ${i}
				# if the file is valid, then we can print the results
				name=$(basename ${names[$i]})
				write_to_file "${name}" "${num_parts[$i]}" "${num_vertices[$i]}" "${num_edges[$i]}" "${part_times[$i]}" "${io_times[$i]}" "${comm_vols[$i]}" "${edge_cuts[$i]}" "${max_mem_used_list[$i]}" "${rusage_ru_maxrss_list[$i]}" $output_file "${vol_or_cut}" "${reporting_times[$i]}" "${all_perf_times[$i]}" "${perf_total_time}" "${perf_std}"
		done
		# exit 1 # TDOO: remove this
}

write_to_file () {
		# write to a csv file
		# 11: is the file output
		# 12: "obj_type": is the vol or cut
		# 13.. : "report_time" "perf_time" "perf_total_time" "perf_std"
		echo "$1, ${12}, $2, $3, $4, $5, ${13}, ${14}, ${15}, ${16}, $6, $7, $8, $9, ${10}" >> ${11}
}

# 1st arg is the output file name 
# last arg is -y to overwrite
if [ $# -lt 2 ]; then
		LOG "Given args: $@"
		BASENAME=$(basename $0)
		LOG "Usage: $BASENAME <search directory for metis output> <output file name> [-y]"
		LOG "Example: $BASENAME /home/user/metis_output/ /home/user/metis_output/metis_output.csv: searches for all .output files, which are Metis output data, and writes them to the output file"
		exit 1
fi

TARGET_DIR=$1
if [ ! -d $TARGET_DIR ]; then
		LOG "Target directory $TARGET_DIR does not exist" ${RED}
		exit 1
fi

# if there are no files, then exit
if [ -z "$(find $TARGET_DIR -name "*.output")" ]; then
		LOG "No .output files found in $TARGET_DIR" ${RED}
		exit 1
fi

LOG "Total number of files to process: $(find $TARGET_DIR -name "*.output" | wc -l)"
LOG "Processing files in $TARGET_DIR directory"

# if -y is given, then overwrite the output file 
output_file=$2
if [ $# -eq 3 ] && [ $3 == "-y" ]; then 
		LOG "Overwriting output file $output_file (-y given)"
		true > $output_file
else
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
fi

LOG "Starting processing"
write_to_file "name" "num_parts" "num_vertices" "num_edges" "part_time" "io_time"  "comm_vol" "edge_cut" "max_mem_used" "rusage_ru_maxrss" $output_file "obj_type" "report_time" "perf_time" "perf_total_time" "perf_std"
# find $TARGET_DIR -name "*.output"  | while read file; do process_file "$file"; done
for file in $(find $TARGET_DIR -name "*.output"); do
		echo "Processing file $file"
		process_file "$file" 
done


LOG "Successfully processed files." ${GREEN}
LOG "Output file: $output_file" ${GREEN}
