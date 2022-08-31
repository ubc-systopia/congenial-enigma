#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

LOG () {
		# the second arg is the color of the message if it is present
		if [ -z "$2" ]; then
			echo -e "$(date +"%Y-%m-%d %H:%M:%S") $1"
		else
			echo -e "$2$(date +"%Y-%m-%d %H:%M:%S") $1$NC"
		fi
}

# 1st arg is the name of the input file
if [ $# -eq 0 ]; then
		LOG "No input file specified"
		LOG "USAGE: $0 <input file: edge list file> [delete previous outputs: -y | -n(default)]"
		LOG "EXAMPLE: $0 input.txt -y"
		exit 1
fi
INPUT_FILE=$1
REMOVE_FILES=false 
if [ $# -eq 2 ]; then
		if [ $2 == "-y" ]; then
				REMOVE_FILES=true
		fi
fi
LOG "Input file: $INPUT_FILE"
if [ ! -f $INPUT_FILE ]; then
		LOG "Input file does not exist"
		exit 1
fi

# confirm that the files are removed before running the script

if [ $REMOVE_FILES == true ]; then
		LOG "Removing all .sorted and .flattened and .uniq files"
		read -p "Are you sure you want to remove all .sorted and .flattened and .uniq files? (y/n) " -n 1 -r
		echo    # (optional) move to a new line
		# if 2nd arg is y, then remove the files

		if [[ $REPLY =~ ^[Yy]$ ]]; then
				rm -f *.source_sorted *.flattened *.uniq
		fi
fi

# remove all lines that starts with # or % if $INPUT_FILE.tmp doesn't exist
if [ ! -f $INPUT_FILE.tmp ]; then
		LOG "Removing all lines that starts with # or %"
		sed '/^#/d;/^%/d' $INPUT_FILE > $INPUT_FILE.tmp
else 
		LOG "File $INPUT_FILE.tmp already exists" $YELLOW
		# print the size of the file 
		LOG "Size of $INPUT_FILE.tmp: `du -h $INPUT_FILE.tmp | awk '{print $1}'`" $YELLOW
		LOG "Skipping removing all lines that starts with # or % (already done)" $YELLOW
fi

LOG "Cleaned input file: $INPUT_FILE.tmp"
CLEAN_INPUT_FILE=$INPUT_FILE.tmp

# flatten the input file into a single line | sort the lines by spliting the line by space
LOG "Flattening input file"
cat $CLEAN_INPUT_FILE | tr '\n' ' ' | tr -s "[\t ]" '\n'  > $CLEAN_INPUT_FILE.flattened

LOG "Sorting the flattened input file"
sort --parallel=28 -nuo $CLEAN_INPUT_FILE.node_ids.uniq $CLEAN_INPUT_FILE.flattened


LOG "Sorted flattened input file: $CLEAN_INPUT_FILE.node_ids.uniq"

LOG "Checking if $CLEAN_INPUT_FILE.node_ids.uniq is sorted" $YELLOW
sort -nc $CLEAN_INPUT_FILE.node_ids.uniq 
if [ $? -eq 0 ]; then
		LOG "$CLEAN_INPUT_FILE.node_ids.uniq is sorted" $GREEN
else
		LOG "$CLEAN_INPUT_FILE.node_ids.uniq is not sorted, try again" $RED
		exit 1
fi

LOG "Check if the $CLEAN_INPUT_FILE is sorted based on the source id (first column)" $YELLOW
sort -c -n -k1,1 $CLEAN_INPUT_FILE 

if [ $? -eq 0 ]; then
		LOG "The $CLEAN_INPUT_FILE is sorted based on the source id (first column)" $GREEN
		mv $CLEAN_INPUT_FILE $CLEAN_INPUT_FILE.source_sorted
else
		LOG "The $CLEAN_INPUT_FILE is not sorted based on the source id (first column)" $RED
		LOG "Sorting the $CLEAN_INPUT_FILE based on the source id (first column)" $YELLOW
		sort -n -k1,1 -o $CLEAN_INPUT_FILE.source_sorted $CLEAN_INPUT_FILE
fi

SOURCE_SORTED_FILE=$CLEAN_INPUT_FILE.source_sorted

LOG "Sort the $SOURCE_SORTED_FILE based on the target id (second column)" $YELLOW 
sort -n -k2,2 -o $SOURCE_SORTED_FILE.reverse $SOURCE_SORTED_FILE

LOG "Replacing tabs and spaces with a single space"
tr -s "[\t ]" " " < $SOURCE_SORTED_FILE > $SOURCE_SORTED_FILE.cleaned 
rm $SOURCE_SORTED_FILE 
mv $SOURCE_SORTED_FILE.cleaned $SOURCE_SORTED_FILE 

LOG "Removing the flattened input file"
rm -f $CLEAN_INPUT_FILE.flattened

LOG "Removing the tmp input file"
rm -f $CLEAN_INPUT_FILE

LOG "Three files are created: $SOURCE_SORTED_FILE (edge list sorted by src), $SOURCE_SORTED_FILE.reverse (edge list sorted by dst), $SOURCE_SORTED_FILE.node_ids.uniq (all node ids sorted)" $GREEN

LOG "DONE"
