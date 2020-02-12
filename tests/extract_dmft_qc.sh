#!/bin/bash

# Script for redirecting testcode to the ed_solver_serial.log file
# rather than looking in test.out.<date>
# 
# Prints a list of keywords and values which testcode then compares
#
# Written by Edward Linscott Feb 2020

file=$1
# If 'test' in filename, redirect it to check the logfile
if [[ $file == 'test'* ]]; then
   # Find the logfile from the last iteration
   file=$(ls -d atom1/dir_green_output* | sort --numeric-sort -t 'r' -k 4 | tail -1)
   file=$file/ed_solver_serial.log
fi

# Get list of unique keywords
keys=()
IFS=$'\n'
for line in $(grep 'QC' $file);
do
   key=${line%]*}
   key=${key#*[}
   if [[ ! " ${keys[@]} " =~ " ${key} " ]]; then
      keys+=( $key )
   fi
done

# Find all instances of each keyword and print in table format
for key in ${keys[@]};
do
   header=$(echo $key | sed 's/ /_/g')
   echo $header
   grep '<QC>' $file | grep $key | awk ' { print $(NF) } '
done
