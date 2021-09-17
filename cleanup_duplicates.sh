#!/bin/bash

tested_list=$( find output -maxdepth 1 -name 'Z*_out.pdbqt' | awk -F_ '{print $1".pdbqt"}')
count=0
for file in $tested_list; do
    curr_file=$(basename $file);
    if [[ -e $curr_file ]] ; then
        rm $curr_file
        count=$((count+1))
        echo "$curr_file" removed. 
    fi
done

tested_list2=$( find output -maxdepth 1 -name 'Z*.*_out.pdbqt' | awk -F. '{print $1".pdbqt"}')
for file in $tested_list2; do
    curr_file=$(basename $file);
    #echo $curr_file
    if [[ -e $curr_file ]] ; then
        rm $curr_file
        count=$((count+1))
        echo "$curr_file" removed. 
    fi
done

echo $count duplicate files removed.

