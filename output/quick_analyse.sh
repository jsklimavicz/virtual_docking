#!/bin/bash


echo "$(date) Commencing pdbqt output analysis." >> /tmp/analysis.log

curr_count=$(wc -l output.csv | awk '{ print $1 }')

find  . -name '[A-Z]*_out.pdbqt' -mtime -0.4 > recent_pdbqt.txt
#find  -maxdepth 1 -name 'Z*_out.pdbqt' -mtime -2.1 > $out_path/recent_pdbqt.txt
#find  -maxdepth 1 -name 'Z*_out.pdbqt' > $out_path/recent_pdbqt.txt

python3 -u analyse_output.py recent_pdbqt.txt -r DmelTyrR2 -o output -n 50 -v -s energy -L 0.4 -E -8.5 --filter --recent 24

newlines=$(cat recent_pdbqt.txt | wc -l)
totallines=$(cat output.csv | wc -l) 

new_count=$(wc -l output.csv | awk '{ print $1 }')

echo "$(date) Analysis successful. $newlines files submitted for potential analysis; $((new_count-curr_count)) new files read. Current output file is $((totallines-1)) lines." >> /tmp/analysis.log
