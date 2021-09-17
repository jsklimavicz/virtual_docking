#!/bin/bash



/usr/local/bin/obabel substances.smi -c -p 7.0 --minimize --gen3d --ff GAFF --better -opdb -m -O "${fileid}dummy.pdb"; 
for file2 in dummy*.pdb; do #for each dummy file...
    zinc=$(grep 'JSK\|ZINC' $file2 |  tr -s ' ' | cut -d ' ' -f2); #get the ZINC name
    i=1; #set i to 1
    if [[ -f "$zinc.pdb" ]] ; then #if the zinc.pdb file already exists, then the current file is a protomer.
        i=$((i+1)); #increment i (minimum here is i=2)
        mv "$zinc.pdbqt" "$zinc.1.pdbqt" #move the pdbqt file
        while [[ -f "$zinc.$i.pdb" ]]; do #if zinc.2+.pdb is a file already, further increment i.
            i=$((i+1));
        done
        file3="$zinc.$i.pdb"; #name for the file
    else 
        file3="$zinc.pdb"; #name for the file
    fi;
    mv $file2 $file3; #move the dummy name to  the new zinc name
    
    found=0
    if [[ -e  ../"$file3"qt ]] ; then 
        found=$((found+1))
        echo "$file3"qt found in priority folder. Doing nothing.
    fi
    if [[ -e  ../../"$file3"qt ]] ; then 
        rm ../../"$file3"qt .
        echo "$file3"qt found in ligands folder. Replacing...
    fi
    if [[ $found -lt 1 ]] ; then 
        pythonsh /opt/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l $file3; #make pdbqt file 
    fi
    
done
rm *.pdb 
cp ../../output/output.csv backup_output.csv

shopt -s nullglob
for file in [A-Z]*.pdbqt; do
    check=$(echo "${file}" | awk -F. '{print $1}')
    #echo $check
    if test -n "$(find ../../output/ -maxdepth 1 -name $check'*_out.pdbqt' -print -quit)"; then
        rm ../../output/"$check"*_out.pdbqt
    fi
    shortname=$(echo $check | grep -o -E '[A-Z]+[0]+[^ ]+' | sed 's/[A-Z]*[0]*//' | awk -F. '{print $1}')
    if test -n "$(find ../../output/images/ -maxdepth 1 -name 'ZINC'$shortname'.*' -print -quit)"; then
        rm ../../output/images/ZINC"$shortname".*
    fi
    if test -n "$(find ../../output/images/ -maxdepth 1 -name 'JSK'$shortname'.*' -print -quit)"; then
        rm ../../output/images/JSK"$shortname".*
    fi
    for csv in ../../output/*.csv ; do
        if grep -q $check $csv ; then
            grep -v $check $csv > $csv.temp  && mv $csv.temp $csv
            echo Removing $file from $csv
        fi
    done
done
shopt -u nullglob

mv *.pdbqt ..
