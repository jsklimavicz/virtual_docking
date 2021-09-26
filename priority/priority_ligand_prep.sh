#!/bin/bash
    mv nohup.out nohup.out.prev
    /usr/local/bin/obabel substances.mol2 -c -p 7.0 --minimize --ff GAFF -opdb -m -O "${fileid}dummy.pdb"; 
    for file2 in dummy*.pdb; do #for each dummy file...
        zinc=$(grep 'JSK\|ZINC' $file2 |  tr -s ' ' | cut -d ' ' -f2); #get the cmpd name
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
            mv ../../"$file3"qt .
            found=$((found+1))
            echo "$file3"qt found in ligands folder. Moving file to here.
        fi
        if [[ $found -lt 1 ]] ; then 
            pythonsh /opt/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l $file3; #make pdbqt file 
        fi
        
    done
    rm *.pdb 
    for file in [A-Z]*.pdbqt; do
        check=$(echo "${file%.pdbqt}")
        if [[ -e  ../../output/"$check"_out.pdbqt ]] ; then 
            rm $file
            echo Removing $file because it has already been tested.
        fi  
done
    mv *.pdbqt ..
