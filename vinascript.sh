#!/bin/bash 

batchsize=96 #size of batches done in parallell before checking for priorities
find  -maxdepth 1 -name '[A-Z]*0[0-9]*[0-9].pdbqt' | sort -R > list.txt  #list of current [A-Z]* ligands in the folder 
ligand_count=$(wc -l list.txt |awk '{ print $1 }') #umber of ligands in the list
priority_count=$(ls -1 ./priority/*.pdbqt | wc -l)
ligand_count=$((ligand_count+priority_count))
# echo $batchsize 
# echo $ligand_count

while [ "$ligand_count" -gt 0 ] 
do
    #counts the number of pdbqt files in the priority folder
    priority_count=$(ls -1 ./priority/*.pdbqt | wc -l)
    if [ "$priority_count" -gt 0 ] ; then 
        #if there are any, they get moved to the front of the line
        find ./priority -maxdepth 1 -name '[A-Z]*0[0-9]*[0-9].pdbqt' | sort -R > list.txt
        #followed by everything else in the normal folder...
        find  -maxdepth 1 -name '[A-Z]*0[0-9]*[0-9].pdbqt' | sort -R >> list.txt
    else
        #just the normal folder makes the whole list. 
        find  -maxdepth 1 -name '[A-Z]*0[0-9]*[0-9].pdbqt' | sort -R > list.txt
    fi

    head -n $batchsize list.txt > temp.txt #get a batchsize of ligands in a temp folder!
    ls $(cat temp.txt) -lS | awk -F" " '{print $9}' > currbatch.txt #sort these files by size
    rm temp.txt
    #make output names for the ligands in the output folder!
    cat currbatch.txt | sed 's/.*\///' | sed 's/\.[^.]*$//' | sed "s/$/_out.pdbqt/g" | sed 's/^/.\/output\//' > currbatch_out.txt

    #performs the docking in parallel in 24 processes until the batchsize is hit.
    parallel -j24  --link  vina --config DmelTyrR2_conf.txt --ligand {1} --out {2} :::: currbatch.txt :::: currbatch_out.txt 

    #moves files that were completed to the completed folder. 
    while read -r filename
    do
        movename=$(basename "$filename")
        mv -- "$filename" ./completed/"$movename"
    done < currbatch.txt

    #recount the ligands!
    ligand_count=$(wc -l list.txt |awk '{ print $1 }')
done
