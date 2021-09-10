#! /bin/bash

output_folder=./output_test/
home_dir=`pwd`

#remove old directory if it exists, and make it again.
if [ -d "$output_folder" ]; then rm -Rf $output_folder; fi
mkdir ${output_folder}

#iterated over the receptors
for receptor in ./flex_receptors/*
do 
	folder_name="$(dirname $receptor)"
	folder_name="$(basename $receptor)"
	mkdir ${output_folder}/${folder_name}
	working_dir="${output_folder}/${folder_name}"
	cp $receptor/*.pdbqt $working_dir
	cp $receptor/conf.txt $working_dir
	for ligand in ./test/*.pdbqt
	#iterate over the ligands for each receptor. 
	do
		cp $ligand $working_dir
		echo Processing $ligand with receptor $folder_name.
		cd $working_dir
		b=$(basename $ligand) #Ligand file 
		c=${b%.pdbqt} #Ligand ID
		#now actually run vina
		vina --config conf.txt --ligand $b --out ${c}_out.pdbqt --log log.txt
		rm $b #remove input .pdbqt file now that we're done with it
		cd $home_dir
	done
done
