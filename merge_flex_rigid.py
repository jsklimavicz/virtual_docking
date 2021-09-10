#!/usr/bin/env python3
#merge_rigid_flex.py
import sys
import argparse
import os, errno
import glob
import re
import math
import getopt

parser = argparse.ArgumentParser(description="Merges a the rigid protein portion from a pdbqt file with the flexible residue portion from autodock/vina.")
parser.add_argument("rigid", 
	help="Protein rigid file, in .pdbqt format")
parser.add_argument("vina_output", 
	help="Output file from vina with ligand and flexible residues, in .pdbqt format")
parser.add_argument("-b", metavar="binding_modes", type=str, 
	help="Binding modes of interest, with numbers seperated by '-' or ',', without spaces, e.g. 1-3,5,7. Default is all available binding modes.", 
	default='all')
parser.add_argument("-s", "--split", 
	help="Whether each binding mode is split into a different file. Default is to keep all binding modes in one file as different models.", 
	default=False, action="store_true")
parser.add_argument("-E", "--energy_out",
	help="Print the energy of each ligand/receptor model to stdout.", 
	default=False, action="store_true")
parser.add_argument("-e", metavar="energy", type=float, 
	help="Energy cutoff for inclusion of a protein-ligand interaction to be included in the output. Setting this parameter may prevent some modes specified with --modes from being written.", 
	default=0)
parser.add_argument("-o", metavar="output", type=str, 
	help="Basename of the output file. The default is to use the name of the rigid receptor before th first '_', and the name of the output file.", 
	default='')


args = parser.parse_args()


#check input args
if not os.path.isfile(args.rigid):
	print("rigid file does not exist. Please check filename or use merge_rigid_flex.py --help for options. Exiting gracefully...")
	exit()
if not os.path.isfile(args.rigid):
	print("vina_output file does not exist. Please check filename or use merge_rigid_flex.py --help for options. Exiting gracefully...")
	exit()

#read rigid file.
rigid_file = open(args.rigid, 'r')
rigid_lines = rigid_file.readlines()
rigid_file.close()

#read flex output file. 
lig_file = open(args.vina_output, 'r')
lig_lines = lig_file.readlines()
lig_file.close()

num_models = 1
for line in lig_lines:
	linetype = line.split()[0].strip()
	if linetype=="MODEL":
		num_models = num_models + 1

if args.modes == "all":
	models = range(1,num_models)
else:
	model_chunks = args.modes.split(",")
	models = []
	for x in model_chunks:
		if x.isdigit() and int(x)>0:
			models.append(int(x))
		elif "-" in x:
			pieces = x.split("-")
			if len(pieces) != 2 or not pieces[0].isdigit() or not pieces[1].isdigit():
				print("Specified binding modes were not positive integers. Please use merge_rigid_flex.py --help for options. Exiting gracefully...")
			else: 
				for i in range(int(pieces[0]), int(pieces[1])+1):
					models.append(i)
		else:
			print("Specified binding modes were not positive integers. Please use merge_rigid_flex.py --help for options. Exiting gracefully...")
			exit()




receptor_name = args.rigid.split("_")[0]
ligand_name = args.vina_output.split("_")[0]
if args.output == "":
	out_basename = receptor_name+"_"+ligand_name
else:
	out_basename = args.output


if not args.split:
	out_name = out_basename + ".pdbqt"
	out_file = open(out_name, 'w')

at_least_one_model_written = False
above_e_cutoff = False
models_written = 0
for model in models:

	curr_ligand = []
	curr_aa_list = []	
	on_model = False
	on_ligand = False
	on_aa = False
	aa_numbers = []

	#read flex output file. 
	for line in lig_lines:
		linetype = line.split()[0].strip()
		if linetype=="MODEL" and int(line.split()[1])==model:
			on_model = True
		elif linetype=="MODEL" and not int(line.split()[1])==model:
			on_model = False
		elif on_model:
			if linetype=="REMARK" and line.split()[1]=="VINA":
				on_ligand = True
				if float(line.split()[3]) > args.energy+0.002:
					above_e_cutoff = True
				if args.energy_out:
					print(ligand_name + " in " + receptor_name + " energy: " + line.split()[3])
			elif linetype=="BEGIN_RES":
				on_ligand = False
				on_aa = True
			elif on_ligand and linetype == "ATOM":
				curr_ligand.append(line)
			elif on_aa and linetype == "ATOM":
				curr_aa_list.append(line)
				aa_numbers.append(int(line.split()[4]))
	if not at_least_one_model_written:
		flex_list = []
		for flex_aa in set(aa_numbers):
			flex_list.append(curr_aa_list[aa_numbers.index(flex_aa)].split()[3] + str(flex_aa))
		print("Flexible residues are " + "_".join(flex_list))
	if args.split:
		out_name = out_basename + ".mode" + str(model) + ".pdbqt"
		out_file = open(out_name, 'w')
	out_file.write("MODEL\n")
	# print(curr_ligand)
	# print(aa_numbers)
	prev_atom_num = 0
	aa_lines_added = 0
	out_file = open(out_name, 'a')
	curr_atom_num = 0
	for line in rigid_lines:
		curr_atom_num = curr_atom_num + 1
		curr_aa_num = int(line.split()[4])
		while curr_aa_num-1 in aa_numbers:
			new_line = curr_aa_list[aa_lines_added][0:7] + '{:4d}'.format(curr_atom_num) + curr_aa_list[aa_lines_added][11:len(curr_aa_list[aa_lines_added])]
			new_line = new_line[0:20] + " A" + new_line[22:len(line)]
			out_file.write(new_line)
			aa_lines_added = aa_lines_added + 1
			curr_atom_num = curr_atom_num + 1
			aa_numbers.remove(curr_aa_num-1)
		new_line = line[0:7] + '{:4d}'.format(curr_atom_num) + line[11:len(line)]
		new_line = new_line[0:20] + " A" + new_line[22:len(line)]
		out_file.write(new_line)
		prev_atom_num = curr_atom_num

	for line in curr_ligand:
		new_line = line[0:7] + '{:4d}'.format(prev_atom_num+1) + line[11:17] + "ZIN B" + line[22:len(line)]
		prev_atom_num = prev_atom_num+1
		out_file.write(new_line)
	out_file.write("ENDMDL\n")
	if args.split:
		out_file.close()

	if above_e_cutoff:
		if at_least_one_model_written: 
			print(str(models_written) + " model(s) written. Remaining models are of higher energy than supplied cuttoff and are not written to file. See --help for details.")
			exit()
		else:
			print("The first model is higher energy than the cutoff supplied. Writing it anyway and then exiting. See --help for details.")
	at_least_one_model_written = True
	models_written = models_written + 1
	if above_e_cutoff:
		exit()

if not args.split:
	out_file.close()


