#!/usr/bin/env python3
#analyze_output.py
import sys
import argparse
import os, errno
import glob
import re


parser = argparse.ArgumentParser(description="Analyses the output files from autodock or vina to make a .csv file of binding modes and energies.")
parser.add_argument("pdbqt_out_list", nargs='+', type = str,
	help="""List of pdbqt files. The names of the .pdbqt files will be used to derive the names of the entries in the output table. 
	Assuming the output files follow the format <ligandID>_out.pdbqt, the entry name will be <ligandID>.""")
parser.add_argument("-r", metavar="receptor_name", dest= "receptor_name", type = str,
	help="Optional receptor name. If supplied, a column will be added to the output file containing the name of the receptor." ,
	default="")
parser.add_argument("-o", metavar="output", dest = "output",
	help="Name (without extension) of the output file.",
	default = "output.csv")
parser.add_argument("-s", metavar="sort", dest="sort", type = str,
	help="""Sort all output. Options are 'none' for no sorting, 'id' for sorting by <ligandID>, 'energy' for sorting by energy, and 'le' for sorting by 
	ligand efficiency. Default behaviour is to sort by energy. """ ,
	choices=['none', 'id', 'energy', 'le'], default="energy")
parser.add_argument("-a", "--append", 
	help="""Appends output to an existing file, which should be specified using [-o output]. Basic checking is performed to make sure there 
	are no exact matches between new lines and lines already in the file. This checking may be slow for very large files.""",
	default=False, action="store_true")
parser.add_argument("--silent",
	help="Silently overrides the overwriting of the output file (if it exists) instead of throwing an error.",
	default=False, action="store_true")
parser.add_argument("-m", metavar="modes", dest = "modes", type = int, 
	help="Number of binding modes to include in the file.",
	default=1)
parser.add_argument("-e", metavar="energy_cutoff", type=float, nargs='?', dest = "energy",
	help="""Energy cutoff for inclusion of a protein-ligand interaction to be included in the output. Setting this parameter may 
	prevent some modes specified with --modes from being written.""", 
	default=0, const=-7)
parser.add_argument("-le", metavar="le_cutoff", type=float, nargs='?', dest = "le_cutoff",
	help="""Ligand efficiency cutoff for inclusion of a protein-ligand interaction to be included in the output. Setting this parameter may prevent some 
	modes specified with --modes from being written.""", 
	default=0, const=0.6)
parser.add_argument("-f", "--include_first", 
	help="Includes the lowest-energy binding mode of an output file even if the [-e [energy_cutoff]] or [-le [le_cutoff]] values would normally exclude the mode.", 
	default=False, action="store_true")
parser.add_argument("-N", metavar="N_best", dest="nprint", type=int, nargs='?',
	help="""Prints the info for the top [N_best] ligans to stdout. 'Best' here is defined by whatever was chosen for the [--sort] method. 
	If [--sort 'none'] was chosen, then sorting method for this flag defaults to 'energy'. Using this flag without an argument prints the top 10.""",
	default=0, const=10)
parser.add_argument("-v", "--verbose", 
	help="Prints out log info to stdout.",
	default=False, action="store_true")



args = parser.parse_args()


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def get_unique(oldlines, newlines):
	#function to check old file for lines that match new data. 
    novellines = []
    for line in newlines:
        if line in oldlines:
            continue
        else:
            novellines.append(line)
    return novellines

#pdbqt_out_folder = sys.argv[1]
zinc_table = sys.argv[2]
#analysis_table = sys.argv[3]
most_recent_count = 10

pdbqt_files = args.pdbqt_out_list
analysis_table = args.output

#remove old file if it exists, or warn of potential overwrite.
if args.silent and not args.append:
	silentremove(analysis_table)
else:
	if os.path.isfile(args.output) and not args.append:
		print("The specified output file already exists, and the option to append was not selected. Please use option --help for help, or use option --silent to override. Exiting gracefully...")
		exit()


use_receptor = True #flag for using a receptor column in the output data.
if args.receptor_name=="": #turn off receptor column if current input does not have a receptor.
	use_receptor = False

update_old = False #flag for if appending data, if the old data does NOT have a receptor column and will need one. 

#possible headers from this file. Used for determining if file for appending has a receptor column.
header = "ZincID,nHeavyAtoms,Model,Energy,Efficiency\n"
alt_header = "ZincID,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency\n"	
if args.append: 
	data_out = open(analysis_table, 'r')
	oldlines = data_out.readlines()
	curr_header = oldlines[0] #first line of file to append to.
	data_out.close()
	if curr_header not in [header, alt_header]:
		print("Supplied output file does not appear to be the correct format for appending. Please see option --help for help. Exiting gracefully...")
		exit()
	if curr_header == alt_header:
		use_receptor = True #old header in file for appendin has a receptor column, so we need to use it.
	if curr_header == header and use_receptor: #old data does not have a receptor column but needs one!
		update_old = True
if args.verbose and use_receptor: print("Using receptor column.")
if args.verbose and update_old: print("Old file does not have a receptor column, and will be updated.")

#initialize values
cmpd_nmb = 0
data_lines =[]

#main loop for reading pdbqt files.
for file in pdbqt_files:
	cmpd_nmb = cmpd_nmb + 1
	zinc_id = file.split("\\")[-1].split("_")[0] 
	#find match in zinc library table
	curr_file = open(file, 'r')


	n_heavy_atoms = 0 #initiate heavy atom count
	curr_model = -1 #current model ID
	on_ligand = False #flag for whether the line belongs to a ligand (versus an amino acid residue)
	curr_energy = 0
	#start reading files
	for line in curr_file:

		linetype = line.split()[0].strip()
		if linetype=="MODEL":
			curr_model = int(line.split()[1])
		if curr_model in range(1,args.modes+1):
			if linetype=="REMARK" and line.split()[1]=="VINA":
				on_ligand = True
				curr_energy = float(line.split()[3])
				if curr_energy > args.energy+0.002:
					above_e_cutoff = True
			elif linetype=="BEGIN_RES":
				on_ligand = False
			elif on_ligand and linetype in ["ATOM","HETATM"] and curr_model==1 and "H" not in line.split()[2]:
				n_heavy_atoms = n_heavy_atoms+1
				
		else:
			break
		if linetype=="ENDMDL":
			#first model has been processed, so we can no analyze the  model appropriately. 
			eff = -(curr_energy)/float(n_heavy_atoms)

			if use_receptor:
				out_string = "{},{},{:n},{:n},{:.1f},{:.3f}\n".format(zinc_id, args.receptor_name, n_heavy_atoms, curr_model, curr_energy, eff)
			else:
				out_string = "{},{:n},{:n},{:.1f},{:.3f}\n".format(zinc_id, n_heavy_atoms, curr_model, curr_energy, eff)

			if (args.include_first and curr_model==1) or (curr_energy <= args.energy and eff >= args.le_cutoff):
				data_lines.append(out_string)

	curr_file.close()
if args.verbose: print("{} .pdbqt files read!".format(cmpd_nmb))
if args.append:
	oldlines = oldlines[1:len(oldlines)]
	if update_old: #update old info by inserting an extra comma to make a new column.
		updated_lines = []
		for line in oldlines:
			updated_lines.append(line.replace(",", ",,", 1))
		oldlines = updated_lines


	newlines = get_unique(oldlines, data_lines)
	if args.verbose: print("{} new lines of data to be added to the output file.".format(len(newlines)))
	if args.sort == "none":
		data_out = open(analysis_table, 'a')
		data_out.write("".join(data_lines))
		data_out.close()
		if args.nprint>0:
			data_lines.sort(key=lambda x:float(x.split(",")[-2]))
			print("Top ligand:") if args.nprint==1 else print("Top {:n} ligands:".format(args.nprint))
			print("".join(data_lines[0:args.nprint])[0:-1])
		exit()
	else:
		oldlines.extend(newlines)
		data_lines = oldlines

if args.sort == "energy":
	data_lines.sort(key=lambda x:float(x.split(",")[-2]))
elif args.sort == "le":
	data_lines.sort(key=lambda x:float(x.split(",")[-1]), reverse=True)
elif args.sort == "id":
	data_lines.sort(key=lambda x:(x.split(",")[0]))
if args.verbose: print("Files sorted by {}.".format(args.sort))

data_out = open(analysis_table, 'w')
data_out.write(alt_header) if use_receptor else data_out.write(header)
data_out.write("".join(data_lines))
data_out.close()

if args.nprint>0:
	if args.sort == "none":
		data_lines.sort(key=lambda x:float(x.split(",")[-2]))
	print("Top ligand:") if args.nprint==1 else print("Top {:n} ligands:".format(args.nprint))
	print("    "+"    ".join(data_lines[0:args.nprint])[0:-1])

exit()
