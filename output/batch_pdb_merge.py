#!/usr/bin/env python3
#batch_pdb_merge.py
import sys
import argparse
import os, errno
import glob
import re
import math
import getopt

parser = argparse.ArgumentParser(description="""Selects the ligands from an output file from ""analyse_output.py"" and makes a .pdb file for viewing in pymol.
	By default, the program makes a .pdb file containing only the flexible residues, and a .mol2 file containing the ligands.""")
parser.add_argument("lig_data", 
	help="Ligand output data file from ""analyse_output.py"", in .csv format")
parser.add_argument("rigid", type=str, 
	help="Path to rigid file.")
parser.add_argument("-o", metavar=["residues", "ligands"], dest="output", nargs=2, type = str,
	help="Output file basenames for the output .pdb file and .mol2 file. Defaults are ""residues"" and ""ligands"" ",
	default=["residues", "ligands"])
parser.add_argument("--merge", dest="merge", 
	help="Outputs the residues and ligands in a single pdb file.",
	default=False, action="store_true")
parser.add_argument("-L", metavar="ligand_file_type", dest="lfiletype", type = str,
	help="Output file basenames for the output .pdb file and .mol2 file. Defaults are ""residues"" and ""ligands"" ",
	choices=["sdf","mol","mol2","xyz"],
	default="mol2")
parser.add_argument("--suffix", metavar="suffix", type = str,
	help="Suffix for pdbqt output files. The default is ""_out.pdbqt"".",
	default="_out.pdbqt")
parser.add_argument("--full", 
	help="Include all residues from the rigid file. Default behaviour is to only include atoms for residues that were flexible.",
	default=False, action="store_true")
parser.add_argument("-n", type=int,
	help="Number of files to include in the output. Default=50.", 
	default=50)
parser.add_argument("--filter", metavar=('atoms', 'eff', 'energy'), dest="filters", nargs=3, type=float,
	help="""Filters the output to meet the specified filter criteria. Filters are based on the number of heavy atoms, the efficiency, and the energy. This is
	done by default with the values [8, 0.4, -7], but these numbers can be changed here.""",
	default=[8, 0.4, -7])
parser.add_argument("-C", "--carbon", metavar=('low', 'high'), dest="carbon", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of nitrogen atoms. Default is [6,30].""",
	default=[6,30])
parser.add_argument("-N", "--nitrogen", metavar=('low', 'high'), dest="nitrogen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of nitrogen atoms. Default is [0,6].""",
	default=[0,6])
parser.add_argument("-O", "--oxygen", metavar=('low', 'high'), dest="oxygen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of oxygen atoms. Default is [0,6].""",
	default=[0,6])
parser.add_argument("-F", "--fluorine", metavar=('low', 'high'), dest="fluorine", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of fluorine atoms. Default is [0,9].""",
	default=[0,9])
parser.add_argument("-S", "--sulfur", metavar=('low', 'high'), dest="sulfur", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of sulfur atoms. Default is [0,2].""",
	default=[0,2])
parser.add_argument("--halogens", metavar=('low', 'high'), dest="halogen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of halogens (F, Cl, Br, I). Default is [0,9].""",
	default=[0,9])
parser.add_argument("--Hbond", metavar=('low', 'high'), dest="hbond", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of hydrogen-bond forming atoms (N and O). Default is [2,7].""",
	default=[2,7])
parser.add_argument("-MW", metavar=('low', 'high'), dest="mw", nargs=2, type=float,
	help="""Lower and upper bounds (inclusive) for the molecular weight of the ligand. Default is [100,450].""",
	default=[100.0,450.0])
parser.add_argument("--rings", metavar=('low', 'high'), dest="rings", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of rings in the ligand. Default is [1,6].""",
	default=[1,6])
parser.add_argument("-s", metavar="sort", dest="sort", type = str,
	help="""Sort all output. Options are: \n
	'none' for no sorting,\n
	'id' for sorting by <ligandID>,\n 
	'energy' for sorting by energy, \n
	'natoms' for number of heavy atoms,\n 
	'le' for sorting by ligand efficiency. \n
	'james' for james's metric, e^2/n_atoms.\n
	Default behaviour is to sort by energy. """ ,
	choices=['none', 'id', 'energy', 'natoms', 'le', 'james'], default="energy")
parser.add_argument("-R", "--reverse",
	help="Reverses the order of normal sort.",
	default=False, action="store_true")
parser.add_argument("-v", "--verbose", 
	help="Prints out log info to stdout.",
	default=False, action="store_true")
parser.add_argument("-ID", 
	help="Changes the ID name in the pdb file.", type=str,
	default="LIG")



args = parser.parse_args()


def formula_sort(formula):
	s = re.findall('([A-Z][a-z]?)([0-9]*)',formula)
	rev_form = ""
	for elem, count in s:
		if count == "": count = "1"
		rev_form = rev_form + "{}{:3n}".format(elem, int(count))
	return rev_form

def allowed_formula(formula, args):
	s = re.findall('([A-Z][a-z]?)([0-9]*)',formula)
	rev_form = ""
	for elem, count in s:
		if count == "": count = "1"
		rev_form = rev_form + "{:s}{:n}".format(elem, int(count))
	
	carbon_count = 0
	oxygen_count = 0
	nitrogen_count = 0
	fluorine_count = 0
	halogen_count = 0
	hb_count = 0
	sulfur_count = 0

	carbons = re.search('C[0-9]*', rev_form)
	if carbons: carbon_count = int(carbons.group(0)[1:])
	if (args.carbon[0] > carbon_count or args.carbon[1] < carbon_count): return False

	oxygens = re.search('O[0-9]*', rev_form)
	if oxygens:
		oxygen_count = int(oxygens.group(0)[1:])
		hb_count = hb_count + oxygen_count
	if (args.oxygen[0] > oxygen_count or args.oxygen[1] < oxygen_count): return False

	nitrogens = re.search('N[0-9]*', rev_form)
	if nitrogens:
		nitrogen_count = int(nitrogens.group(0)[1:])
		hb_count = hb_count + nitrogen_count
	if (args.nitrogen[0] > nitrogen_count or args.nitrogen[1] < nitrogen_count): return False
	if (args.hbond[0] > hb_count or args.hbond[1] < hb_count): return False

	sulfurs = re.search('S[0-9]*', rev_form)
	if sulfurs: sulfur_count = int(sulfurs.group(0)[1:])
	if (args.sulfur[0] > sulfur_count or args.sulfur[1] < sulfur_count): return False

	fluorines = re.search('F[0-9]*', rev_form)
	if fluorines:
		fluorine_count = int(fluorines.group(0)[1:])
		halogen_count = halogen_count + fluorine_count
	if (args.fluorine[0] > fluorine_count or args.fluorine[1] < fluorine_count): return False

	chlorines = re.search('Cl[0-9]*', rev_form)
	if chlorines: halogen_count = halogen_count + int(chlorines.group(0)[2:])
	bromines = re.search('Br[0-9]*', rev_form)
	if bromines: halogen_count = halogen_count + int(bromines.group(0)[2:])
	iodines = re.search('I[0-9]*', rev_form)
	if iodines: halogen_count = halogen_count + int(iodines.group(0)[1:])
	if (args.halogen[0] > halogen_count or args.halogen[1] < halogen_count): return False

	return True

def allowed_line(formula, args, header):
	curr_form = line[header.split(",").index('Formula')]
	if not allowed_formula(curr_form, args): return False
	curr_eff = float(line[header.split(",").index('Efficiency')])
	if curr_eff < args.filters[1]: return False 
	curr_e = float(line[header.split(",").index('Energy')])
	if curr_e > args.filters[2]: return False 
	curr_natoms = int(line[header.split(",").index('nHeavyAtoms')])
	if curr_natoms < args.filters[0]: return False 
	curr_MW = float(line[header.split(",").index('MW')])
	if curr_MW < args.mw[0] or curr_MW > args.mw[1]: return False 
	curr_nrings = int(line[header.split(",").index('nRings')])
	if curr_nrings < args.rings[0] or curr_nrings > args.rings[1]: return False 
	return True

def flex_only_pdb(flex_out_file, rigid_file, model_number):
	full_out = full_pdb(flex_out_file, rigid_file, model_number)
	model = full_out[0]
	res = full_out[1]

	flex=[]
	for line in model:
		split_line=line.split()
		if split_line[0] in ["HETATM", "ATOM"]:
			if int(split_line[5]) in res or split_line[3] == args.ID: flex.append(line)
		elif split_line[0] == "CONECT": flex.append(line)
	
	return [flex, res]

def full_pdb(flex_out_file, rigid_file, model_number):
	flex = open(flex_out_file, 'r')
	lig_lines = flex.readlines()
	flex.close()

	rigid = open(rigid_file, 'r')
	rigid_lines = rigid.readlines()
	flex.close()

	curr_ligand = []
	curr_aa_list = []	
	on_model = False
	on_ligand = False
	on_aa = False
	aa_numbers = []

	output=[]

	#read flex output file. 
	for line in lig_lines:
		linetype = line.split()[0].strip()
		if linetype=="MODEL" and int(line.split()[1])==model_number:
			on_model = True
		elif linetype=="MODEL" and int(line.split()[1]) > model_number:
			on_model = False
			break
		elif on_model:
			if linetype=="REMARK" and line.split()[1]=="VINA":
				on_ligand = True
			elif linetype=="BEGIN_RES":
				on_ligand = False
				on_aa = True
			elif on_ligand and linetype in ["ATOM","HETATM"]:
				curr_ligand.append(line)
			elif on_aa and linetype == "ATOM":
				curr_aa_list.append(line)
				aa_numbers.append(int(line.split()[4]))
	flex_list = []
	aa_numbers.sort()
	aa_numbers_store = list(set(aa_numbers))
	for flex_aa in set(aa_numbers):
		flex_list.append(curr_aa_list[aa_numbers.index(flex_aa)].split()[3] + str(flex_aa))

	prev_atom_num = 0
	aa_lines_added = 0
	curr_atom_num = 0
	for line in rigid_lines:
		curr_atom_num = curr_atom_num + 1
		curr_aa_num = int(line.split()[4])
		while curr_aa_num-1 in aa_numbers:
			new_line = curr_aa_list[aa_lines_added][0:7] + '{:4d}'.format(curr_atom_num) + curr_aa_list[aa_lines_added][11:len(curr_aa_list[aa_lines_added])]
			new_line = new_line[0:20] + " A" + new_line[22:len(line)]
			output.append(new_line)
			aa_lines_added = aa_lines_added + 1
			curr_atom_num = curr_atom_num + 1
			aa_numbers.remove(curr_aa_num-1)
		new_line = line[0:7] + '{:4d}'.format(curr_atom_num) + line[11:len(line)]
		new_line = new_line[0:20] + " A" + new_line[22:len(line)]
		output.append(new_line)
		prev_atom_num = curr_atom_num


	temp_file = open(".temp.pdbqt","w")
	temp_file.write("".join(output))
	temp_file.close()
	os.system('obabel -ipdbqt .temp.pdbqt -opdb -O .temp.pdb > log.txt 2>&1')

	temp_file = open(".temp.pdb","r")
	output = temp_file.readlines()
	temp_file.close()
	for i, s in enumerate(output):
		if "CONECT" in s:
			rm_ind = i
			break

	output = output[0:rm_ind]

	os.system('rm .temp.pdb*')


	temp_file = open(".temp.pdbqt","w")
	atom_num = 0
	for line in curr_ligand:
		new_line = line[0:7] + '{:4d}'.format(atom_num+1) + line[11:17] + args.ID + " B" + line[22:len(line)]
		atom_num = atom_num+1
		temp_file.write(new_line)
	temp_file.close()

	os.system('obabel -ipdbqt .temp.pdbqt -opdb -O .temp.pdb > log.txt 2>&1')
	#os.system('obabel -ipdb .temp.pdb -h -opdb -O .temp2.pdb > log.txt')

	temp_file = open(".temp.pdb","r")
	curr_ligand = temp_file.readlines()
	temp_file.close()
	#os.system('rm .temp*.pdb*')

	for line in curr_ligand:
		if line.split()[0] in ["ATOM","HETATM"]:
			curr_num=line.split()[1]
			new_line = "HETATM " + '{:4d}'.format(int(curr_num)+prev_atom_num) + line[11:17] + args.ID + " A" + line[22:len(line)]
			output.append(new_line)
		if line.split()[0]=="CONECT":
			new_line = line[0:6]
			for i in range(1,len(line.split())):
				curr_num=line.split()[i]
				new_line = new_line + '{:5d}'.format(int(curr_num)+prev_atom_num) 
			new_line = new_line + "\n"
			output.append(new_line)

	return [output,aa_numbers_store]

if __name__ == "__main__":

	#check input args
	if os.path.isfile(args.lig_data):
		lig_file = open(args.lig_data, 'r')
		lig_lines = lig_file.readlines()
		lig_file.close()
		lig_lines = [x.strip() for x in lig_lines]
	else:
		print("ligand output file does not exist. Please check filename or use merge_rigid_flex.py --help for options. Exiting gracefully...")
		exit()

	if len(args.ID) < 3:
		print("Identifier length is not long enough. Setting to ""LIG"".")
		args.ID = "LIG"
	elif len(args.ID) >3 :
		print("Identifier length is not long enough. Truncating to {}.".format(args.ID[0:3].upper()))
		args.ID = args.ID[0:3].upper()
	if not args.ID.isalnum():
		print("Identifier length is not alphanumeric. Setting to ""LIG"".")
		args.ID = "LIG"
	if args.ID in ["ABU","ACD","ALA","ALB","ALI","ARG","ASN","ASP","ASX","BAS","CYS","CYH","CSH","CSS","CYX","GLN","GLU","GLX","GLY","HIS","HYP","ILE","LEU","LYS","MET","PCA","PHE","PRO","PR0","PRZ","SAR","SER","THR","TRP","TYR","VAL"]:
		print("Attempted to set ligand identified to the disallowed value {}. Setting to ""LIG"".".format(args.ID))
		args.ID = "LIG"



	header = lig_lines[0]
	lig_lines = lig_lines[1:len(lig_lines)]
	if len(lig_lines) < args.n:
		print(f"There are only {len(lig_lines):n} data lines in {args.lig_data}, but {args.n:n} were requested. Defaulting to {len(lig_lines):n} lines.")
		args.n = len(lig_lines)


	#special metrics
	lines_sorted = False
	sort_rev = args.reverse

	if args.sort.lower() == "james":
		lig_lines.sort(key=lambda x:float(x.split(",")[header.split(",").index('Efficiency')])*float(x.split(",")[header.split(",").index('Energy')]), reverse=args.reverse)
		sortval = "James's metric (E^2/n_heavy_atoms)"
		lines_sorted = True
	elif args.sort.lower() in ["form","formula"]:
		lig_lines.sort(key=lambda x:formula_sort(x.split(",")[header.split(",").index('Formula')]), reverse=sort_rev)
		sortval = "formula"
		lines_sorted = True
		#sort my chemical formula
	data_float = True
	if not lines_sorted:
		#header = "ZincID,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,logP,PSA,MR,SMILES\n"
		if args.sort.lower() in ["energy","e"]:
			sortval = "Energy"
		elif args.sort.lower() in ["le", "eff", "efficiency"]:
			sortval = "Efficiency"
			sort_rev = not args.reverse
		elif args.sort.lower() in ["id", "zincid", "zinc", "z"]:
			sortval = "ZincID"
			data_float = False
		elif args.sort.lower() in ["mw", "weight", "mol_weight"]:
			sortval = "MW"
		elif args.sort.lower() == "logp":
			sortval = "logP"
		elif args.sort.lower() in ["natoms", "n_atoms", "atoms"]:
			sortval = "nHeavyAtoms"
		elif args.sort.lower() in ["psa", "sa"]:
			sortval = "PS"

		# print(header)
		ind = header.split(",").index(sortval)
		if data_float:
			lig_lines.sort(key=lambda x:float(x.split(",")[ind]), reverse=sort_rev)
		else:
			lig_lines.sort(key=lambda x:(x.split(",")[ind]), reverse=sort_rev)

	reverse_string = " in reverse of the typical order" if args.reverse else ""
	if args.verbose: print("Files sorted by {}{}.".format(sortval,  reverse_string))





	model_ind = header.split(",").index("Model")

	inc_count = 0 
	count = 0

	aa_file = open(args.output[0]+".pdb",'w')
	if not args.merge: lig_file = open(args.output[1]+".pdb",'w')

	ligand_titles = []
	while inc_count < args.n and count < len(lig_lines):
		line = lig_lines[count].split(",")
		count = count + 1
		curr_id=line[0].strip()
		curr_model=int(line[model_ind])
		# print(line)
		if allowed_line(line, args, header):
			ligand_titles.append(curr_id)
			inc_count = inc_count + 1
			aa_file.write("MODEL  {}  \nTITLE    {}\n".format(inc_count, curr_id))
			if not args.merge: lig_file.write("MODEL  {}  \nTITLE    {}\n".format(inc_count, curr_id))
			curr_lines = full_pdb(curr_id+args.suffix, args.rigid, curr_model) if args.full else flex_only_pdb(curr_id+args.suffix, args.rigid, curr_model)
			aa_data = []
			if not args.merge:
				lig_data = []
				# print(curr_lines[0])
				for mod_line in curr_lines[0]:
					if mod_line.split()[0]=="ATOM": aa_data.append(mod_line)
					elif mod_line.split()[0] in ["HETATM","CONECT"]: lig_data.append(mod_line)
				lig_file.write("".join(lig_data))
				lig_file.write("ENDMDL\n\n")
			else:
				aa_data = "".join(curr_lines[0])
			aa_file.write("".join(aa_data))
			aa_file.write("ENDMDL\n\n")

	aa_file.close()

	if not args.merge: 
		lig_file.close()
		
		if args.lfiletype != "pdb": 
			if args.verbose: print(f"Converting ligands to {args.lfiletype:s} format.")
			s = f'obabel {args.output[1]+".pdb":s} -p 7.4 -o{args.lfiletype:s} -O {args.output[1]:s}.{args.lfiletype:s}  > /dev/null 2>&1'
			os.system(s)
			# os.system(f'rm {args.output[1]+".pdb":s}')

			if args.lfiletype in ['mol','mol2']:
				ligand_file = open(f'{args.output[1]:s}.{args.lfiletype:s}','r')
				ligands = ligand_file.readlines()
				ligand_file.close()
				new_ligand_file = []
				zincid_counter = 0 
				for line in ligands:
					if "ligands.pdb" in line:
						line = ligand_titles[zincid_counter] + "\n"
						zincid_counter = zincid_counter +1
					new_ligand_file.append(line)
				ligand_file = open(f'{args.output[1]:s}.{args.lfiletype:s}','w')
				ligand_file.write("".join(new_ligand_file))
				ligand_file.close()




