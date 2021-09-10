#!/usr/bin/env python3
#analyze_output.py
import sys
import argparse
import os, errno
import glob
import re
import time
# import subprocess
import numpy as np
import pathlib
import time
import math
from datetime import datetime, timedelta


parser = argparse.ArgumentParser(description="""Analyses the output files from autodock or vina to make a .csv file of binding modes and energies. This program also produces
	a file called .ligands_analysed that keeps track of analyzed ligands so that they're not checked again when the --append option is used.""")
parser.add_argument("pdbqt_to_read", nargs='+', type = str,
	help="""List of pdbqt files. The names of the .pdbqt files will be used to derive the names of the entries in the output table. 
	Assuming the output files follow the format <ligandID>_out.pdbqt, the entry name will be <ligandID>. A txt file can also be supplied that contains a list of
	all of the files to be used, with one filename per line. This is assumed if the input has a .txt extension, or if the --txtfile flag is used.""")
parser.add_argument("--txtfile",
	help="Specifies that the input list is a text file containing a list of output files, with one filename per line.",
	default=False, action="store_true")
parser.add_argument("-r", metavar="receptor_name", dest= "receptor_name", type = str,
	help="Optional receptor name." ,
	default="")
parser.add_argument("-o", metavar="output_filename", dest = "output",
	help="Name (without extension) of the output csv file. Default is ""output""[.csv].",
	default = "output")
parser.add_argument("-s", "--sort", metavar="sort", dest="sort", type = str,
	help="""Sort all output. Options are: \n
	'none' for no sorting,\n
	'id' for sorting by <ligandID>,\n 
	'energy' for sorting by energy, \n
	'natoms' for number of heavy atoms,\n 
	'le' for sorting by ligand efficiency. \n
	'james' for james's metric, e^2/n_atoms.\n
	Default behaviour is to sort by energy. """ ,
	default="energy")
parser.add_argument("-R", "--reverse",
	help="Reverses the order of normal sort.",
	default=False, action="store_true")
parser.add_argument("-m", metavar="modes", dest = "modes", type = int, 
	help="Number of binding modes to include in the file. Default behavior is to include the best binding mode only.",
	default=1)
parser.add_argument("-n", metavar="N_best", dest="nprint", type=int, nargs='?',
	help="""Prints the info for the top [N_best] ligands to stdout. 'Best' here is defined by whatever was chosen for the [--sort] method. 
	If [--sort 'none'] was chosen, then sorting method for this flag defaults to 'energy'. Using this flag without an argument prints the top 10.""",
	default=0, const=10)
parser.add_argument("--filter",
	help="Produces a filtered output csv with the filename <output_filename>_filtered.csv.",
	default=False, action="store_true")
parser.add_argument("-E", "--energy", dest="energy", nargs='?', type=float,
	help="""Upper bound (inclusive) for the binding energy. Default is -7.0.""",
	default=-7, const=-7)
parser.add_argument("-A", "--n_heavy_atoms", metavar=('low', 'high'), dest="natoms", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of non-hydrogen atoms. Default is [6,35].""",
	default=[8,35])
parser.add_argument("-L", "--ligeff", dest="le", nargs='?', type=float,
	help="""Lower bound (inclusive) for the ligand efficiency. Default is 0.4.""",
	default=0.4, const=0.4)
parser.add_argument("-C", "--carbon", metavar=('low', 'high'), dest="carbon", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of nitrogen atoms. Default is [6,30].""",
	default=[8,30])
parser.add_argument("-N", "--nitrogen", metavar=('low', 'high'), dest="nitrogen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of nitrogen atoms. Default is [0,6].""",
	default=[0,6])
parser.add_argument("-O", "--oxygen", metavar=('low', 'high'), dest="oxygen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of oxygen atoms. Default is [0,6].""",
	default=[0,6])
parser.add_argument("-F", "--fluorine", metavar=('low', 'high'), dest="fluorine", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of fluorine atoms. Default is [0,7].""",
	default=[0,7])
parser.add_argument("-S", "--sulfur", metavar=('low', 'high'), dest="sulfur", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of sulfur atoms. Default is [0,3].""",
	default=[0,3])
parser.add_argument("--halogens", metavar=('low', 'high'), dest="halogen", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of halogens (F, Cl, Br, I). Default is [0,9].""",
	default=[0,9])
parser.add_argument("--Hbond", metavar=('low', 'high'), dest="hbond", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of hydrogen-bond forming atoms (N and O). Default is [2,10].""",
	default=[2,10])
parser.add_argument("-MW", metavar=('low', 'high'), dest="mw", nargs=2, type=float,
	help="""Lower and upper bounds (inclusive) for the molecular weight of the ligand. Default is [100,500].""",
	default=[100.0,500.0])
parser.add_argument("--rings", metavar=('low', 'high'), dest="rings", nargs=2, type=int,
	help="""Lower and upper bounds (inclusive) for the number of rings in the ligand. Default is [1,5].""",
	default=[1,5])
parser.add_argument("--recent", metavar=('hours'), dest="recent",  nargs='?', type=float,
	help="""For filtering, specifies how recent data should be. Value in units of hours""",
	default = "nan", const=24)
parser.add_argument("-v", "--verbose", 
	help="Prints out log info to stdout.",
	default=False, action="store_true")
parser.add_argument("--nosave", 
	help="""No output file is saved, despite any options chosen for output. Filter file is still made if 
	the --filter option is chosen, and name is still based on output name.""",
	default=False, action="store_true")

global args
global header
global date_format
global curr_datetime
args = parser.parse_args()
header = "ZincID,Time,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,nRotors,logP,PSA,MR,SMILES\n"
date_format = "%H:%M:%S %d %b %Y"
curr_datetime = datetime.now()


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def get_unique(oldlines, newlines):
	#function to check old file for lines that match new data. (lines in new lines not in old lines)
    novellines = []
    for line in newlines:
        if line in oldlines:
            continue
        else:
            novellines.append(line)
    return novellines


class PDBQTModelInterpreter:	
	def __init__(self, file_path, model = 1):
		self.filepath = file_path
		self.model = model 
		self.model_array = []
		self.line_dict = {}
		self.get_model()
		self.calculate_parameters()
		# print(self.line_dict)

	def get_model(self):
		self.line_dict["ZincID"] = self.filepath.split("\\")[-1].split("_")[0] 

			#find match in zinc library table
		curr_file = open(self.filepath, 'r')

		self.line_dict["nHeavyAtoms"] = 0 #initiate heavy atom count
		self.line_dict["Model"] = -1 #current model ID
		on_ligand = False #flag for whether the line belongs to a ligand (versus an amino acid residue)
		self.line_dict["Energy"] = 0
		self.line_dict["nRotors"] = 0
		#start reading files
		self.model_array=[]
		for line in curr_file:
			linetype = line.split()[0].strip()
			if linetype=="MODEL": curr_model = int(line.split()[1]) 
			if curr_model==self.model:
				if linetype=="REMARK" and line.split()[1]=="VINA": self.line_dict["Energy"] = float(line.split()[3])
				elif linetype in ["ATOM","HETATM"]:
					self.model_array.append(line)
					if "H" not in line.split()[2]: self.line_dict["nHeavyAtoms"] += 1
				elif linetype in ["REMARK"] and "A    between atoms" in line:
					self.line_dict["nRotors"] += 1
				elif linetype=="BEGIN_RES":
					self.line_dict["Model"] = curr_model
					#first model has been processed, so we can no analyze the  model appropriately. 
					self.line_dict["Efficiency"] = round(-(self.line_dict["Energy"])/float(self.line_dict["nHeavyAtoms"]),4)
					break
		curr_file.close()

	def calculate_parameters(self):
		temp_mol_file=open("./.temp_mol.pdbqt","w")
		temp_mol_file.write("".join(self.model_array))
		temp_mol_file.close()
		#convert to pdb and then read
		os.system('obabel -ipdbqt ./.temp_mol.pdbqt -opdb -O ./.temp_mol.pdb > /dev/null 2>&1') #does not add hydrogens properly even with -h switch
		os.system('obabel -ipdb ./.temp_mol.pdb -h -opdb -O ./.temp_mol2.pdb > /dev/null 2>&1') #does add hydrogens properly
		#now we have to fix charges in the pdb file because .pdbqt to pdb sucksssss
		#fuucK thisss
		temp_mol_file = open("./.temp_mol2.pdb","r")
		temp_mol = temp_mol_file.readlines()
		temp_mol_file.close()
		#take care of nitrogens
		N_array = []
		flagged_N = []
		O_array = []
		flagged_O = []
		temp_mol2 = []
		for line in temp_mol:
			linetype = line.split()[0].strip()
			if linetype in ["ATOM","HETATM"]:
				if "N" in line.split()[2]: 
					N_array.append(line.split()[1].strip()) #add index of atom
				elif "O" in line.split()[2]: 
					O_array.append(line.split()[1].strip()) #add index of atom
			if linetype in ["CONECT"]: #all atoms already read in.
				if line.split()[1].strip() in N_array and len(line.split()) == 6:
					flagged_N.append(line.split()[1].strip())
				if line.split()[1].strip() in flagged_N and len(line.split()) % 2==1 :
					flagged_N.remove(line.split()[1].strip())
				if line.split()[1].strip() in O_array and len(line.split()) == 3:
					flagged_N.append(line.split()[1].strip())
		for line in temp_mol:
			if line.split()[0].strip() in ["ATOM","HETATM"]: 
				if line.split()[1].strip() in flagged_N:
					line = line.strip()+"1+\n"
				elif line.split()[1].strip() in flagged_N:
					line = line.strip()+"1-\n"
				temp_mol2.append(line)
			else:
				temp_mol2.append(line)
		temp_mol_file = open("./.temp_mol2.pdb","w")
		temp_mol_file.write("".join(temp_mol2))
		temp_mol_file.close()

		os.system('obprop ./.temp_mol2.pdb > ./.temp_mol.txt 2>&1')


		temp_mol_out=open("./.temp_mol.txt","r")
		temp_mol_data=temp_mol_out.readlines()
		temp_mol_out.close()
		# print(f"Getting OB info for {self.filepath:s}")
		for line in temp_mol_data:
			if line in ["","\n"]: continue
			if line.split()[0]=="formula":
				self.line_dict["Formula"]=line.split()[1]
			elif line.split()[0]=="canonical_SMILES":
				self.line_dict["SMILES"]=line.split()[1]
			elif line.split()[0]=="num_atoms":
				self.line_dict["nAtoms"]=int(line.split()[1])
			elif line.split()[0]=="num_bonds":
				self.line_dict["nBonds"]=int(line.split()[1])
			elif line.split()[0]=="num_rings":
				self.line_dict["nRings"]=int(line.split()[1])
			elif line.split()[0]=="mol_weight":
				self.line_dict["MW"]=float(line.split()[1])
			elif line.split()[0]=="logP":
				self.line_dict["logP"]=float(line.split()[1])
			elif line.split()[0]=="PSA":
				self.line_dict["PSA"]=float(line.split()[1])
			elif line.split()[0]=="MR":
				self.line_dict["MR"]=float(line.split()[1])


class PDBQTDataset:
	
	class PDBQTLine:

		def __init__(self, dictionary = {}):
			self.data = dictionary
			# print(self.data)
			self.data["filename"] = self.data["ZincID"] + "_out.pdbqt"
			if "ReceptorID" not in self.data: self.data["ReceptorID"] = args.receptor_name
			if "nRings" not in self.data: 
				mol_info_finder = PDBQTModelInterpreter(self.data["filename"], self.data["Model"])
				self.data.update(mol_info_finder.line_dict)
			if "Time" not in self.data or self.data["Time"] == "": self.missing_time()
			if "nRotors" not in self.data or self.data["nRotors"] == "": self.missing_rotors()


		def missing_time(self):
			fname = pathlib.Path(self.data["filename"])
			self.data["Time"] = datetime.fromtimestamp(fname.stat().st_mtime).strftime(date_format)

		def missing_rotors(self):
			self.data["nRotors"] = 0
			curr_file = open(self.data["filename"], 'r')
			for line in curr_file:
				linetype = line.split()[0].strip()
				if linetype=="MODEL": curr_model = int(line.split()[1]) 
				if curr_model==self.data["Model"]:
					if linetype in ["REMARK"] and "A    between atoms" in line:
						self.data["nRotors"] += 1
					elif linetype=="BEGIN_RES":
						curr_file.close()
						return
			curr_file.close()
			return

		def __str__(self):
			#header = "ZincID,Time,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,logP,PSA,MR,SMILES\n"
			string = ""
			# print(self.data)
			for item in header.strip().split(","): 
				if item in ["Model", "nHeavyAtoms", "nRings", "nRotors"]:
					string += f"{self.data[item]:n}" + ","
				elif item in ["Energy"]:
					string += f"{self.data[item]:.1f}" + ","
				elif item in ["Efficiency", "logP"]:
					string += f"{self.data[item]:.3f}" + ","
				elif item in ["MW", "PSA", "MR"]:
					string += f"{self.data[item]:.2f}" + ","
				else:
					string += str(self.data[item]) + ","
			return string[:-1] # + "\n"

		def print_all(self):
			string = ""
			for item in self.data: string += str(self.data[item]) + ","
			return string[:-1] # + "\n"

		def summary(self):
			string = f"{self.data['ZincID']:18s}, {self.data['Time']:20s}, {self.data['Energy']:5.1f}, {self.data['nHeavyAtoms']:2n}, {self.data['Efficiency']:5.3f}, {self.data['MW']:6.2f}, {self.data['Formula']:s}, {self.data['SMILES']:s}"
			return string

		def allowed_formula(self):
			s = re.findall('([A-Z][a-z]?)([0-9]*)',self.data["Formula"])
			self.data["revForm"] = ""
			for elem, count in s:
				if count == "": count = "1"
				self.data["revForm"] = self.data["revForm"] + "{:s}{:n}".format(elem, int(count))
			carbon_count = 0
			oxygen_count = 0
			nitrogen_count = 0
			fluorine_count = 0
			halogen_count = 0
			hb_count = 0
			sulfur_count = 0

			carbons = re.search('C[0-9]*', self.data["revForm"])
			if carbons: carbon_count = int(carbons.group(0)[1:])
			if (args.carbon[0] > carbon_count or args.carbon[1] < carbon_count): return False

			oxygens = re.search('O[0-9]*', self.data["revForm"])
			if oxygens:
				oxygen_count = int(oxygens.group(0)[1:])
				hb_count = hb_count + oxygen_count
			if (args.oxygen[0] > oxygen_count or args.oxygen[1] < oxygen_count): return False

			nitrogens = re.search('N[0-9]*', self.data["revForm"])
			if nitrogens:
				nitrogen_count = int(nitrogens.group(0)[1:])
				hb_count = hb_count + nitrogen_count
			if (args.nitrogen[0] > nitrogen_count or args.nitrogen[1] < nitrogen_count): return False
			if (args.hbond[0] > hb_count or args.hbond[1] < hb_count): return False

			sulfurs = re.search('S[0-9]*', self.data["revForm"])
			if sulfurs: sulfur_count = int(sulfurs.group(0)[1:])
			if (args.sulfur[0] > sulfur_count or args.sulfur[1] < sulfur_count): return False

			fluorines = re.search('F[0-9]*', self.data["revForm"])
			if fluorines:
				fluorine_count = int(fluorines.group(0)[1:])
				halogen_count = halogen_count + fluorine_count
			if (args.fluorine[0] > fluorine_count or args.fluorine[1] < fluorine_count): return False
			#remaining halogens
			chlorines = re.search('Cl[0-9]*', self.data["revForm"])
			if chlorines: halogen_count = halogen_count + int(chlorines.group(0)[2:])
			bromines = re.search('Br[0-9]*', self.data["revForm"])
			if bromines: halogen_count = halogen_count + int(bromines.group(0)[2:])
			iodines = re.search('I[0-9]*', self.data["revForm"])
			if iodines: halogen_count = halogen_count + int(iodines.group(0)[1:])
			if (args.halogen[0] > halogen_count or args.halogen[1] < halogen_count): return False

			return True

		def allowed_line(self, time = False):
			if time and datetime.strptime(self.data['Time'], date_format) < curr_datetime - timedelta(hours=args.recent): return False
			if not self.allowed_formula(): return False
			if self.data["Efficiency"] < args.le: return False 
			if self.data["Energy"] > args.energy: return False 
			if self.data["nHeavyAtoms"] < args.natoms[0] or self.data["nHeavyAtoms"] > args.natoms[1]: return False 
			if self.data["MW"] < args.mw[0] or self.data["MW"] > args.mw[1]: return False 
			if self.data["nRings"] < args.rings[0] or self.data["nRings"] > args.rings[1]: return False 
			
			return True

		def set_filename(self):
			self.data["filename"] = self.data["ZincID"] + "_out.pdbqt"

	def __init__(self):
		self.data_lines = []
		self.new_files_read = 0
		self.file_list = []


	def read_old_data(self, old_data_lines):
		self.old_header = old_data_lines[0].strip()
		split_header =  [x.strip() for x in self.old_header.split(",")] 
		for line in old_data_lines[1:]:
			split_line = line.strip().split(",")
			split_line = [x.strip() for x in split_line]
			ind = 0
			curr_dict = {}
			for ind in range(0,len(split_line)): 
				try : curr_dict[split_header[ind]] = float(split_line[ind])
				except: curr_dict[split_header[ind]] = split_line[ind]
			curr_line = self.PDBQTLine(curr_dict)
			self.data_lines.append(curr_line)
			self.file_list.append(curr_line.data["filename"])

	def read_new_data(self, new_file_list):
		tic = time.perf_counter()
		for file in new_file_list:
			self.new_files_read += 1
			toc = time.perf_counter()
			if args.verbose: print(f"On file {self.new_files_read:n}. {toc - tic:0.4f} seconds have elapsed. ({((toc - tic)/self.new_files_read):0.4f} s/file)",end='\r')
			mol_info_finder = PDBQTModelInterpreter(file, 1)
			curr_line = self.PDBQTLine(mol_info_finder.line_dict)
			self.data_lines.append(curr_line)

	def __iter__(self):
		for x in self.data_lines:
			yield x

	def sort_lines(self):
		#special metrics
		tic = time.perf_counter()
		lines_sorted = False
		sort_rev = args.reverse

		if args.sort.lower() == "james":
			self.data_lines.sort(key=lambda x:x.data['Efficiency']*x.data['Energy'], reverse=args.reverse)
			lines_sorted = True
			sortval = "James's metric"
		elif args.sort.lower() in ["form","formula"]:
			self.data_lines.sort(key=lambda x:self.formula_sort(x.data['Formula']), reverse=sort_rev)
			lines_sorted = True
			sortval = "formula"
		elif args.sort.lower() in ["time","date", "day", "d", "t", "datetime"]:
			self.data_lines.sort(key=lambda x:datetime.strptime(x.data['Time'], date_format), reverse = not sort_rev)
			lines_sorted = True
			sortval = "file creation date"
			#sort my chemical formula
		if not lines_sorted:
			#header = "ZincID,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,logP,PSA,MR,SMILES\n"
			if args.sort.lower() in ["energy","e"]:
				sortval = "Energy"
			elif args.sort.lower() in ["le", "eff", "efficiency"]:
				sortval = "Efficiency"
				sort_rev = not args.reverse
			elif args.sort.lower() in ["id", "zincid", "zinc", "z", "name"]:
				sortval = "ZincID"
			elif args.sort.lower() in ["mw", "weight", "mol_weight"]:
				sortval = "MW"
			elif args.sort.lower() == "logp":
				sortval = "LogP"
			elif args.sort.lower() in ["natoms", "n_atoms", "atoms"]:
				sortval = "nHeavyAtoms"
			elif args.sort.lower() in ["psa", "sa"]:
				sortval = "PSA"

			self.data_lines.sort(key=lambda x:x.data[sortval], reverse=sort_rev) 

		reverse_string = " in reverse of the typical order" if args.reverse else ""
		toc = time.perf_counter()
		if args.verbose: print(f"Files sorted by {sortval:s}{reverse_string:s}. Sorting performed in {toc - tic:0.4f} seconds")

	def formula_sort(self, formula):
		s = re.findall('([A-Z][a-z]?)([0-9]*)',formula)
		rev_form = ""
		for elem, count in s:
			if count == "": count = "1"
			rev_form = rev_form + "{}{:3n}".format(elem, int(count))
		return rev_form

	def filter(self, time = False, prev_array = []):
		self.filtered_data = []
		if prev_array == []: 
			for line in self.data_lines: 
				if line.allowed_line(time): self.filtered_data.append(line)
		else:
			for line in prev_array : 
				if line.allowed_line(time): self.filtered_data.append(line)


class DataProcessor:
	def __init__(self):
		#add .csv file designater to output file if needed
		if (not args.output[0][-4:]==".csv"):
			self.out_filename = (args.output + ".csv") 
		else:
			self.output_filename = args.output
		self.data = PDBQTDataset()
		self.filtered = False

		if os.path.exists(self.out_filename): 
			data_out = open(self.out_filename, 'r')
			self.oldlines = data_out.readlines()
			self.data.read_old_data(self.oldlines)
			data_out.close()


	def get_new_file_list(self):
		if args.txtfile or args.pdbqt_to_read[0][-4:]==".txt": self.find_new_cmpd_list()
		new_list = np.setdiff1d(args.pdbqt_to_read,self.data.file_list) # yields the elements in `pdbqt_files` that are NOT in `alreay_read`
		pdbqt_old_ind = len(args.pdbqt_to_read)
		self.files_to_read=new_list
		

		if args.verbose: 
			print("{:n} files listed for analysis.".format(pdbqt_old_ind))
			print("{:n} files found in the old file and removed from the input list. {:n} files to be read.".format(pdbqt_old_ind-len(self.files_to_read), len(self.files_to_read)))
		else:
			self.alreay_read=[]

	def find_new_cmpd_list(self):
		# print(args.pdbqt_to_read)
		pdbqt_txt_file = open(args.pdbqt_to_read[0], 'r')
		args.pdbqt_to_read = pdbqt_txt_file.readlines()
		pdbqt_txt_file.close
		args.pdbqt_to_read.sort()
		args.pdbqt_to_read = [x.strip() for x in args.pdbqt_to_read] #remove any extra spaces
		if "./" in args.pdbqt_to_read[0]: args.pdbqt_to_read = [x.split("/")[1] for x in args.pdbqt_to_read] #remove any extra spaces

	def read_new_data(self):
		self.data.read_new_data(self.files_to_read)

	def save(self):
		file_name = args.output+ ".csv"
		self.write_helper(file_name, self.data)

	def filter(self):
		self.data.filter()
		filter_table_file = args.output + "_filtered.csv"
		self.write_helper(filter_table_file, self.data.filtered_data)
		self.filter = True

	def write_helper(self, filename, PDBQTDataset_pointer):
		output = open(filename, 'w')
		output.write(header)
		for line in PDBQTDataset_pointer: 
			output.write(str(line) + "\n")
		output.close()

	def print_to_stdout(self):
		if not math.isnan(args.recent): 
			self.data.filter(time = True, prev_array = self.data.filtered_data)
		elif not self.filtered: 
			self.data.filter(time = False)

		num_print = min(args.nprint, len(self.data.filtered_data))
		if num_print == 0:
			print("No ligands meet filtering requirements. Try running again with looser print requirements.")
		else: 
			print("Top ligand:") if args.nprint==1 else print("Top {:n} ligands:".format(num_print))
			count = 0
			for line in self.data.filtered_data:
				if count < num_print: print(line.summary())
				count += 1


if __name__ == '__main__':
	tic1 = time.perf_counter()
	if args.verbose: print(f"Beginning the process of reading the previous file...")
	tic = time.perf_counter()
	datawriter = DataProcessor()
	toc = time.perf_counter()
	if args.verbose: print(f"Read previous file in {toc - tic:0.4f} seconds")

	tic = time.perf_counter()
	datawriter.get_new_file_list()
	toc = time.perf_counter()
	if args.verbose: print(f"New list of files found in {toc - tic:0.4f} seconds")

	if args.verbose: print(f"Reading new ligand files...")
	tic = time.perf_counter()
	datawriter.read_new_data()
	toc = time.perf_counter()
	if args.verbose: print(f"New files read and analyzed in {toc - tic:0.4f} seconds")
	if args.verbose: print("{} .pdbqt files read!".format(datawriter.data.new_files_read))
	datawriter.data.sort_lines()
	if not args.nosave: datawriter.save()
	tic = time.perf_counter()
	if args.verbose: print(f"Filtering data...")
	if args.filter: datawriter.filter()
	toc = time.perf_counter()
	if args.verbose: print(f"Data filtered and written to file in {toc - tic:0.4f} seconds")
	if args.nprint>0: datawriter.print_to_stdout()
	toc = time.perf_counter()
	toc1=time.perf_counter()
	if args.verbose: print(f"Total program runtime was {toc1 - tic1:0.4f} seconds")



	exit()
