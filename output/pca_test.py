#!/usr/bin/env python3
#pca_param_prep.py
import sys
import argparse
import os, errno
import glob
import re
import time
import subprocess
import numpy as np
import math
import pandas as pd
from sklearn import mixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

parser = argparse.ArgumentParser(description="""Produces adequate parameters for PCA analysis""")
parser.add_argument("lig_data", 
	help="Ligand output data file from ""analyse_output.py"", in .csv format")
parser.add_argument("-o", metavar="output_filename", dest = "output",
	help="Name (without extension) of the output csv file with extension. Default is bayesMM_ligands.csv.",
	default = "bayesMM_ligands.csv")
parser.add_argument("-s", metavar="suffix", dest = "suffix",
	help="Suffix of the output pdbqt files. Default is ""_out.pdbqt"".",
	default = "_out.pdbqt")
parser.add_argument("-n", metavar="number", dest = "n", type = int,
	help="Number of lines to include in the output file. Default is to print all.",
	default = -1)
parser.add_argument("-i", metavar="iter", dest = "iter", type = int,
	help="Number of initializations to peform for the Bayesian Gaussian Mixture Model. Default is 100.",
	default = 100)
parser.add_argument("-v", "--verbose", 
	help="Prints out log info to stdout.",
	default=False, action="store_true")
parser.add_argument("-V", "--var",  dest = "var", type = float,
	help="Ratio of variance explained by PCA.",
	default = 0.90)


args = parser.parse_args()

class Pdbqt_Parser:
	def __init__(self, zinc_id, suffix, model):
		self.zinc_id = zinc_id
		self.__path = zinc_id + suffix 
		self.__model = int(model)
		# self.ligand 
		self.residues = Residue_List()
		self.__get_ligand_and_residues()

	def __get_ligand_and_residues(self):
		pdbqt_out_file = open(self.__path, "r")
		pdbqt_lines = pdbqt_out_file.readlines()
		pdbqt_out_file.close()

		curr_ligand = []
		curr_aa_list = []	
		on_model = False
		on_ligand = False
		on_aa = False
		aa_numbers = []

		#read flex output file. 
		for line in pdbqt_lines:
			linetype = line.split()[0].strip()
			if linetype=="MODEL" and int(line.split()[1])==self.__model:
				on_model = True
			elif linetype=="MODEL" and int(line.split()[1]) > self.__model:
				on_model = False
				break
			elif on_model:
				if linetype=="REMARK" and line.split()[1]=="VINA":
					on_ligand = True
				elif linetype=="BEGIN_RES":
					on_ligand = False
					on_aa = True
				elif on_ligand and linetype in ["ATOM","HETATM"] :
					curr_ligand.append(line)
				elif on_aa and linetype in ["ATOM","HETATM"]:
					curr_aa_list.append(line)
				elif on_aa and linetype == "END_RES":
					res = Residue(curr_aa_list)
					curr_aa_list = []
					self.residues.add_res(res)
		self.ligand = Ligand(curr_ligand, zinc_id)
		# print(str(self.ligand))
		# for res in self.residues:
		# 	print(str(res.name) + str(res.num) + "  " + str(res.center) + "  "  )
		# 	print(res.min_polar_dist(self.ligand))

	def __str__(self):
		return self.__path

class Residue_List:
	def __init__(self):
		self.res_list = []

	def add_res(self, residue):
		self.res_list.append(residue)

	def __iter__(self):
		for x in self.res_list:
			yield x

	def get_res(self, res_num):
		for res in self.res_list:
			if res.num == int(res_num):
				return res
		return None


class Residue:
	def __init__(self, atom_lines):
		self.atom_list = []
		for line in atom_lines:
			self.atom_list.append(Atom(line))
		self.num = self.atom_list[0].res_seq
		self.name = self.atom_list[0].res
		self.__center()
		self.__H_bonds()

	def __add__(self, atom_line):
		self.atom_list.append(Atom(atom_line))

	def __str__(self):
		str_out = []
		for atom in self.atom_list:
			str_out.append(str(atom))
		return "\n".join(str_out)

	def __iter__(self):
		for x in self.atom_list:
			yield x

	def __center(self):
		center_list = []
		if self.name in ["TYR"]: 
			for atom in self.atom_list:
				if atom.name in ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"] and atom.element != "H":
					center_list.append(atom)
		else: 
			for atom in self.atom_list:
				if "CA" not in atom.name and "CB" not in atom.name and atom.element != "H":
					center_list.append(atom)
		center = Coord(0,0,0)
		for atom in center_list: center += atom.coord
		self.center = center/len(center_list)

#[ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",PHE","PRO","SER","THR","TRP","TYR","VAL"]:
	def __H_bonds(self):
		self.H_bonds=[]
		if self.name in ["ARG", "ASN", "ASP", "GLN", "CYS", "MET" "GLU", "HIS", "LYS", "SER", "THR", "TRP", "TYR"]:
			h_accept = ["OH", "NE", "NG", "NE1", "NE2", "ND", "ND1", "ND2", "OG", "OG1" "OD","OE","OE1", "OE2", "OD1", "OD2", "NH", "NH1", "NH2", "SG", "SD"]
			for atom in self.atom_list:
				if atom.name in h_accept:
					self.H_bonds.append(atom)

	def min_polar_dist(self, ligand_list):
		dist = math.inf
		for atom1 in ligand_list:
			if atom1.element != "H":
				for atom2 in self.H_bonds:
					if atom1.dist(atom2) < dist:
						dist = atom1.dist(atom2)
		return dist

class Ligand:
	def __init__(self, atom_lines, name):
		self.atom_list = []
		self.PolarAtoms = []
		for line in atom_lines:
			self.atom_list.append(Atom(line))
			if Atom(line).element in ["O","N"]: self.PolarAtoms.append(Atom(line))
		self.name = name
		try:
			self.__center()
		except:
			print("Error!")
			print(f"{name:s}")

	def __add__(self, atom_line):
		self.atom_list.append(Atom(atom_line))

	def __str__(self):
		str_out = []
		for atom in self.atom_list:
			str_out.append(str(atom))
		return "\n".join(str_out)

	def __center(self):
		center = Coord(0,0,0)
		for atom in self.atom_list: 
			if atom.element != "H":	center += atom.coord
		self.center = center/len(self.atom_list)

	def min_dist(self, coord):
		dist = math.inf
		for atom in self.atom_list:
			if atom.element != "H":
				if coord.dist(atom.coord) < dist:
					dist = coord.dist(atom.coord)

		return dist

	def furthest_h_bonders(self):
		
		if len(self.PolarAtoms) == 0: 
			return [self.center, self.center]
		elif len(self.PolarAtoms) == 1: 
			return [self.center.min(self.PolarAtoms[0].coord), 
				self.center.max(self.PolarAtoms[0].coord)]
		# elif len(self.PolarAtoms) == 2: 
		# 	return [self.PolarAtoms[1].coord.min(self.PolarAtoms[0].coord), 
		# 		self.PolarAtoms[1].coord.max(self.PolarAtoms[0].coord)]
		else: 
			dist = []
			for ind1 in range(0,len(self.PolarAtoms)):
				for ind2 in range(0,len(self.PolarAtoms)):
					if ind1 == ind2:
						polar_dist = self.PolarAtoms[ind1].dist_coord(self.center)
					else:
						polar_dist = self.PolarAtoms[ind1].dist(self.PolarAtoms[ind2])
					dist.append(polar_dist)
			max_value = max(dist)
			max_index = dist.index(max_value)
			if max_index // len(self.PolarAtoms) == max_index % len(self.PolarAtoms):
				return [self.center.min(self.PolarAtoms[max_index // len(self.PolarAtoms)].coord), 
					self.center.max(self.PolarAtoms[max_index // len(self.PolarAtoms)].coord)]
			else:
				atom1 = self.PolarAtoms[max_index // len(self.PolarAtoms)]
				atom2 = self.PolarAtoms[max_index % len(self.PolarAtoms)]
				return [atom1.coord.min(atom2.coord), atom1.coord.max(atom2.coord)]

	def __iter__(self):
		for x in self.atom_list:
			yield x

class Atom:
	def __init__(self, pdbqt_line):
		self.name = pdbqt_line[12:16].strip()
		self.ser = int(pdbqt_line[6:11].strip())
		self.res = pdbqt_line[17:20].strip()
		self.res_seq = int(pdbqt_line[22:26].strip())
		self.coord = Coord(float(pdbqt_line[30:38]), float(pdbqt_line[38:46]), float(pdbqt_line[46:54]))
		self.id = pdbqt_line[12:17].strip()
		if "Cl" in self.id or "CL" in self.id:
			self.element = "Cl"
		elif "Br" in self.id or "BR" in self.id:
			self.element = "Br"
		elif "O" in self.id:
			self.element = "O"
		elif "N" in self.id:
			self.element = "N"
		elif "C" in self.id:
			self.element = "C"
		elif "P" in self.id:
			self.element = "P"
		elif "S" in self.id:
			self.element = "S"
		elif "F" in self.id:
			self.element = "F"
		elif "I" in self.id:
			self.element = "I"
		elif "H" in self.id:
			self.element = "H"
		else:
			print(f"Unknown atom type. {pdbqt_line:s}")
			self.element = "H"

	def __str__(self):
		return f"{self.element:s} {self.ser:n} {self.res:s} {self.res_seq:n} {self.name:s} {str(self.coord):s}"

	def dist(self, other):
		return self.coord.dist(other.coord)

	def dist_coord(self, coord):
		return self.coord.dist(coord)

class Coord:
	def __init__(self, x=0, y=0, z=0):
		self.x=x
		self.y=y
		self.z=z

	def __add__(self, other):
		x = self.x + other.x
		y = self.y + other.y
		z = self.z + other.z
		return Coord(x,y,z)

	def __iadd__(self, other):
		x = self.x + other.x
		y = self.y + other.y
		z = self.z + other.z
		return Coord(x,y,z)

	def __sub__(self, other): 
		x = self.x - other.x
		y = self.y - other.y
		z = self.z - other.z
		return Coord(x,y,z)

	def __truediv__(self, float_in): 
		x = self.x/float_in
		y = self.y/float_in
		z = self.z/float_in
		return Coord(x,y,z)

	def dist(self, other):
		dx = self.x - other.x
		dy = self.y - other.y
		dz = self.z - other.z
		return math.sqrt(dx**2+dy**2+dz**2)

	def __str__(self):
		return f"{self.x:.3f},{self.y:.3f},{self.z:.3f}"

	def min(self, other):
		if self.z < other.z: return self
		else: return other

	def max(self, other):
		if self.z > other.z: return self
		else: return other


out_data_file = open(args.lig_data, "r")
out_data = out_data_file.readlines()
out_data_file.close()

header = out_data[0] #remove header
out_data = out_data[1:]
#find various useful indices from header
#header = "ZincID,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,logP,PSA,MR,SMILES\n"
zinc_id_ind = header.split(",").index("ZincID")
model_ind = header.split(",").index("Model")


allowed_center = [166,170,189,193,196,277,538,542,545,573]
allowed_polar= [166,170,189,193,196,277,542,573]


new_header_labels= "lig_x,lig_y,lig_z"
for res in allowed_center:
	new_header_labels += ",c" + str(res)
for res in allowed_polar:
	new_header_labels += ",p" + str(res)

polor_coord = ["polar1x","polar1y","polar1z","polar2x","polary2y","polar2z"]


new_header = header.strip() + "," + new_header_labels +  "," + ",".join(polor_coord) +"\n"

pca_data = []

n = 0
for line in out_data:
	n = n + 1
	if args.verbose and n % 1000 == 0:
		print(f"On ligand {n:n} of {len(out_data):n}")
	if n > args.n and args.n > 0:
		break
	split_line = line.split(",")
	split_line = [x.strip() for x in split_line]
	zinc_id = split_line[zinc_id_ind]
	model_id = split_line[model_ind]

	curr_model = Pdbqt_Parser(zinc_id, args.suffix, model_id)

	lig_x = curr_model.ligand.center.x
	lig_y = curr_model.ligand.center.y
	lig_z = curr_model.ligand.center.z

	add_info = f"{lig_x:.3f},{lig_y:.3f},{lig_z:.3f}"

	#print("center")
	for res in curr_model.residues:
		if res.num in allowed_center:
			dist = curr_model.ligand.min_dist(res.center)
			add_info = f"{add_info:s},{dist:.3f}"

	#print("polar")
	for res in curr_model.residues:
		if res.num in allowed_polar:
			dist = res.min_polar_dist(curr_model.ligand)
			add_info = f"{add_info:s},{dist:.3f}"

	sites = curr_model.ligand.furthest_h_bonders()

	add_info = add_info + "," + str(sites[0]) + "," + str(sites[1])


	pca_data.append(",".join(split_line) + "," + add_info + "\n")


file_out = open("temp.csv", "w")
file_out.write(new_header)
file_out.write("".join(pca_data))
file_out.close()

df=pd.read_csv("temp.csv", sep=',')

print(len(df.columns))

df1 = df.iloc[:,np.r_[2,4,7,9,10,13:len(df.columns)]]

if args.verbose: print(f"Starting PCA.")

pipeline = make_pipeline(StandardScaler(), PCA(n_components=args.var, copy=True, whiten=True, svd_solver='full', tol=0.0, iterated_power='auto', random_state=None))

X = pipeline.fit_transform(df1.to_numpy())


if args.verbose:
	print(f"Number of PCA components: {pipeline['pca'].n_components_:n}")
	print(f"Explained variance ratios: {' '.join('%.4f' % x for x in pipeline['pca'].explained_variance_ratio_):s}")

	print("Components:")
	components = ["Headers"]
	for i in range(0,pipeline['pca'].n_components_):
		components.append(f" PC{i:n}")
	print("".join('%8s' % x for x in components))



	headers_used = df1.columns.tolist()

	for feature in range(0, pipeline['pca'].n_features_):
		data_row = headers_used[feature].ljust(8)
		data_row = data_row[0:8]
		for row in pipeline['pca'].components_ :
			data_row += f"{row[feature]:8.4f}"
		print(data_row)

	verbose_int = 2
	print(f"Starting Bayes BayesianGaussianMixture with {args.iter:n} iterations.")
else:
	verbose_int = 0

dpgmm = mixture.BayesianGaussianMixture(n_components=16,
                                        covariance_type='full',
                                        max_iter=2000,
                                        n_init=args.iter,
                                        verbose=verbose_int,
                                        verbose_interval = 100,
                                        weight_concentration_prior_type='dirichlet_process',
                                        weight_concentration_prior=1e-8).fit(X)




df['group'] = dpgmm.predict(X)+1

df['PC1'] = X[:,0]
df['PC2'] = X[:,1]

group_file = open(args.output, 'w')
df.to_csv(group_file, index = False, float_format='%.3f')
group_file.close()
# print(groups.max())
groups = dpgmm.predict(X)+1
for i in range(1,groups.max()+1): print('{} has occurred {} times'.format(i, groups.tolist().count(i))) 

