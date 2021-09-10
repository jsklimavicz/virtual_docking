#!/usr/bin/env python3
#fpocket_bounds.py
import sys
import argparse
import os, errno
import glob
import re
import math


parser = argparse.ArgumentParser(description="Finds the amino acid residues surrounding one or more pockets found by fpocket.")
parser.add_argument("fpocket_folder", 
	help="Diretory of output from fpocket")
parser.add_argument("pocket", 
	help="A string of one or more pockets, seperated by commas or dashes with no spaces")
parser.add_argument("-p", "--pymol", 
	help="If flag is included, a string of residue numbers seperated by '+' is output. This is useful for selecting residues in pymol.", 
	default=False, action="store_true")
parser.add_argument("-v", "--vina",
	help="If flag is included, residues bordering the pocket are listed in the format RES1_RES2_..._RESN, aas required for autodock's prepare_flexreceptor4.py for vina", 
	default=False, action="store_true")
parser.add_argument("-a", "--autodock",
	help="This is an alias for --vina", 
	default=False, action="store_true")
parser.add_argument("-b", "--bounds", 
	help="Print the x, y, and z centers and sizes for the pocket to stdout.", 
	default=False, action="store_true")
parser.add_argument("-B", metavar="box_padding",type=float, 
	help="Padding for the size of the boxes for an autodock config file.", 
	default=1)
parser.add_argument("-c", "--conf", 
	help="Setting this flag makes an autodock vina config file containing x, y, and z centers and sizes.", 
	default=False, action="store_true")
parser.add_argument("-r", metavar="repeats", type=int, 
	help="Number of atoms from a residue that must be included in a pocked boundary to be included in the output lists. Default is 2.", 
	default=2)
parser.add_argument("-o", metavar="output",type=float, 
	help="Name (and extensionm, if desired) of the output file.", 
	default="")

args = parser.parse_args()


#gets unique items in a list without converting to a set
def get_unique(numbers):
    unique = []
    for number in numbers:
        if number in unique:
            continue
        else:
            unique.append(number)
    return unique

#removes one instance of an element from each unique item in a list
def listed_more_than_once(numbers):
    unique = get_unique(numbers)
    for number in unique:
    	numbers.remove(number)
    numbers = get_unique(numbers)
    return numbers


fpocket_folder = args.fpocket_folder

if(fpocket_folder[-1]=="/"):
	fpocket_folder = fpocket_folder[0:(len(fpocket_folder)-1)]
protein_name = os.path.split(fpocket_folder)[-1].split("_")[0]

#pockets argument should be of form pocket1+pocket2+....+pocketn, i.e. a list of integers seperated by plus signs
pocket_list = args.pocket.split("+")
pockets = []
for x in pocket_list:
	if x.isdigit() and int(x)>0:
		pockets.append(int(x))
	else:
		print("Specified binding modes were not positive integers. Please use merge_rigid_flex.py --help for options. Exiting gracefully...")
		exit()


minx = math.inf 
miny = math.inf 
minz = math.inf 
maxx = -math.inf 
maxy = -math.inf 
maxz = -math.inf 

res_list = list()

for pocket in pockets:
	pocket_file = open(fpocket_folder+"/pockets/pocket"+str(pocket)+"_atm.pdb", 'r')
	pocket_file_lines = pocket_file.readlines()
	pocket_file.close()
	for line in pocket_file_lines:
		#we need ATOM lines in which the atom name is not something in the backbone (e.g. not N, O, C, CA) or in prolines
		if(line[0:6].strip()=="ATOM" and (line[12:16].strip() not in ["N","O","C","CA"]) and (line[17:20].strip() not in ["PRO"])):
			res_num = int(line[23:27].strip())
			res_list.append(res_num)
			currx = float(line[30:38])
			curry = float(line[38:46])
			currz = float(line[46:54])
			minx = (currx if currx < minx else minx)
			maxx = (currx if currx > maxx else maxx)
			miny = (curry if curry < miny else miny)
			maxy = (curry if curry > maxy else maxy)
			minz = (currz if currz < minz else minz)
			maxz = (currz if currz > maxz else maxz)

center_x = round((minx+maxx)/2,1)
center_y = round((miny+maxy)/2,1)
center_z = round((minz+maxz)/2,1)

size_x = round(maxx - minx,1) + args.box_padding
size_y = round(maxy - miny,1) + args.box_padding
size_z = round(maxz - minz,1) + args.box_padding

size_str = "size_x = {:.1f}\nsize_y = {:.1f}\nsize_z = {:.1f}".format(size_x,size_y,size_z)
center_str = "center_x = {:.1f}\ncenter_y = {:.1f}\ncenter_z = {:.1f}".format(center_x,center_y,center_z)


if args.bounds:
	print(size_str)
	print(center_str)


res_list.sort()
unit_res = get_unique(res_list)
n_times_included = 1

mult_list = unit_res
while args.rep > n_times_included:
	mult_list = listed_more_than_once(res_list)
	n_times_included = n_times_included+1
if args.pymol:
	#print(mult_list)
	print('+'.join(str(x) for x in mult_list))

if args.autodock or args.vina:
	vina_flex = []
	pdb_file_name = fpocket_folder+"/"+protein_name+"_out.pdb"
	pdb_file = open(pdb_file_name, 'r')
	pdb_file_lines = pdb_file.readlines()
	pdb_file.close()
	for x in mult_list:
		for line in pdb_file_lines:
			if(int(line[23:27].strip()) == x):
				vina_flex.append(line[17:20].strip()+str(x))
				break

	print('_'.join(vina_flex))


if args.conf:
	if args.output=="":
		vina_conf_file = open(protein_name + "_conf.txt", 'w')
	else:
		vina_conf_file = open(args.output, 'w')
	vina_conf_file.write("receptor = " + protein_name.split(".")[0] + "_rigid.pdbqt\n")
	vina_conf_file.write("flex = " + protein_name.split(".")[0] + "_flex.pdbqt\n\n")
	vina_conf_file.write(size_str+"\n\n")
	vina_conf_file.write(center_str+"\n\n")
	vina_conf_file.write("exhaustiveness = 8\n")
	vina_conf_file.close()


