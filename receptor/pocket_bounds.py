#!/usr/bin/env python3
#pocket_bounds.py
import sys
import argparse
import os, errno
import glob
import re
import math

#making flexible ligands requires residues to be made flexible are of the form ARG8_ILE84_..._THR443 
#i.e. residue and numbers, separated by an underscore.

parser = argparse.ArgumentParser(description="Usage: pocket_bounds.py pdb_file aa_list [--autodock] [--bounds] [-B box_padding] [--conf] [-o output]")
parser.add_argument("pdb_file", help="Protein pdb file")
parser.add_argument("aa_list", help="A string of amino acid numbers, seperated by + or , only")
parser.add_argument("-v", "--vina",
	help="If flag is included, residues bordering the pocket are listed in the format RES1_RES2_..._RESN, aas required for autodock's prepare_flexreceptor4.py for vina", 
	default=False, action="store_true")
parser.add_argument("-a", "--autodock",
	help="This is an alias for --vina", 
	default=False, action="store_true")
parser.add_argument("-b", "--bounds", 
	help="Print the x, y, and z centers and sizes for the pocket to stdout.", 
	default=False, action="store_true")
parser.add_argument("-B", "--box_padding",type=float, 
	help="Padding for the size of the boxes for an autodock config file.", 
	default=1)
parser.add_argument("-c", "--conf", 
	help="Setting this flag makes an autodock vina config file containing x, y, and z centers and sizes.", 
	default=False, action="store_true")
parser.add_argument("-o", "--output",
	help="Name (and extensionm, if desired) of the output file.", 
	default="")


args = parser.parse_args()

#gets unique items in a list without converting to a set
def get_unique(aa_list):
    unique = []
    for aa in aa_list:
        if aa in unique:
            continue
        else:
            unique.append(aa)
    return unique

protein_name = os.path.basename(args.pdb_file)
aa_list = re.split("\+|,", args.aa_list)

minx = math.inf 
miny = math.inf 
minz = math.inf 
maxx = -math.inf 
maxy = -math.inf 
maxz = -math.inf 


vina_flex = []
pdb_file = open(args.pdb_file, 'r')
pdb_file_lines = pdb_file.readlines()
pdb_file.close()
for aa in aa_list:
	for line in pdb_file_lines:
		if(int(line[23:27].strip()) == int(aa)):
			vina_flex.append(line[17:20].strip()+str(aa))
			currx = float(line[30:38])
			curry = float(line[38:46])
			currz = float(line[46:54])
			minx = (currx if currx < minx else minx)
			maxx = (currx if currx > maxx else maxx)
			miny = (curry if curry < miny else miny)
			maxy = (curry if curry > maxy else maxy)
			minz = (currz if currz < minz else minz)
			maxz = (currz if currz > maxz else maxz)
		elif(int(line[23:27].strip()) > int(aa)):
			break


vina_flex = get_unique(vina_flex)
if args.vina or args.autodock:
	print('_'.join(vina_flex))

center_x = round((minx+maxx)/2,1)
center_y = round((miny+maxy)/2,1)
center_z = round((minz+maxz)/2,1)

size_x = round(maxx - minx + args.box_padding ,1)
size_y = round(maxy - miny + args.box_padding ,1)
size_z = round(maxz - minz + args.box_padding ,1)

size_str = "size_x = {:.1f}\nsize_y = {:.1f}\nsize_z = {:.1f}".format(size_x,size_y,size_z)
center_str = "center_x = {:.1f}\ncenter_y = {:.1f}\ncenter_z = {:.1f}".format(center_x,center_y,center_z)


if args.conf or args.output != "":
	if args.output=="":
		vina_conf_file = open(protein_name.split(".")[0] + "_conf.txt", 'w')
	else:
		vina_conf_file = open(args.output, 'w')
	vina_conf_file.write("receptor = " + protein_name.split(".")[0] + "_rigid.pdbqt\n")
	vina_conf_file.write("flex = " + protein_name.split(".")[0] + "_flex.pdbqt\n\n")

	vina_conf_file.write(size_str+"\n\n")
	vina_conf_file.write(center_str+"\n\n")
	vina_conf_file.write("exhaustiveness = 8\n")
	vina_conf_file.close()

if args.bounds:
	print(center_str)
	print(size_str)
