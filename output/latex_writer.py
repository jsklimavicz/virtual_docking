#!/usr/bin/env python3
#late_writer.py
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
import time


parser = argparse.ArgumentParser(description="""Produces pdf file for virtual screening results""")
parser.add_argument("grouped_data", 
	help="Clustered data file from ""bayesMM.py"", in .csv format")
parser.add_argument("-o", metavar="output_filename", dest = "output",
	help="Name (with extension) of the output pdf file without extension. Default is compound_structures.",
	default = "compound_structures")
parser.add_argument("-p", metavar="pages", dest = "p", type = int,
	help="Number of lines to include in the output file. Default is to print all.",
	default = 5)
parser.add_argument("-i", metavar="imagepath", dest = "imagepath",
	help="Directory of where .png images are/should be stored. Default is ""./images"" ",
	default = "./images")
parser.add_argument("-v", "--verbose", 
	help="Prints out log info to stdout.",
	default=False, action="store_true")


class Latex_Writer:
	gridx = 4
	gridy = 6
	bs = '\\'

	def __init__(self, filename, data, args, title = "Chemical Structures of Ligand Binding Compounds"):
		self.verbose = args.verbose
		self.image_folder = args.imagepath
		self.output = args.output + ".tex"
		self.df = data
		self.page_max = args.p
		self.make_header()
		self.make_title(title)
		self.make_header()
		self.section_count = 0
		self.images_made = 0
		self.section_bodies = []
		self.section_driver()
		self.make_end()
		self.write_file()


	def make_header(self):
		self.header = ["\\documentclass{article}"]
		self.header.append("\\usepackage[utf8]{inputenc}")
		self.header.append("\\usepackage[letterpaper,left=0.75in, right=0.75in, top=0.75in, bottom=0.75in]{geometry}")
		#Packages for imagesPackages for images
		self.header.append("\\usepackage{graphicx}")
		self.header.append("\\usepackage[dvipsnames]{xcolor}")
		self.header.append("\\definecolor{TTgreen}{cmyk}{0.56, 0, 82, 0.26}")
		# self.header.append("\\usepackage{svg}")
		self.header.append("\\usepackage[font=scriptsize]{subcaption}")
		self.header.append("\\captionsetup[subfigure]{format=hang,justification=raggedright,singlelinecheck=false}")
		#Packages for symbols and chemistry text
		self.header.append("\\usepackage{amsmath,amsthm,amsfonts,amssymb,mathtools}")
		self.header.append("\\usepackage[version=4]{mhchem}")
		self.header.append(f"\\graphicspath{{{{{self.image_folder:s}}}}}")
		self.header.append("\\setcounter{topnumber}{8}")
		self.header.append("\\setcounter{bottomnumber}{8}")
		self.header.append("\\setcounter{totalnumber}{8}")
		self.header.append("\\usepackage{fancyhdr}")
		self.header.append("\\pagestyle{fancy}")
		self.header.append("\\fancyhf{}")
		self.header.append("\\fancyhead[L]{\\rightmark}")
		self.header.append("\\fancyfoot[C]{\\thepage}")
		self.header.append("\\fancyhead[R]{\\today}")
		self.header.append("\\renewcommand{\\headrulewidth}{0pt}")
		self.header.append("\n\n")




		self.header.append("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")


	def make_title(self, title):
		self.title = [f"\\title{{{title:s}}}"]
		self.title.append("\\author{James Klimavicz}")
		self.title.append("\\date{\\today}")
		self.title.append("\\begin{document}")
		self.title.append("\n\n")


	def make_end(self):
		self.end = ["\n\n"]
		self.end.append("\\end{document}")


	def zinc_png_checker(self, zincID, smiles):
		mod_id = "_".join(zincID.split("."))
		image_name = self.image_folder + "/" + mod_id 
		if not os.path.exists(image_name+".eps"): 
			self.images_made += 1
			if self.verbose and (self.images_made % 10 == 0) : print(f"Making image {self.images_made:n}")
			os.system(f'obabel -:"{smiles:s}" -O {image_name:s}.svg -xp 144  -xC -xs > /dev/null 2>&1')
			with open(image_name+".svg", 'r') as file_in:
				with open("temp.svg", 'w') as file_out:
					for line in file_in:
						line = line.replace('width="144px" height="144px" viewBox="0 0 100 100">', 'width="144px" height="96px" viewBox="0 0 144 96">')
						line = line.replace('width="100" height="100"', 'width="144px" height="96px"')
						file_out.write(line)
			os.system(f'cp temp.svg {image_name:s}.svg')
			os.system(f'inkscape {image_name:s}.svg -E {image_name:s}.eps > /dev/null 2>&1')
		return image_name + ".eps"

	def mod_zinc_id(self, zincID):
		id = zincID.split(".")[0]
		res = re.split('[ZINC|JSK]+[0]+', id)
		if "JSK" in id: 
			return "JSK" + res[1]
		else: 
			return "ZINC" + res[1]
		

	def make_section(self, section_dict):
		section = [f"\\section*{{Group {self.section_count + 1:n}}}"]
		section.append("\\vspace{-0.7cm}")
		img_count = 0
		pages = 0
		row_end = False
		image_width = 1/self.gridx
		for index, (row, row_dict) in enumerate(section_dict.items()):
			smiles = row_dict["SMILES"]
			zincID = self.mod_zinc_id(row_dict["ZincID"])
			image_name = self.zinc_png_checker(zincID, smiles)
			zincID = self.mod_zinc_id(row_dict["ZincID"])
			if "JSK" in zincID:
				caption = f'\\caption*{{\\textcolor{{TTgreen}}{{TT-{zincID:s}}} '
			else: 
				caption = f'\\caption*{{{zincID:s} '
			caption += f'E: {float(row_dict["Energy"]):.1f}, F: {float(row_dict["Efficiency"]):.3f}\\\\ \\ce{{{row_dict["Formula"]:s}}}, MW: {float(row_dict["MW"]):.2f} \\\\ {row_dict["Time"]:s}}}'
			line = f"\\includegraphics[width = {{0.95\\textwidth}}]{{{image_name:s}}}"
			if img_count % (self.gridx * self.gridy) == 0:  section.append("\\begin{figure}[thp!]")
			section.append(f"   \\begin{{subfigure}}{{{image_width:.3f}\\textwidth}}")
			section.append("      \\centering")
			section.append("      "+line)
			section.append("      \\vspace{-0.35cm}")
			section.append("      "+caption)
			section.append("   \\end{subfigure}%")
			img_count += 1
			if img_count % (self.gridx * self.gridy) == 0: 
				section.append("\\end{figure}")
				img_count = 0
				pages += 1
				if pages >= self.page_max:
					break
				section.append("\\clearpage")
				section.append("\\pagebreak")
			if img_count % self.gridx == 0: section.append("\\vspace{-0.1cm}")

		if img_count % (self.gridx * self.gridy) != 0:
			section.append("\\end{figure}")
		section.append("\\pagebreak")
		section.append("\n\n")
		self.section_bodies.append(section)
		self.section_count += 1	


	def write_file(self):
		file = open(self.output, 'w')
		file.write("\n".join(self.header))
		file.write("\n".join(self.title))
		for section in self.section_bodies:
			file.write("\n".join(section))
		file.write("\n".join(self.end))
		file.close()


	def section_driver(self):
		unique_groups = self.df.group.unique()
		for group in unique_groups:
			# "ZincID,Time,ReceptorID,nHeavyAtoms,Model,Energy,Efficiency,Formula,MW,nRings,logP,PSA,MR,SMILES\n"
			df = self.df.loc[self.df['group'] == group]
			df = df[["ZincID", "Time", "Energy", "Efficiency", "Formula", "MW", "logP", "SMILES"]]
			section_dict = df.to_dict('index')
			if self.verbose: print(f"On section {self.section_count+1:n}")
			self.make_section(section_dict)
			# if self.section_count >= 1 : break #for debug



if __name__ == "__main__":
	args = parser.parse_args()
	df=pd.read_csv(args.grouped_data, sep=',')
	# print(args)
	writer = Latex_Writer(filename = "temp.tex", data = df, args = args)
	if args.verbose: print("Writing pdf file...")
	# subprocess.call(["zsh","-c",f'pdflatex {args.output:s}.tex > /dev/null'])
	# subprocess.call(["zsh","-c",f'pdflatex {args.output:s}.tex'])
	os.system(f'pdflatex --shell-escape {args.output:s}.tex > pdf_tex.log 2>&1')
	if args.verbose: print("Done!")




