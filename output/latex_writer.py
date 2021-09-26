#!/usr/bin/env python3
#late_writer.py
import sys
import argparse
import os, errno
import glob
import re
import time
import numpy as np
import math
import pandas as pd
import time



class LatexWriter:
	gridx = 2
	gridy = 3
	bs = '\\'

	def __init__(self, img_folder, title = "Merlin Bioassay Results"):
		self.image_folder = img_folder
		self.make_header()
		self.make_title(title)
		self.cmpd_graph_list = []
		self.graph_count = 0
		self.make_end()

	def make_header(self):
		self.header = ["\\documentclass{article}"]
		self.header.append("\\usepackage[letterpaper,left=0.75in, right=0.75in, top=0.75in, bottom=0.75in]{geometry}")
		self.header.append("\\usepackage{graphicx}")
		self.header.append("\\usepackage[font=scriptsize]{subcaption}")
		self.header.append("\\captionsetup[subfigure]{format=hang,justification=raggedright,singlelinecheck=false}")
		self.header.append("\\usepackage{amsmath,amsthm,amsfonts,amssymb,mathtools}")
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

	def make_cmpd_caption(self, name, lc50, lc50CI, R2, reps, lcTrue = True):
		'''
		<Compound Name>
		LC50: <LC50> <CI>; Curve r^2: <R2> 
		Biological Replicates: <Bio reps>
		'''
		caption = f'\\textbf{{{name}}}' 
		if lcTrue: caption += ' LC$_{50}$: ' + f'{lc50} ppm {lc50CI} \\\\ \n'
		else: caption += ' LC$_{50}$ cannot be estimated with the given data. \\\\ \n'

		if reps == 1: caption += f'{reps} biological replicate; R$^2$: {R2}'
		else: caption += f'{reps} biological replicates; R$^2$: {R2}'
		#TODO: add reference compound rel. potency
		return caption

	def make_cmpd_graph(self, image_dir, name, lc50, lc50CI, R2, reps, lcTrue = True):
		graph_lines = []
		image_width = 1/self.gridx
		img_line = f"\\includegraphics[width = {{0.95\\textwidth}}]{{{image_dir:s}}}"
		name, alpha_name = self.latexify_name(name)
		print(alpha_name, name)
		caption = f'\\caption*{{{self.make_cmpd_caption(name, lc50, lc50CI, R2, reps, lcTrue = lcTrue)}}}'
		graph_lines.append(f"   \\begin{{subfigure}}{{{image_width:.3f}\\textwidth}}")
		graph_lines.append("      \\centering")
		graph_lines.append("      "+img_line)
		graph_lines.append("      \\vspace{-0.05cm}")
		graph_lines.append("      "+caption)
		graph_lines.append("      \\vspace{0.1cm}")
		graph_lines.append("   \\end{subfigure}%")
		self.graph_count += 1
		self.cmpd_graph_list.append((alpha_name, graph_lines))


	def make_body(self):
		self.body = []
		self.graph_count = 0
		row_end = False
		
		# print(self.cmpd_graph_list)
		self.cmpd_graph_list.sort(key=lambda x:x[0])
		# print(self.cmpd_graph_list)
		for index, (alpha_name, graph_lines) in enumerate(self.cmpd_graph_list):
			if self.graph_count % (self.gridx * self.gridy) == 0:  self.body.append("\\begin{figure}[thp!]")
			for i in graph_lines: self.body.append(i)
			# self.body.append("\n")
			self.graph_count += 1
			if self.graph_count % (self.gridx * self.gridy) == 0: 
				self.body.append("\\end{figure}")
				self.graph_count = 0
				self.body.append("\\clearpage")
				self.body.append("\\pagebreak")
			if self.graph_count % self.gridx == 0: self.body.append("\\vspace{-0.1cm}")

		if self.graph_count % (self.gridx * self.gridy) != 0: self.body.append("\\end{figure}")
		self.body.append("\\pagebreak")
		self.body.append("\n\n")



	def write_file(self, out_path):
		self.make_body()
		file = open(out_path, 'w')
		file.write("\n".join(self.header))
		file.write("\n".join(self.title))
		file.write("\n".join(self.body))
		file.write("\n".join(self.end))
		file.close()
		self.make(out_path)

	def make(self, out_path):
		# os.system(f'pdflatex --shell-escape {out_path:s} > pdf_tex.log 2>&1')
		# os.system(f'pdflatex {out_path:s} > pdf_tex.log 2>&1')
		os.system(f'pdflatex -interaction=nonstopmode {out_path:s} --shell-escape  --enable-write18  > pdf_tex.log 2>&1')


	def latexify_name(self, name):
		greek_letters = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'lambda']
		mod = False
		alpha_name = name
		if name[:2].lower() == "s-": 
			alpha_name = name[2:]
			name = "\\textit{S}-" + alpha_name
			mod = True
		elif name[:2].lower() == "r-": 
			alpha_name = name[2:]
			name = "\\textit{R}-" + alpha_name
			mod = True
		else:
			for letter in greek_letters:
				if name[:len(letter)+1].lower() == letter+'-': 
					alpha_name = name[len(letter)+1:]
					name = "$\\" + letter + '$' + alpha_name
					mod = True
					break
		if len(name)< 5:
			return name.upper(), alpha_name.lower()
		elif mod:
			return name, alpha_name.lower()
		else:
			return name[0].capitalize() + name[1:], alpha_name.lower()
