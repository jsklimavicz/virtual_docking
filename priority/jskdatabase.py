#!/usr/bin/env python3
#jskdatabase.py
import sys
import argparse
import os, errno
import glob
import re
import time

parser = argparse.ArgumentParser(description="""Adds new JSK-made compounds to the JSK database fo later reference.""")
parser.add_argument("old_database", type = str,
	help="JSK database filename. Should include the .csv extension.")
parser.add_argument("smi_file", type = str,
	help="New smi file for assignment/addition of SMILES to JSK database.")


class JSKDataset:
	def __init__(self, old_database, new_smi_file):
		self.old = old_database
		self.new = new_smi_file
		self.maxID = 0
		self.curr_list = []
		self.pendingID_list = []
		self.smi_list = []
		self.read_old_file()
		self.read_new_file()
		self.process_pending_line()

	def read_old_file(self):
		smi_file = open(self.old, 'r')
		jsk_lines = smi_file.readlines()
		smi_file.close()
		for line in jsk_lines: self.process_old_line(line)
			
	def process_old_line(self, line):
		line_split = line.split(",")
		line_split = [x.strip() for x in line_split]
		if line_split == "" or line_split == []: return
		self.process_complete_line(line_split)

	def id_num(self, cmpd_id):
		return int(re.findall('([0-9]+)',cmpd_id)[0])

	def read_new_file(self):
		smi_file = open(self.new, 'r')
		smi_lines = smi_file.readlines()
		smi_file.close()
		for line in smi_lines: 
			if line == "" or line == []: continue
			self.process_new_line(line)
		#now process pending lines

	def process_new_line(self, line):
		line_split = line.split()
		# print(line)
		# print(line_split)
		if line == "" or line_split == []: return
		smiles = self.reduce_smiles(line_split[0])
		if not self.new_smiles(smiles): return
		if len(line_split) == 1:
			self.pendingID_list.append(smiles)
		else: 
			self.process_complete_line(line_split, new = True)

	def process_complete_line(self, line_split, new = False):
		cmpd_id = line_split[1]
		smiles = self.reduce_smiles(line_split[0])
		cmpd_num = self.id_num(cmpd_id)
		if "ZINC" in cmpd_id: 
			self.smi_list.append([smiles,cmpd_id,cmpd_num])
			return
		if cmpd_num > self.maxID: self.maxID = cmpd_num
		self.curr_list.append([smiles,cmpd_id,cmpd_num])
		if new: self.smi_list.append([smiles,cmpd_id,cmpd_num])

	def new_smiles(self, smiles):
		for line in self.curr_list:
			if smiles == line[0]: return False
		return True

	def process_pending_line(self):
		for line in self.pendingID_list:
			if not self.new_smiles(line): break
			smiles = line
			self.maxID += 1
			cmpd_num = self.maxID
			cmpd_id = "JSK" + str(cmpd_num).zfill(7)
			self.curr_list.append([smiles,cmpd_id,cmpd_num])
			self.smi_list.append([smiles,cmpd_id,cmpd_num])

	def print_list(self):
		self.curr_list.sort(key = lambda x:int(re.findall('([0-9]+)',x[1])[0]))
		smi_file = open(self.old, 'w')
		for line in self.curr_list:
			smi_file.write(f"{line[0]:s},{line[1]:s}\n")
		smi_file.close()

	def update_smi_file(self):
		smi_file = open(self.new, 'w')
		for line in self.smi_list:
			smi_file.write(f"{line[0]:s} {line[1]:s}\n")
		smi_file.close()

	def reduce_smiles(self,smiles):
		# breaks = list(set(re.findall('([%][0-9]+)', smiles)))
		# break_dict = {}
		# ring_count = 0
		# for item in breaks: 
		# 	ring_count += 1
		# 	break_dict[item] = ring_count
		# for key in break_dict: smiles = smiles.replace(key, str(break_dict[key]))
		return smiles






if __name__ == "__main__":
	args = parser.parse_args()
	writer = JSKDataset("jsk_database.csv", "substances.smi")
	writer.print_list()
	writer.update_smi_file()



