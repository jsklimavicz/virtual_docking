#analyze_pockets.py
import sys
import argparse
import os, errno
import glob
import re

fpocket_folder = sys.argv[1]
if(fpocket_folder[-1]=="/"):
	fpocket_folder = fpocket_folder[1:(len(fpocket_folder)-1)]
protein_name = os.path.split(fpocket_folder)[-1].split("_")[0]


print(fpocket_folder)
print(protein_name)

