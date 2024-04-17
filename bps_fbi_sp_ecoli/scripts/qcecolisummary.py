#!/usr/bin/env python3

import glob
import sys
import argparse
import os.path
import re
import os
from os import walk
import subprocess
from random import randint


parser = argparse.ArgumentParser(description='Summarizes the ecoli reports after running bifrostecolipostkma.py.')
parser.add_argument("-i", "--input_dir", default='input')
parser.add_argument("-o", "--output_dir", default='output')
args = parser.parse_args()


# Input variables
input_dir=args.input_dir
output_dir=args.output_dir


# Functions
def make_folder_if_not_exists(folder_name):
	"""This function creates output folders if they don't exists."""
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)

def find_dirs(path):
	"""This function finds dirs in a path,
	excluding files."""
	dirs = [os.path.join(path, dir) for dir in os.listdir(path) if os.path.isdir(os.path.join(path, dir))]
	return dirs

def find_files(path):
	"""This function finds files in a path,
	excluding dirs."""
	files = [os.path.join(path, file) for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))]
	return files


# Main
# Create outputs
make_folder_if_not_exists(output_dir)


sampledirs = find_dirs(input_dir)
for sampledir in sampledirs:
	outfile_path = os.path.join(output_dir, os.path.basename(input_dir))
	print(outfile_path)
	outfile=open(f"{outfile_path}.tsv",'w')
	header="false"
	files = find_files(sampledir)
	for file in files:
		if header=="false":
			proc=subprocess.Popen("".join(["head -1 ", file]), stdout=subprocess.PIPE, shell=True)
			(out, err)=proc.communicate()
			outfile.write(out.decode())
			header="true"
		proc=subprocess.Popen("".join(["tail -n +2 ", file]), stdout=subprocess.PIPE, shell=True)
		(out, err)=proc.communicate()
		outfile.write(out.decode())
