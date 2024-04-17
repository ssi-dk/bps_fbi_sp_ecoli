#!/usr/bin/env python3
# bifrostecolipostkma.py


import glob
import sys
import argparse
import os.path
import re
import os
from os import walk
import subprocess


parser = argparse.ArgumentParser(description='Determines the serotype and virulence in ecoli through kmer alignment')
parser.add_argument("-i","--sampleid", help='ecoli1')
parser.add_argument("-R1", "--read1", help="ecoli1_R1.fastq.gz")
parser.add_argument("-R2", "--read2", help="ecoli2_R2.fastq.gz")
parser.add_argument("-stbit", help="ST string", default='ST:NA,NA')
parser.add_argument("-db", "--db_path", default="db")
parser.add_argument("-k", "--kma_path", default='usr/bin/kma')
parser.add_argument("-o", "--output_dir", default='output')
args=parser.parse_args()


# Input variables
sampleid=args.sampleid
read1 = args.read1
read2 = args.read2
stbit = args.stbit
db_path = args.db_path
kma_path=args.kma_path
output_dir=args.output_dir


# Functions
def get_rundir(read1) -> str:
	"""This function gets the rundir from a dir path."""
	rundir = os.path.basename(os.path.dirname(read1))
	return rundir

def make_folder_if_not_exists(folder_name):
	"""This function creates output folders if they don't exists."""
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)

def convert_number_to_date(number):
	"""This function converts a number to a date from 221117 to a 2022-11-17."""
	year=int(str(number)[:2])
	month=int(str(number)[2:4])
	day=int(str(number)[4:6])
	# Create a date string in the "YYYY-MM-DD" format
	date_string=f"20{year:02d}-{month:02d}-{day:02d}"
	return date_string

def get_wgs_date_and_number(rundir):
	"""This function separates the rundir into date and experiment_name/wgsnumber."""
	#231006_NB551234_0051_N_WGS_743_AHNLHHAFX5
	rundir_list=rundir.split("_") 	# ['231006', 'NB551234', '0051', 'N', 'WGS', '743', 'AHNLHHAFX5']
	wgsdate=convert_number_to_date(rundir_list[0]) # 231006 -> 2023-10-06
	wgsnumber="_".join(rundir_list[3:6]) # N_WGS_743
	return wgsdate, wgsnumber

def print_header_to_output(OUTFILE):
	header="isolate\twzx\twzy\tfliC\tOH\tstx\teae\tehx\tother\twgsrun\twgsdate\tST\tSTgenes\tverbose\n"
	outfile=open(OUTFILE, 'w')
	print(header, end='', file=outfile)
	return header


# Main
# Hardcoded paths
ECOLIGENESDB=os.path.join(db_path, "ecoligenes")
SPECOLIFBIDIR=os.path.join(output_dir, sampleid, "sp_ecoli_fbi") # KMA results
OUTFILE=os.path.join(output_dir, sampleid, f"{sampleid}.tsv")
					 
# Create outputs
make_folder_if_not_exists(SPECOLIFBIDIR)
header = print_header_to_output(OUTFILE)

# TODO: Remember to check if the db exists if not give a warning and index it
#kma_index -i ecoligenes.fsa -o ecoligenes

rundir = get_rundir(read1)
wgsdate, wgsnumber = get_wgs_date_and_number(rundir)

ST_list=args.stbit.split(",") # ['ST:11', 'adk:12', 'fumC:12', 'gyrB:8', 'icd:12', 'mdh:15', 'purA:2', 'recA:2']
ST=ST_list.pop(0).split(":")[1] # 11
STgenes=",".join(ST_list) # adk:12,fumC:12,gyrB:8,icd:12,mdh:15,purA:2,recA:2

# Run KMA
command=f"{kma_path} -ipe {args.read1} {args.read2} -matrix -t_db {ECOLIGENESDB} -o {SPECOLIFBIDIR}/colipost"
print(f"Command {command}") 
os.system(command)


# KMA results processing
# Settings
fliflag="false"
locations={"OH": 4, "stx":5, "wzx":1, "wzy":2, "wzt":1, "wzm":2, "fliC":3, "fli":3, "eae": 6, "ehxA":7, "wgsrun":8, "wgsdate":9, "ST": 10, "STgenes":11, "other":12}
thresholds={"stx":[98,98], "wzx":[98,98], "wzy":[98,98], "wzt":[98,98], "wzm":[98,98], "fliC":[90,90], "fli":[90,90], "eae": [95,95], "ehxA":[95,95], "other":[98,98]}
stxbadthreshold=[30,90]

# Processing of .res file
hitdict={}
hitdict[sampleid]={}
res_file=open(os.path.join(SPECOLIFBIDIR, "colipost.res"), 'r')
for line in res_file:
	if line.startswith("#"):
		continue
	#line[0]=the specific hit, line[9]=pc_idt for the hit to the ref, line[18]=the mutation of the specific line
	line=line.split("\t") 
	template=line[0]
	template_cov=line[5]
	query_ident=line[6]
	gene_list=template.split("__")
	gene_name=gene_list[1]
	serotype_or_virulence=gene_list[2]
	if gene_name.startswith("stx"):
		gene_name="stx"
	if not gene_name in hitdict[sampleid]:
		hitdict[sampleid][gene_name]=[]#{"lencov":lencov, "cov":line[5], "SNP":0, "DEL":0, "INS":0, ".":0}
	hit_results_list=[serotype_or_virulence, template_cov, query_ident]
	hitdict[sampleid][gene_name].append(hit_results_list)
print('HITDICT',hitdict)


# Processing of hits
results_dict={}
for sampleid in hitdict:
	foundgoodstx="false"
	results_dict[sampleid]=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]# Remember to increase this if "locations" is made longer  [[], [], [],[],[],[],[],[],[],[],[],,[],[][]]
	Opick="-"	
	# First the mapped hits are parsed and types are decided upon based on a logic from KNJ and SSC
	for gene_name in hitdict[sampleid]:					
		for hit_results_list in hitdict[sampleid][gene_name]:
			serotype_or_virulence=hit_results_list[0]
			template_cov=hit_results_list[1]
			query_ident=hit_results_list[2]
			if gene_name.startswith("fli") and not "fliC" in gene_name and fliflag=="false":
				results_dict[sampleid][locations["fliC"]]=[]
				fliflag="true"
			if not gene_name in thresholds:
				gene_name="other"
			if fliflag=="true":
				gene_name="fli"
			if "stx" in gene_name:
				serotype_or_virulence=serotype_or_virulence.lower()
			if "100.00" in template_cov and "100.00" in query_ident:
				if "wzt" in gene_name or "wzm" in gene_name or fliflag=="true":
					results_dict[sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				if gene_name.startswith("e"):
					serotype_or_virulence="positive"
				if gene_name.startswith("stx"):
					foundgoodstx="true"
			elif float(template_cov)> thresholds[gene_name][0] and float(query_ident)>thresholds[gene_name][1]:
				if gene_name.startswith("wz") or fliflag=="true":
					results_dict[sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				else:
					results_dict[sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				if gene_name.startswith("e"):
					serotype_or_virulence="positive"
				if gene_name.startswith("stx"):
					foundgoodstx="true"
			elif gene_name.startswith("stx") and float(template_cov)>stxbadthreshold[0] and float(query_ident)>stxbadthreshold[1]:
				results_dict[sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				serotype_or_virulence=serotype_or_virulence[:4]
			else:
				if gene_name.startswith("wz"):
					results_dict[sampleid][-1].append(":".join(["_".join([gene_name, serotype_or_virulence]), template_cov.lstrip(), query_ident.lstrip()]))
				else:
					results_dict[sampleid][-1].append(":".join([serotype_or_virulence, template_cov.lstrip(), query_ident.lstrip()]))
				continue
			if "fliC" in gene_name and fliflag=="true":
				continue
			if not serotype_or_virulence in results_dict[sampleid][locations[gene_name]]:
				results_dict[sampleid][locations[gene_name]].append(serotype_or_virulence)
	print(foundgoodstx)
	if foundgoodstx=="true":
		for stxcase in results_dict[sampleid][locations["stx"]]:
			if len(stxcase)==4:
				results_dict[sampleid][locations["stx"]].pop(results_dict[sampleid][locations["stx"]].index(stxcase))
	results_dict[sampleid][locations["stx"]]="; ".join(results_dict[sampleid][locations["stx"]]).replace("-", "")
	for gene in results_dict[sampleid]:
		if type(gene)==list:
			results_dict[sampleid][results_dict[sampleid].index(gene)]=":".join(gene)
	if len(results_dict[sampleid][locations["wzx"]]) > 1 and not "*" in results_dict[sampleid][locations["wzx"]]:
		Opick=results_dict[sampleid][locations["wzx"]]
	if len(results_dict[sampleid][locations["wzy"]]) > 1 and not "*" in results_dict[sampleid][locations["wzy"]]:
		Opick=results_dict[sampleid][locations["wzy"]]
	
	if len(results_dict[sampleid][locations["fliC"]]) < 2:
		results_dict[sampleid][locations["fliC"]]="-"
	results_dict[sampleid][locations["OH"]]=":".join([Opick, results_dict[sampleid][locations["fliC"]]])

	# Here the results_dict is curated with dashes
	for element in results_dict[sampleid]:
		if len(element)<2:
			results_dict[sampleid][results_dict[sampleid].index(element)]="-"
	results_dict[sampleid].pop(0)
	results_dict[sampleid][locations["wgsrun"]]=wgsnumber
	print(wgsnumber)
	results_dict[sampleid][locations["wgsdate"]]=wgsdate
print(results_dict)	


# Here the results_dict is curated with dashes
for sampleid in results_dict.keys():
	outfile=open(OUTFILE, 'w')
	results_dict[sampleid][locations["ST"]]=ST # as defined on line 34 KRKI 02-05-2019
	results_dict[sampleid][locations["STgenes"]]=STgenes # as defined on line 35 KRKI 02-05-2019
	lineelements=[]
	for element in results_dict[sampleid]:
		if type(element) is dict:
			for gene in ["O", "H", "stx", "eae", "ehx"]:
				if gene in element:
					lineelements.append(";".join(element[gene]))
		elif len(element)<1:
			print('element', element)
			lineelements.append("-")
		else:
			lineelements.append(element)

	
	# Print results to outfile
	print("".join([header, "\t".join([sampleid, "\t".join(lineelements)])]), file=outfile)
