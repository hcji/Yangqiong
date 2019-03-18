# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 16:33:18 2019

@author: hcji
"""

from scipy.stats.stats import pearsonr
from subprocess import call
import numpy as np
import math
import sys 
import glob 
import subprocess 
import os 
import shutil 
import collections
import time
import platform
import csv

spartan_path = 'C:/Program Files/Wavefunction/Spartan14v114'

header = [
	"C OPT B3LYP 6-31G* FINDBEST=MMFF FREQ NMR POLAR\n",
	"C CONVERGE SCFTOLERANCE=VERYHIGH GRADIENTTOLERANCE=0.000005\n",
	"DISTANCETOLERANCE=0.00002 BIGGRID\n"
]
footer = [
	"NOCONFORMER\n",
	"BEGINPROPIN\n",
	"C DOQSAR PRINTMO THERMO PRINTFREQ PRINTCHG\n", 
	"C QSAR NBO DIPOLE POP PROPPRINTLEV=2 PROPRERUN\n",
	" PRINTVIBCOORDS PROP:IGNORE_WARN\n",
	"ENDPROPIN\n",
	"BEGINPREFERENCES\n", 
	" MM:CONF_SELECTION_RULE=2\n",
	"ENDPREFERENCES\n"
]
label = "*"
atomic_weights = {
	'H': 1.0079, 
	'C': 12.0107, 
	'N': 14.0067, 
	'O': 15.9994, 
	'F': 18.9984,
	'P': 30.9737,
	'S': 32.065, 
	'Cl': 35.453, 
	'Br': 79.904,
	'I': 126.9045,
}

def in_between(text_list, start_keyword, end_keyword): 
	"""Takes a list of text and returns the list containing entries between
	start_keyword and end_keyword (inclusive).
	"""
	output = []
	parsing = False
	for line in text_list:
		if start_keyword in line: 
			parsing = True 
		if parsing:
			output.append(line) 
		if end_keyword in line: 
			parsing = False
	return output

def spartan_calculation(fpath):
	head, tail = os.path.split(fpath)  
	# if output file exists, calculation is not run
	if os.path.exists(fpath + ".spardir\\M0001\\output"): 
		print('Spartan calculation not run; output file already exists for ' + tail)
	else:
		# make file.spardir directory 
		if os.path.exists(fpath + ".spardir"):
			shutil.rmtree(fpath + ".spardir") 
		os.mkdir(fpath + ".spardir")
		# make empty _spartandir file 
		open(fpath + '.spardir\\_spartandir', 'w').close()
		# make empty M0001 directory
		os.mkdir(fpath + '.spardir\\M0001')
		# make empty _spartan file 
		open(fpath + '.spardir\\M0001\\_spartan', 'w').close()
		# extract molecular information from .spinput file created by user
		with open(fpath + '.spinput', 'r') as f:
			file_text = f.readlines() 
		molecular_data = in_between(file_text, 'M0001', 'ENDHESS')
		# make new input file in M0001 directory with calculation parameters
		with open(fpath + '.spardir\\M0001\\input', 'w') as f:
			f.writelines(header)
			f.writelines(molecular_data) 
			f.writelines(footer) 
		# run calculation via command line 
		if platform.system()=='Windows':
			# spartan_path = 'C:\\Program Files\\Wavefunction\\Spartan14v114'
			spartanx_path = spartan_path + '\\spartanx' 
			command = '"' + spartanx_path + '" --path "' + spartan_path + \
				'" --foreground-submit "' + fpath + '.spardir"'
			print('Calculation running for ' + tail + '...')
			print(command)
			call(command)
			print('Calculation complete')
		else:
			raise RuntimeError('Only Windows is currently supported')

def run_spartan(dpath):
    # dpath = 'Spartan/data'
    fs = os.listdir(dpath)
    ms = list(set([f.split('.')[0] for f in fs]))
    try:
        for m in ms:
            fpath = dpath + '/' + m
            spartan_calculation(fpath)
        return 1
    except:
        return 0
    