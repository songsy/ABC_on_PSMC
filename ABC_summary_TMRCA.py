#!/usr/bin/env python
# python calc_TMRCA_from_macs.py
# Shiya Song
# 17th April 2015

import argparse
from NGS_utils import *
import numpy as np
import pickle
import math
import glob
import os
import sys
COPY = 4
G=30


def summary_calc_TMRCA(file):
	chunk_percentage_bin=[[0 for x in range(63)] for x in range(COPY)] 
	genome_percentage_bin=[[0 for x in range(63)] for x in range(COPY)]
	for i in range(args.nchr):
		f=open(file+'_TMRCA_%s.txt' %(i),'r')
		for line_number,line in enumerate(f):
			line=line.strip().split('\t')
			for j in range(4):
#				print line_number,j,line[j*2+1],line[j*2+2]
				genome_percentage_bin[j][line_number]+=float(line[j*2+1])
				chunk_percentage_bin[j][line_number]+=float(line[j*2+2])
		print i,'done'
	f=open(file+'_TMRCA_0.txt','r')
	fout=open(args.file+'_TMRCA_test4.txt','w')
	for line_number,line in enumerate(f):
		line=line.strip().split('\t')
		output = '%s' %(line[0])
		for index in range(COPY):
			output+='\t%f\t%f' %(float(genome_percentage_bin[index][line_number])/sum(genome_percentage_bin[index]),float(chunk_percentage_bin[index][line_number])/sum(chunk_percentage_bin[index]))
		print >>fout,output
	for i in range(args.nchr):
		os.popen('rm '+file+'_TMRCA_%s.txt' %(i))
	return genome_percentage_bin,chunk_percentage_bin

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--nchr", default=100,type=int,dest='nchr',help="number of chromosome")
	args = parser.parse_args()

	summary_calc_TMRCA(args.file)
