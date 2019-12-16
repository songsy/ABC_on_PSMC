#!/usr/bin/env python
# python TMRCA_decode.py
# Shiya Song
# 7 April 2015

import argparse
from NGS_utils import *
import math
import glob
import os


def TMRCA_dist_auto(file):
	f = open(file,'r')
#	f_out = open(file.replace('.bed','_TMRCA.txt'),'w')
	f_out = open(file.replace('.bed','.TMRCA'),'w')
	TMRCA=[]
	TMRCA_piece=[]
	bin_l= {}
	bin_r= {}
	for i in range(0,28):
		TMRCA.append(0)
		TMRCA_piece.append(0)
		bin_l[i]=float(0)
		bin_r[i]=float(0)
	first = True
	for line in f:
		line = line.strip().split('\t')
		if first is True:
			first = False
			old_chr = line[0]
			old_chunk = line[3]
			TMRCA_piece[int(line[3])]+=1
#		if old_chr!=line[0] or old_chunk!=line[3]:
		else:
			TMRCA_piece[int(line[3])]+=1
			old_chr = line[0]
			old_chunk = line[3]
		TMRCA[int(line[3])]+=int(line[2])-int(line[1])
		bin_l[int(line[3])]=float(line[4])
		bin_r[int(line[3])]=float(line[5])
	print sum(TMRCA),sum(TMRCA_piece)
	for i in range(0,28):
		if TMRCA_piece[i]!=0:
			print >>f_out,'%f\t%i\t%f\t%i\t%f\t%f' %(bin_l[i]/2/args.mu*args.G,TMRCA[i],float(TMRCA[i])/sum(TMRCA),TMRCA_piece[i],float(TMRCA[i])/TMRCA_piece[i],float(TMRCA_piece[i])/sum(TMRCA_piece))
		else:
			print >>f_out,'%f\t%i\t%f\t%i\t%f\t%f' %(bin_l[i]/2/args.mu*args.G,TMRCA[i],float(TMRCA[i])/sum(TMRCA),TMRCA_piece[i],0,float(TMRCA_piece[i])/sum(TMRCA_piece))

#		print >>f_out,'%f\t%i' %(bin_r[i]/2/args.mu*args.G,TMRCA[i])

def get_TMRCA_bin(file):
	TMRCA_bin=[]
	genome_percentage=[]
	rho = 0
	theta = 0
	if True:
		f=open(file,'r')
		start = False
		for line in f:
			line = line.strip()
			if line[:2]=='RD':
				number=int(line.split('\t')[1])
				if number==20:
					start = True
				else:
					start = False
			elif start is True:
				if line[:2]=='RS':
					line = line.split('\t')
					N0=theta/100/4/args.mu
#					print line[2],N0,float(line[2])*2*N0*G
					TMRCA_bin.append(float(line[2])*2*N0*args.G)
					genome_percentage.append(float(line[5]))
				if line[:2]=='TR':
					line = line.split('\t')
					rho=float(line[2])
					theta=float(line[1])
		print TMRCA_bin,len(TMRCA_bin)
		print genome_percentage,len(genome_percentage)
		print theta,rho
	return TMRCA_bin,genome_percentage,rho

def calc_chunk_percentage(genome_percentage,TMRCA_bin,rho):
	chunk_percentage=[]
	if True:
		for i in range(len(genome_percentage)-1):
			if i==0:
				a=0
			else:
				a=genome_percentage[i]*rho/(math.log(TMRCA_bin[i+1])-math.log(TMRCA_bin[i]))*(TMRCA_bin[i+1]-TMRCA_bin[i])
			chunk_percentage.append(a)
		print chunk_percentage
	return chunk_percentage

def write_output_v2(file):
	if True:
		fout=open(file.replace('.psmc','_TMRCA_calc.txt'),'w')
		for i in range(len(TMRCA_bin)-1):
			print >>fout,'%f\t%f\t%f\t%f' %(TMRCA_bin[i],genome_percentage[i],chunk_percentage[i],float(chunk_percentage[i])/sum(chunk_percentage))

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='decode TMRCA')
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--G", default=30,dest='G',help="generation time")
	
	args = parser.parse_args()

	for file in glob.glob('*decode*.psmc*'):
		cmds='decode2bed.pl %s > %s' %(file,file.replace('.psmc','.bed'))
		print cmds
		os.popen(cmds)

	for file in glob.glob('*decode*.bed*'):
		print file
		TMRCA_dist_auto(file)
	'''
	for file in glob.glob('*decode*.psmc*'):
		print file
		print 'get TMRCA bin'
		TMRCA_bin,genome_percentage,rho=get_TMRCA_bin(file)
		print 'get chunk percentage'
		chunk_percentage=calc_chunk_percentage(genome_percentage,TMRCA_bin,rho)
		print 'write TMRCA calculation output'
		write_output_v2(file)
	'''