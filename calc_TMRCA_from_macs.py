#!/usr/bin/env python
# python calc_TMRCA_from_macs.py
# Shiya Song
# 17th April 2015

import argparse
from NGS_utils import *
import numpy as np
from ete2 import Tree
import pickle
import math
import glob
COPY = 4
G=30

def get_TMRCA_from_macs(file,TMRCA_bin):
	TMRCA=[[] for x in range(COPY)]
	prev_index_i=[0 for x in range(COPY)]
	chunk_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)] 
	genome_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	average_length_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	print len(genome_percentage_bin),len(genome_percentage_bin[0])
	length = []
	first = [True for x in range(COPY)]
	for i in range(args.nchr):
		f=open(file+'_%s_sites.txt' %(args.start+i),'r')
		for line in f:
				line=line.strip().split('\t')
				if line[0]=='NEWICK_TREE:':
					info=line[1]
					info=info.replace(']','[')
					info=info.split('[')
					size = info[1]
					newick = info[2]
#					print size,newick
					t=Tree(newick)
					length.append(size)
					index = 0
					dist_all = []
					for j in range(int(args.nsample)):
						for hap in range(int(args.nhap)):
							a=j
							b=hap+int(args.nhap)
							dist=t.get_distance(str(a),str(b))*2*args.N0				# 2 times tree length, but here it's scaled to 4Ne, so no need to /2
							TMRCA[index].append(dist)
							dist_all.append(dist)
							index_i = find_interval(dist,TMRCA_bin[index])
							if index_i > len(TMRCA_bin[index]):
								print index_i,dist,line
							else:
								genome_percentage_bin[index][index_i]+=int(size)
								chunk_percentage_bin[index][index_i]+=1
							index +=1
#					print dist_all
		print i,'done',max(TMRCA[0]),max(TMRCA[1]),max(TMRCA[2]),max(TMRCA[3])
	for i in range(COPY):
		for j in range(len(TMRCA_bin[0])):
			if chunk_percentage_bin[i][j]!=0:
				average_length_bin[i][j]=genome_percentage_bin[i][j]/chunk_percentage_bin[i][j]
			else:
#				print i,j,TMRCA_bin[i][j],genome_percentage_bin[i][j],chunk_percentage_bin[i][j]
				average_length_bin[i][j]=0
			print genome_percentage_bin[i][j],chunk_percentage_bin[i][j],average_length_bin[i][j]
	return genome_percentage_bin,chunk_percentage_bin,average_length_bin,TMRCA,length

def get_TMRCA_from_macs_v2(file,TMRCA_bin,start):
	TMRCA=[[] for x in range(COPY)]
	prev_index_i=[0 for x in range(COPY)]
	size_chunk=[0 for x in range(COPY)]
	chunk_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)] 
	genome_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	average_length_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	print len(genome_percentage_bin),len(genome_percentage_bin[0])
	first = [True for x in range(COPY)]
	for i in range(args.nchr):
		f=open(file+'_%s_sites.txt' %(start+i),'r')
		print start+i
		for line in f:
				line=line.strip().split('\t')
				if line[0]=='NEWICK_TREE:':
					info=line[1]
					info=info.replace(']','[')
					info=info.split('[')
					size = info[1]
					newick = info[2]
#					print size,newick
					t=Tree(newick)
					index = 0
					dist_all = []
					for j in range(int(args.nsample)):
						for hap in range(int(args.nhap)):
							a=j
							b=hap+int(args.nhap)
							dist=t.get_distance(str(a),str(b))*2*args.N0				# 2 times tree length, but here it's scaled to 4Ne, so no need to /2
							dist_all.append(dist)
							index_i = find_interval(dist,TMRCA_bin[index])
							if first[index] is True:
								prev_index_i[index]=index_i
								size_chunk[index]+=int(size)
								first[index]=False
							elif index_i==prev_index_i[index]:
								size_chunk[index]+=int(size)
							else:
								genome_percentage_bin[index][prev_index_i[index]]+=size_chunk[index]
								chunk_percentage_bin[index][prev_index_i[index]]+=1
								size_chunk[index] =int(size)
								prev_index_i[index]=index_i
							index +=1
	for i in range(COPY):
		for j in range(len(TMRCA_bin[0])):
			if chunk_percentage_bin[i][j]!=0:
				average_length_bin[i][j]=genome_percentage_bin[i][j]/chunk_percentage_bin[i][j]
			else:
#				print i,j,TMRCA_bin[i][j],genome_percentage_bin[i][j],chunk_percentage_bin[i][j]
				average_length_bin[i][j]=0
#			print genome_percentage_bin[i][j],chunk_percentage_bin[i][j],average_length_bin[i][j]
	print start, 'done'
	return genome_percentage_bin,chunk_percentage_bin,average_length_bin

def get_TMRCA_from_macs_v3(file,TMRCA_bin,start,COPY):
	TMRCA=[[] for x in range(COPY)]
	prev_index_i=[0 for x in range(COPY)]
	prev_TMRCA=[0 for x in range(COPY)]
	size_chunk=[0 for x in range(COPY)]
	chunk_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)] 
	genome_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	average_length_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	print len(genome_percentage_bin),len(genome_percentage_bin[0])
	first = [True for x in range(COPY)]
	done = [True for x in range(COPY)]
	for i in range(args.nchr):
		f=open(file+'_%s_sites.txt' %(start+i),'r')
		print start+i
		for line in f:
				line=line.strip().split('\t')
				if line[0]=='NEWICK_TREE:':
					info=line[1]
					info=info.replace(']','[')
					info=info.split('[')
					size = info[1]
					newick = info[2]
#					print size,newick
					t=Tree(newick)
					index = 0
					dist_all = []
					for j in range(int(args.nsample)):
						for hap in range(int(args.nhap)):
							a=j
							b=hap+int(args.nhap)
							dist=t.get_distance(str(a),str(b))*2*args.N0				# 2 times tree length, but here it's scaled to 4Ne, so no need to /2
							dist_all.append(dist)
							index_i = find_interval(dist,TMRCA_bin[index])
							if first[index] is True:
								prev_index_i[index]=index_i
								size_chunk[index]+=int(size)
								prev_TMRCA[index]=dist
								first[index]=False
							elif index_i==prev_index_i[index] and dist==prev_TMRCA[index]:
								size_chunk[index]+=int(size)
								done[index] = False
							else:
								genome_percentage_bin[index][prev_index_i[index]]+=size_chunk[index]
								chunk_percentage_bin[index][prev_index_i[index]]+=1
								size_chunk[index] =int(size)
								prev_index_i[index]=index_i
								prev_TMRCA[index]=dist
								done[index] = True
							if COPY==1:
								break
							index +=1
						if COPY==1:
							break
	for index in range(COPY):
		if done[index]==False:
			assert size_chunk[index]!=0
			genome_percentage_bin[index][prev_index_i[index]]+=size_chunk[index]
			chunk_percentage_bin[index][prev_index_i[index]]+=1
	for i in range(COPY):
		for j in range(len(TMRCA_bin[0])):
			if chunk_percentage_bin[i][j]!=0:
				average_length_bin[i][j]=genome_percentage_bin[i][j]/chunk_percentage_bin[i][j]
			else:
#				print i,j,TMRCA_bin[i][j],genome_percentage_bin[i][j],chunk_percentage_bin[i][j]
				average_length_bin[i][j]=0
#			print genome_percentage_bin[i][j],chunk_percentage_bin[i][j],average_length_bin[i][j]
	print start, 'done'
	return genome_percentage_bin,chunk_percentage_bin,average_length_bin

def find_interval(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = (right + left)/2
		if pos < list[midpoint]:
			right = midpoint
		elif pos > list[midpoint]:
			left = midpoint+1
		elif pos == list[midpoint]:
			left = midpoint
			find = midpoint
			break
		if pos<list[midpoint] and midpoint!=0:
			midpoint=midpoint-1
	return midpoint

def get_TMRCA_bin(file):
	TMRCA_bin=[[] for x in range(COPY)]
	genome_percentage=[[] for x in range(COPY)]
	rho = [0 for x in range(COPY)]
	theta = [0 for x in range(COPY)]
	for i in range(COPY):
		f=open('%s.psmc.%i' %(file,i+2),'r')
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
					N0=theta[i]/100/4/args.mu
#					print line[2],N0,float(line[2])*2*N0*G
					TMRCA_bin[i].append(float(line[2])*2*N0)
					genome_percentage[i].append(float(line[5]))
				if line[:2]=='TR':
					line = line.split('\t')
					rho[i]=float(line[2])
					theta[i]=float(line[1])
		print TMRCA_bin[i],len(TMRCA_bin[i])
		print genome_percentage[i],len(genome_percentage[i])
		print theta[i],rho[i]
		print np.sum(genome_percentage[i])
	return TMRCA_bin,genome_percentage,rho

def get_TMRCA_bin_single(file,copy):

	TMRCA_bin=[[] for x in range(copy)]
	genome_percentage=[[] for x in range(copy)]
	rho = [0 for x in range(copy)]
	theta = [0 for x in range(copy)]
	for i in range(copy):
		f=open('%s' %(file),'r')
		start = False
		for line in f:
			line = line.strip()
			if line[:2]=='RD':
				try:
					number=int(line.split('\t')[1])
				except ValueError:
					print 'Error:',line
				if number==25:
					start = True
				else:
					start = False
			elif start is True:
				if line[:2]=='RS':
					line = line.split('\t')
					N0=theta[i]/100/4/args.mu
#					print line[2],N0,float(line[2])*2*N0*G
					TMRCA_bin[i].append(float(line[2])*2*N0)
					genome_percentage[i].append(float(line[5]))
				if line[:2]=='TR':
					line = line.split('\t')
					rho[i]=float(line[2])
					theta[i]=float(line[1])
		print TMRCA_bin[i],len(TMRCA_bin[i])
		print genome_percentage[i],len(genome_percentage[i])
		print theta[i],rho[i]
		print np.sum(genome_percentage[i])
	return TMRCA_bin,genome_percentage,rho

def calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,copy):
	chunk_percentage=[[] for x in range(copy)]
	scale_factor=[[] for x in range(copy)]
	for index in range(copy):
		chunk_percentage[index].append(0)
		scale_factor[index].append(0)
		for i in range(1,len(genome_percentage[index])-1):
			if rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i]))==0:
				a=0
				scale=0
			else:
#				a=genome_percentage[index][i]*rho[index]/(math.log(TMRCA_bin[index][i+1])-math.log(TMRCA_bin[index][i]))*(TMRCA_bin[index][i+1]-TMRCA_bin[index][i])
#				scale=1/(rho[index]/(math.log(TMRCA_bin[index][i+1])-math.log(TMRCA_bin[index][i]))*(TMRCA_bin[index][i+1]-TMRCA_bin[index][i]))
				a=genome_percentage[index][i]*rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i]))
				scale=1/(rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i])))
			chunk_percentage[index].append(a)
			scale_factor[index].append(scale)
		print chunk_percentage[index]
		print scale_factor[index]
	return chunk_percentage,scale_factor

def calc_chunk_percentage_v2(genome_percentage,TMRCA_bin,rho,copy):   # edit based on ABC_calc_TMRCA_log.py
	chunk_percentage=[[] for x in range(copy)]
	scale_factor=[[] for x in range(copy)]
	for index in range(copy):
		for i in range(len(genome_percentage[index])-1):
			if rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i]))==0:
				a=0
				scale=0
			else:
#				a=genome_percentage[index][i]*rho[index]/(math.log(TMRCA_bin[index][i+1])-math.log(TMRCA_bin[index][i]))*(TMRCA_bin[index][i+1]-TMRCA_bin[index][i])
#				scale=1/(rho[index]/(math.log(TMRCA_bin[index][i+1])-math.log(TMRCA_bin[index][i]))*(TMRCA_bin[index][i+1]-TMRCA_bin[index][i]))
				a=genome_percentage[index][i]*rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i]))
				scale=1/(rho[index]*(0.5*(TMRCA_bin[index][i+1]+TMRCA_bin[index][i])))
			chunk_percentage[index].append(a)
			scale_factor[index].append(scale)
		a=genome_percentage[index][len(genome_percentage[index])-1]*rho[index]*(TMRCA_bin[index][len(genome_percentage[index])-1])
		scale=1/(rho[index]*TMRCA_bin[index][len(genome_percentage[index])-1])
		chunk_percentage[index].append(a)
		scale_factor[index].append(scale)
		print len(chunk_percentage[index])
		print len(scale_factor[index])
	return chunk_percentage,scale_factor

def write_output():
	for index in range(COPY):
		fout=open(args.file+'_v3_TMRCA_%s.txt' %(index),'w')
		for i in range(len(TMRCA_bin[index])-1):
			print >>fout,'%f\t%f\t%f\t%f\t%f\t%f\t%f' %(TMRCA_bin[index][i],genome_percentage[index][i],float(genome_percentage_bin[index][i])/sum(genome_percentage_bin[index]),average_length_bin[index][i],scale_factor[index][i],float(chunk_percentage[index][i])/sum(chunk_percentage[index]),float(chunk_percentage_bin[index][i])/sum(chunk_percentage_bin[index]))

def write_output_v2():
	for index in range(COPY):
		fout=open(args.file+'_TMRCA_ave_length_%s.txt' %(index),'w')
		for i in range(len(TMRCA_bin[index])-1):
			print >>fout,'%f\t%f\t%f\t%f' %(TMRCA_bin[index][i],genome_percentage[index][i],scale_factor[index][i],float(chunk_percentage[index][i])/sum(chunk_percentage[index]))

def write_output_v3(index_rep,genome_percentage_bin,chunk_percentage_bin):
	TMRCA = []
	fout=open(args.file+'_%s_macs_TMRCA.txt' %(index_rep),'w')
	for i in range(len(TMRCA_bin[0])-1):
		TMRCA = [TMRCA_bin[index][i] for index in range(COPY)]
		output = '%f' %(np.mean(TMRCA))
		for index in range(COPY):
			output+='\t%f\t%f' %(float(genome_percentage_bin[index][i])/sum(genome_percentage_bin[index]),float(chunk_percentage_bin[index][i])/sum(chunk_percentage_bin[index]))
		print >>fout,output

def write_output_v4(index_rep,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,copy):
	TMRCA = []
	fout=open(args.file+'_%s_psmc_TMRCA.txt' %(index_rep),'w')
	for i in range(len(TMRCA_bin_ind[0])-1):
		TMRCA = [TMRCA_bin_ind[index][i] for index in range(copy)]
		output = '%f' %(np.mean(TMRCA))
		for index in range(copy):
			output+='\t%f\t%f' %(genome_percentage_ind[index][i],float(chunk_percentage_ind[index][i])/sum(chunk_percentage_ind[index]))
		print >>fout,output

def write_output_v5(file,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,copy):
	TMRCA = []
	fout=open(file.replace('psmc','TMRCA.txt'),'w')
	for i in range(len(TMRCA_bin_ind[0])):
		TMRCA = [TMRCA_bin_ind[index][i] for index in range(copy)]
		output = '%f' %(np.mean(TMRCA))
		for index in range(copy):
			output+='\t%f\t%f' %(genome_percentage_ind[index][i],float(chunk_percentage_ind[index][i])/sum(chunk_percentage_ind[index]))
		print >>fout,output

def write_output_v6(file,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,scale_factor_ind,copy,genome_percentage_bin,chunk_percentage_bin,average_length_bin):
	TMRCA = []
	fout=open(file+'_TMRCA.txt','w')
	for i in range(len(TMRCA_bin_ind[0])):
		TMRCA = [TMRCA_bin_ind[index][i] for index in range(copy)]
		output = '%f' %(np.mean(TMRCA))
		for index in range(copy):
			output+='\t%f\t%f\t%f\t%f\t%f\t%f' %(genome_percentage_ind[index][i],float(chunk_percentage_ind[index][i])/sum(chunk_percentage_ind[index]),scale_factor_ind[index][i],float(genome_percentage_bin[index][i])/sum(genome_percentage_bin[index]),float(chunk_percentage_bin[index][i])/sum(chunk_percentage_bin[index]),average_length_bin[index][i])
		print >>fout,output

def read_macs(macs_file,TMRCA_bin):  #diploid+4*pseudo-diploid
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			start = i
#			genome_percentage_bin,chunk_percentage_bin,average_length_bin=get_TMRCA_from_macs_v3(args.file,TMRCA_bin,start)
			psmc_file = '%s_%i.0.psmc' %(args.file,index)
#			write_output_v3(index,genome_percentage_bin,chunk_percentage_bin)
#			if os.path.isfile(psmc_file) is True:
			if True:
				TMRCA_bin_ind,genome_percentage_ind,rho_ind=get_TMRCA_bin_single(psmc_file,1)
				chunk_percentage_ind,scale_factor_ind=calc_chunk_percentage(genome_percentage_ind,TMRCA_bin_ind,rho_ind,1)
				write_output_v4(index,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,1)
			index+=1

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--psmc", dest='psmc',help="macs file prefix")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
#	parser.add_argument("--N0", default=5e04*2.36/1.25,dest='N0',help="N0")
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--start", default=0,type=int,dest='start',help="start index")
	parser.add_argument("--nrep", default=300,dest='nrep',help="number of repetition")
	parser.add_argument("--nchr", default=100,type=int,dest='nchr',help="number of chromosome")
	parser.add_argument("--nhap", default=2,type=int,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", type=int,default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	args = parser.parse_args()

	'''
	print 'get TMRCA bin'
	TMRCA_bin,genome_percentage,rho=get_TMRCA_bin(args.file)
	print 'get chunk percentage'
	chunk_percentage,scale_factor=calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,COPY)
	print 'get TMRCA from macs'
#	genome_percentage_bin,chunk_percentage_bin,average_length_bin=get_TMRCA_from_macs_v3(args.file,TMRCA_bin,args.start)

	write_output_v2()

	print 'get TMRCA bin', args.file
	TMRCA_bin,genome_percentage,rho=get_TMRCA_bin(args.psmc)
#	print 'get chunk percentage'
#	chunk_percentage,scale_factor=calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,COPY)
	read_macs(args.file,TMRCA_bin)
	'''

	for psmc_file in glob.glob('*_0.psmc'):
		TMRCA_bin_ind,genome_percentage_ind,rho_ind=get_TMRCA_bin_single(psmc_file,1)
		chunk_percentage_ind,scale_factor_ind=calc_chunk_percentage_v2(genome_percentage_ind,TMRCA_bin_ind,rho_ind,1)
		psmc_file=psmc_file.replace(".psmc","")
		genome_percentage_bin,chunk_percentage_bin,average_length_bin=get_TMRCA_from_macs_v3(psmc_file,TMRCA_bin_ind,args.start,1)
		write_output_v6(psmc_file,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,scale_factor_ind,1,genome_percentage_bin,chunk_percentage_bin,average_length_bin)
#		write_output_v5(psmc_file,TMRCA_bin_ind,genome_percentage_ind,chunk_percentage_ind,1)
		print psmc_file
	