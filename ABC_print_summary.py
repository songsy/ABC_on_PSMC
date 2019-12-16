import argparse
from NGS_utils import *
import numpy as np
#from ete2 import Tree
import pickle
import math
import glob
import os
import subprocess
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys

COPY =4

def get_TMRCA_bin_real(file):
	m=['11','12','21','22']
	TMRCA_bin=[[] for x in range(COPY)]
	genome_percentage=[[] for x in range(COPY)]
	rho = [0 for x in range(COPY)]
	theta = [0 for x in range(COPY)]
	for i in range(COPY):
		f=open('%s_%s.psmc' %(file,m[i]),'r')
		start = False
		for line in f:
			line = line.strip()
			if line[:2]=='RD':
				number=int(line.split('\t')[1])
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
#		print TMRCA_bin[i],len(TMRCA_bin[i])
#		print genome_percentage[i],len(genome_percentage[i])
#		print theta[i],rho[i]
#		print np.sum(genome_percentage[i])
	return TMRCA_bin,genome_percentage,rho

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
#		print TMRCA_bin[i],len(TMRCA_bin[i])
#		print genome_percentage[i],len(genome_percentage[i])
#		print theta[i],rho[i]
#		print np.sum(genome_percentage[i])
	return TMRCA_bin,genome_percentage,rho

def get_TMRCA_bin_single_psmc(file):
	copy=1
	TMRCA_bin=[[] for x in range(copy)]
	genome_percentage=[[] for x in range(copy)]
	rho = [0 for x in range(copy)]
	theta = [0 for x in range(copy)]
	for i in range(copy):
		f=open(file,'r')
		start = False
		for line in f:
			line = line.strip()
			if line[:2]=='RD':
				number=int(line.split('\t')[1])
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
#		print TMRCA_bin[i],len(TMRCA_bin[i])
#		print genome_percentage[i],len(genome_percentage[i])
#		print theta[i],rho[i]
#		print np.sum(genome_percentage[i])
	return TMRCA_bin,genome_percentage,rho

def calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,copy):
	chunk_percentage=[[] for x in range(copy)]
	scale_factor=[[] for x in range(copy)]
	for index in range(copy):
		for i in range(len(genome_percentage[index])-1):
#			print TMRCA_bin[index][i]
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
		print len(genome_percentage[index]),genome_percentage[index]
		a=genome_percentage[index][len(genome_percentage[index])-1]*rho[index]*(TMRCA_bin[index][len(genome_percentage[index])-1])
		scale=1/(rho[index]*TMRCA_bin[index][len(genome_percentage[index])-1])
		chunk_percentage[index].append(a)
		scale_factor[index].append(scale)
	for index in range(copy):
		chunk_percentage_sum=sum(chunk_percentage[index])
		for i in range(len(chunk_percentage[index])):
			chunk_percentage[index][i]=chunk_percentage[index][i]/chunk_percentage_sum
	return chunk_percentage,scale_factor

def read_tm_list(iter):
	f=open('macs_iter%s_stats_weight.txt' %(iter),'r')
#	fout = open('macs_iter%s_stats.txt' %(iter),'w')
	tm_list = []
	diff_genome_list=[]
	diff_chunk_list=[]
	weight_list=[]
	for line in f:
		line=line.rstrip().split('\t')
		index,t,m,genome,chunk,weight=map(float,line)
		tm_list.append([float(t),float(m)])
		diff_chunk_list.append(chunk)
		weight_list.append(weight)
#		print t,m,chunk,weight
	return tm_list,diff_chunk_list,weight_list

def plot_stats(iter,new_tm_list,diff_chunk_list,weight_list):
	t=[i[0] for i in new_tm_list]
	m=[i[1] for i in new_tm_list]
	print len(t),len(m),len(diff_chunk_list),len(weight_list)
	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	ax1.scatter(t,m,c='blue')
	ax1.set_xlabel('split time($*10^4$ years)')
	ax1.set_ylabel('migration rate($*10^-5$)')
	fig.savefig("iter%s_sampled_param.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	cax=ax1.scatter(t,m,c=diff_chunk_list,cmap=plt.cm.coolwarm)
	ax1.set_xlabel('split time($*10^4$ years)')
	ax1.set_ylabel('migration rate($*10^-5$)')
	fig.colorbar(cax,orientation='vertical')
	fig.savefig("iter%s_sampled_param_summary_stats.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	cax=ax1.scatter(t,m,c=weight_list,cmap=plt.cm.coolwarm)
	ax1.set_xlabel('split time($*10^4$ years)')
	ax1.set_ylabel('migration rate($*10^-5$)')
	fig.colorbar(cax,orientation='vertical')
	fig.savefig("iter%s_sampled_param_weight.pdf" %(iter),format='pdf')
	return

def get_summary_stats(iter,tm_list,diff_chunk_list):
	fout = open('ABCestimator/iter%s_input_detail.txt' %(iter),'w')
	for i in range(len(tm_list)):
		t=tm_list[i][0]
		m=tm_list[i][1]
		file='iter%s/sim_macs_%s_%i_%s_%s.psmc.txt' %(iter,iter,i,t,m)
		data=map(lambda sline: sline.rstrip().split('\t'), open(file,'r'))
		a=[data[j][3] for j in [3,6,7,9,11]]
		print >>fout,'%s\t%s\t%s' %(t,m,'\t'.join(map(str,a)))

def write_output_v2():
	fout=open('ABCestimator/obs.txt','w')
	m=[]
	for i in [3,6,7,9,11]:
		a=[float(chunk_percentage[index][i])/sum(chunk_percentage[index]) for index in range(4)]
		m.append(np.mean(a))
	print >>fout,'\t'.join(map(str,m))
	fout.close()

def write_output(TMRCA_bin,genome_percentage,chunk_percentage,scale_factor):
	for index in range(1):
		fout=open(args.psmc+'_TMRCA_bin_%s.txt' %(index),'w')
		for i in range(len(TMRCA_bin[index])):
			print >>fout,'%f\t%f\t%f\t%f' %(TMRCA_bin[index][i]*30,genome_percentage[index][i],scale_factor[index][i],float(chunk_percentage[index][i])/sum(chunk_percentage[index]))
	fout.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='draw plot')
	parser.add_argument("--iter",type=int,default=0,dest='iter',help="iter")
	parser.add_argument("--psmc", dest='psmc',help="psmc file prefix")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--N0", default=5e04*2.36/1.25,dest='N0',help="N0")
	parser.add_argument("--simN0", default=5e04,dest='simN0',help="N0")
	args = parser.parse_args()


	tm_list,diff_chunk_list,weight_list=read_tm_list(args.iter)
	get_summary_stats(args.iter,tm_list,diff_chunk_list)
#	plot_stats(args.iter,tm_list,diff_chunk_list,weight_list)


	TMRCA_bin,genome_percentage,rho=get_TMRCA_bin_real(args.psmc)
	print 'get chunk percentage'

	chunk_percentage,scale_factor=calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,COPY)
	print 'get TMRCA from macs'

	write_output_v2()
	'''
	TMRCA_bin1,genome_percentage1,rho1=get_TMRCA_bin_single_psmc(args.psmc)
	chunk_percentage1,scale_factor1=calc_chunk_percentage(genome_percentage1,TMRCA_bin1,rho1,1)
	write_output(TMRCA_bin1,genome_percentage1,chunk_percentage1,scale_factor1)
	'''