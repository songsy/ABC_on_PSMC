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
	ax1.set_xlim([7, 16])
	ax1.set_ylim([-5, 25])
	fig.savefig("iter%s_sampled_param.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	cax=ax1.scatter(t,m,c=diff_chunk_list,cmap=plt.cm.coolwarm)
	ax1.set_xlabel('split time($*10^4$ years)')
	ax1.set_ylabel('migration rate($*10^-5$)')
	ax1.set_xlim([7, 16])
	ax1.set_ylim([-5, 25])
	fig.colorbar(cax,orientation='vertical')
	fig.savefig("iter%s_sampled_param_summary_stats.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311)
	cax=ax1.scatter(t,m,c=weight_list,cmap=plt.cm.coolwarm)
	ax1.set_xlabel('split time($*10^4$ years)')
	ax1.set_ylabel('migration rate($*10^-5$)')
	ax1.set_xlim([7, 16])
	ax1.set_ylim([-5, 25])
	fig.colorbar(cax,orientation='vertical')
	fig.savefig("iter%s_sampled_param_weight.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(211,title='iter %s stats' %(iter))
	
	bins = 40
	ax1.hist(diff_chunk_list,facecolor='blue',alpha=0.5)
	ax1.set_xlabel('distance between summary statistics')
	ax1.set_ylabel('count')
	ax1.set_xlim([0, 0.03])
	if iter<=1:
		e = np.percentile(diff_chunk_list,50)
	else:
		e = np.percentile(diff_chunk_list,20)
	ax1.axvline(x=e,linewidth=2,color='k')
	accept_indicator=[1 if i<=e else 0 for i in diff_chunk_list]
	fig.savefig("iter%s_summary_stats_distance.pdf" %(iter),format='pdf')
	return

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='draw plot')
	parser.add_argument("--iter",type=int,default=0,dest='iter',help="iter")
	args = parser.parse_args()

	tm_list,diff_chunk_list,weight_list=read_tm_list(args.iter)
	plot_stats(args.iter,tm_list,diff_chunk_list,weight_list)


