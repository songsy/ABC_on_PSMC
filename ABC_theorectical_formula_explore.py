#!/usr/bin/env python
# python ABC_calc_TMRCA_v2.py
# Shiya Song
# 17th April 2015

import argparse,time,itertools,multiprocessing,math,glob,os,subprocess,pickle,sys
import numpy as np
from scipy import linalg

from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class TMRCA_BIN:
	def __init__(self):
		self.time_bin=[]
		self.genome_percentage=[]
		self.genome_percentage_macs=[]
		self.pop_size=[]
		self.rho=0
		self.theta=0
		self.covariance=0

	def read_from_psmc_mulitple(self,file,prefix,copy):
		if prefix is None:
			m=['11','12','21','22']
		else:
			m = prefix
		TMRCA_bin=[[] for x in range(copy)]
		genome_percentage=[[] for x in range(copy)]
		pop_size = [[] for x in range(copy)]
		rho = [0 for x in range(copy)]
		theta = [0 for x in range(copy)]
		for i in range(copy):
			if copy==1:
				f=open(file,'r')
			else:
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
						TMRCA_bin[i].append(float(line[2])*2*N0)
						genome_percentage[i].append(float(line[5]))
						pop_size[i].append(float(line[3])*N0)
					if line[:2]=='TR':
						line = line.split('\t')
						rho[i]=float(line[2])
						theta[i]=float(line[1])
		self.time_bin=list(np.mean(np.array(TMRCA_bin),axis=0))
		self.genome_percentage=list(np.mean(np.array(genome_percentage),axis=0))
		self.pop_size=list(np.mean(np.array(pop_size),axis=0))
		self.rho=np.mean(rho)
		self.theta=np.mean(theta)

	def read_from_msmc(self,file):
		f = open(file,'r')
		line = f.readline()
		for line in f:
			line=line.rstrip().split('\t')
			self.time_bin.append(abs(float(line[1]))/args.G)
			self.pop_size.append(abs(float(line[3]))*10000)
			self.genome_percentage.append(0)

	def print_TMRCA_bin(self):
		print self.theta,self.rho
		for i in range(len(self.time_bin)):
			print '%f\t%f\t%f' %(self.time_bin[i],self.genome_percentage[i],self.pop_size[i])

	def read_from_TMRCA_txt(self,file):
		f=open(file,'r')
		for line in f:
			line = line.rstrip().split('\t')
			self.time_bin.append(float(line[0]))
			self.genome_percentage.append(float(line[1]))
			self.genome_percentage_macs.append(float(line[4]))
		print self.time_bin
		print self.genome_percentage
		print self.genome_percentage_macs

	def calc_covariance(self):
		index=0
		genome_percentage_matrix=[]
		stats1=[]
		stats2=[]
		for prefix in ['11','12','21','22']:
			for i in range(100):
				genome_percentage_matrix.append([])
				f=open('/mnt/EXT/Kidd-scratch/shiya-projects/msmc/bootstrap_v2/%s/%s_true_fosmid_skip_unmask_%s_split_%s.psmc' %(args.pop,args.pop,prefix,i))
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
							genome_percentage_matrix[index].append(float(line[5]))
				index +=1
		genome_percentage_matrix_T=np.array(genome_percentage_matrix).T
#		print genome_percentage_matrix
		self.covariance=np.cov(genome_percentage_matrix_T)
		for i in range(len(genome_percentage_matrix)):
			stats1.append(calc_distance_between_two_dist(self.genome_percentage,genome_percentage_matrix[i],genome_percentage_matrix[i],self.covariance,'chi covariance'))
			stats2.append(calc_distance_between_two_dist(self.genome_percentage,genome_percentage_matrix[i],genome_percentage_matrix[i],self.covariance,'chi'))
		stats1=np.array(stats1)
		stats2=np.array(stats2)
		print 'chi covariance',np.mean(stats1),np.std(stats1),min(stats1),max(stats1),np.median(stats1),stats.percentileofscore(stats1,86)
		print 'chi',np.mean(stats2),np.std(stats2),min(stats2),max(stats2),np.median(stats2),stats.percentileofscore(stats2,0.018)
#		print self.covariance
		for i in range(len(self.covariance)):
			print i,self.covariance[i,i]

class time_bin:
	def __init__(self,time=0.,size=0.,theta1=0.,theta2=0.,theta_anc=0.,m1=0.,m2=0.,ancestral_integral=0.,recent_integral=0.,Q=0.,label=''):
		self.time=time
		self.theta1=theta1
		self.theta2=theta2
		self.theta_anc=theta_anc
		self.ancestral_integral=ancestral_integral
		self.recent_integral=recent_integral
		self.m1=m1
		self.m2=m2
		self.Q=Q
		self.label=label

	def print_time_bin(self):
		print '\t'.join(map(str,[self.time,self.theta1,self.theta2,self.theta_anc,self.m1,self.m2,self.label]))

def sample_parameters(Ntrial,tprior,mprior):
	tprior_list=np.random.uniform(tprior[0],tprior[1],Ntrial)
	mprior_list=np.random.uniform(mprior[0],mprior[1],Ntrial)
	tm_list=zip(tprior_list,mprior_list)
	return tm_list

def calc_f_density(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
	NA=2e04
	NB=1e04
	N_Anc=1e04
	theta1=4*NA*args.mu;
	theta2=4*NB*args.mu;
	theta_Anc=4*N_Anc*args.mu;
	T=t*1e04/args.G*args.mu
	t_bin=np.arange(0,0.004,0.000001)
	Q=np.array([[-2*m1-2/theta1,2*m1,0,2/theta1,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta2,0,2/theta2],[0,0,0,-m1,m1],[0,0,0,m2,-m2]]);
	f_density=np.zeros(len(t_bin))
	M=linalg.expm(Q*T);
	for i in range(len(t_bin)):
	    if t_bin[i]<T:
	        m=linalg.expm(Q*t_bin[i]);
	        f_density[i]=m[1,0]*2/theta1+m[1,2]*2/theta2;
	    else:
	        f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc*np.exp(-(2/theta_Anc)*(t_bin[i]-T));
	return f_density

def calc_f_density_Gravel(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
	NA=14474*2.36/1.25
	NB=1861*2.36/1.25
	N_Anc1=7310*2.36/1.25
	N_Anc2=14474*2.36/1.25
	theta1=4*NA*args.mu;
	theta2=4*NB*args.mu;
	theta_Anc1=4*N_Anc1*args.mu;
	theta_Anc2=4*N_Anc2*args.mu;
	t_change=335300/args.G*args.mu
	T=t*1e04/args.G*args.mu
	t_bin=np.arange(0,0.004,0.000001)
	Q=np.array([[-2*m1-2/theta1,2*m1,0,2/theta1,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta2,0,2/theta2],[0,0,0,-m1,m1],[0,0,0,m2,-m2]]);
	f_density=np.zeros(len(t_bin))
	M=linalg.expm(Q*T);
	for i in range(len(t_bin)):
	    if t_bin[i]<T:
	        m=linalg.expm(Q*t_bin[i]);
	        f_density[i]=m[1,0]*2/theta1+m[1,2]*2/theta2;
	    elif t_bin[i]<t_change:
			f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc2*np.exp(-(2/theta_Anc2)*(t_bin[i]-T));
	    else:
	        f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc1*np.exp(-(2/theta_Anc2)*(t_change-T)-(2/theta_Anc1)*(t_bin[i]-t_change));
	return f_density

def calc_f_density_complex_v5(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12]
	elif args.type=='msmc':
		interval_Af = [0,1,2,3,4,5,6]
	migration_done=False
	if migration_end==0:
		migration_done=True
	print t,m
	Af_index=0
	Eu_index=0
	current_t=0
	while True:
		print Af_index,Eu_index,current_t
		if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]==TMRCA_bin_Eu.time_bin[interval_Af[Eu_index]]:
			if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=float(migration_end*1e04):
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,label='after split'))
		else:
			if current_t*args.G>=float(migration_end*1e04):
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index-1 if Af_index>=1 else Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index-1 if Eu_index>=1 else Eu_index]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if current_t*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=current_t*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				if current_t*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=current_t*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,label='after split'))
		if TMRCA_bin_Af.time_bin[interval_Af[Af_index+1]]<TMRCA_bin_Eu.time_bin[interval_Af[Eu_index+1]]:
			Af_index+=1
			current_t=TMRCA_bin_Af.time_bin[interval_Af[Af_index]]
		else:
			Eu_index+=1
			current_t=TMRCA_bin_Eu.time_bin[interval_Af[Eu_index]]

	'''
	for i in range(len(interval_Af)):
		if TMRCA_bin_Af.time_bin[interval_Af[i]]==TMRCA_bin_Eu.time_bin[interval_Af[i]]:
			print i,'case1'
			if TMRCA_bin_Af.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,label='after split'))
		elif TMRCA_bin_Af.time_bin[interval_Af[i]]<TMRCA_bin_Eu.time_bin[interval_Af[i]]:
			print i,'case2'
			if TMRCA_bin_Af.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,label='after split'))
			if TMRCA_bin_Eu.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,label='after split'))

		else:
			print i,'case3'
			if TMRCA_bin_Eu.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
			else:
				if TMRCA_bin_Af.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,label='after split'))
			if TMRCA_bin_Af.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))	
			else:
				if TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,label='after split'))	
	'''
	if TMRCA_bin_psuedo.time_bin[checkpoint]*args.G>=t*1e04:
		Anc_interval=checkpoint
	else:
		for i in range(len(TMRCA_bin_psuedo.time_bin)):
			if TMRCA_bin_psuedo.time_bin[i]*args.G<=t*1e04 and TMRCA_bin_psuedo.time_bin[i+1]*args.G>=t*1e04:
				Anc_interval=i
				break
	theta_Anc=4*TMRCA_bin_psuedo.pop_size[Anc_interval]*args.mu
	time_bin_list.append(time_bin(time=float(t*1e04)/args.G*args.mu,theta_anc=theta_Anc,label='split time'))
	first = True
	prev_size=TMRCA_bin_psuedo.pop_size[Anc_interval+1]
	for i in range(Anc_interval+1,len(TMRCA_bin_psuedo.time_bin)):
		if first is True:
			first= False
		else:
			if prev_size==TMRCA_bin_psuedo.pop_size[i]:
				continue
		prev_size=TMRCA_bin_psuedo.pop_size[i]
		time_bin_list.append(time_bin(time=TMRCA_bin_psuedo.time_bin[i]*args.mu,theta_anc=4*TMRCA_bin_psuedo.pop_size[i]*args.mu,label='before split'))
	t_bin=np.arange(0,0.0065,0.000001)
	f_density=np.zeros(len(t_bin))
	f_density11=np.zeros(len(t_bin))
	f_density22=np.zeros(len(t_bin))
	for i in range(len(time_bin_list)):
		if time_bin_list[i].label in ['after split','migration_end']:
			m1=time_bin_list[i].m1
			m2=time_bin_list[i].m2
			theta_Af=time_bin_list[i].theta1
			theta_Eu=time_bin_list[i].theta2
			time_bin_list[i].Q=np.mat([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]],dtype=np.float64);
#			print i,"Q matrix:",sum(sum(np.array(time_bin_list[i].Q)))
#			print time_bin_list[i].Q
			if i==0:
				time_bin_list[i].recent_integral=1
			else:
				time_bin_list[i].recent_integral=time_bin_list[i-1].recent_integral*np.mat(linalg.expm(time_bin_list[i-1].Q*(time_bin_list[i].time-time_bin_list[i-1].time)))
#				print time_bin_list[i].recent_integral
#				print i,'matrix sum',sum(sum(np.array(time_bin_list[i].recent_integral)))
		elif time_bin_list[i].label=='split time':
			time_bin_list[i].recent_integral=time_bin_list[i-1].recent_integral*np.mat(linalg.expm(time_bin_list[i-1].Q*(time_bin_list[i].time-time_bin_list[i-1].time)))
			M=time_bin_list[i].recent_integral
#			print 'M:',sum(sum(M))
			time_bin_list[i].ancestral_integral=0
		else:
			time_bin_list[i].ancestral_integral=time_bin_list[i-1].ancestral_integral-(2./time_bin_list[i-1].theta_anc)*(time_bin_list[i].time-time_bin_list[i-1].time)
	time_bin_list.append(time_bin(time=0.0065))
#	print t,m,Anc_interval
	for i in range(1,len(time_bin_list)):
		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong'

	T=float(t*1e04)/args.G*args.mu
#	print len(t_bin),T
	for i in range(len(t_bin)):
		if t_bin[i]<T:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time:
					m_final=time_bin_list[j].recent_integral*np.mat(linalg.expm(time_bin_list[j].Q*(t_bin[i]-time_bin_list[j].time)))
#					print i,t_bin[i],f_density[i]
					f_density[i]=m_final[1,0]*2./time_bin_list[j].theta1+m_final[1,2]*2./time_bin_list[j].theta2
					f_density11[i]=m_final[0,0]*2./time_bin_list[j].theta1+m_final[0,2]*2./time_bin_list[j].theta2
					f_density22[i]=m_final[2,0]*2./time_bin_list[j].theta1+m_final[2,2]*2./time_bin_list[j].theta2
					break
		else:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time and time_bin_list[j].label in ['split time','before split']:
#					if time_bin_list[j].theta_anc==0:
#						print t_bin[i],j,time_bin_list[j].theta_anc
					f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
					f_density11[i]=(M[0,0]+M[0,1]+M[0,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
					f_density22[i]=(M[2,0]+M[2,1]+M[2,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					print i,t_bin[i],f_density[i]
					break
#		print i,t_bin[i],f_density[i]
#	print sum(f_density),sum(f_density)*0.000001
	return f_density,f_density11,f_density22

def calc_f_density_integrate_single(t,m,migration_end,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density,f_density11,f_density22=calc_f_density_complex_v5(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P11=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P22=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				P11[j]=P11[j]+f_density11[i]/sum(f_density11);
				P22[j]=P22[j]+f_density22[i]/sum(f_density22);
				break
	return [f_density,f_density11,f_density22],[P,P11,P22]

def draw_f_density_compare_v3(t,m,migration_end,f_density,P,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy')
	color=['r','g','b','y','c']
#	label=['70k,20','80k,20','90k,20','100k,20']
	label=['pop 1-2','pop 1','pop 2']
	for index in range(len(f_density)):
		ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density[index],color[index])
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax2 = fig.add_subplot(312,title='TMRCA distribution')
	for index in range(len(P)):
		ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P[index],color[index],label=label[index])
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'k',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)
	ax3 = fig.add_subplot(313,title='population size(in 10k)')
	ax3.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r',label='pop 1-2')
	ax3.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],[float(i)/1e04 for i in TMRCA_bin_Af.pop_size],'b',label='pop 1')
	ax3.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],[float(i)/1e04 for i in TMRCA_bin_Eu.pop_size],'g',label='pop 2')
	ax3.legend(prop={'size':10})
	ax3.grid(True)
	ax3.set_ylim([0,10])
	fig.savefig("sim_YRI_CEU_%s_%s_%s_f_density_compare.pdf" %(t,m,migration_end),format='pdf')

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--psmc", dest='psmc',help="psmc file prefix")  
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--Ntrial", default=1000,type=int,dest='Ntrial',help="Number of trials per iteration")
	parser.add_argument("--Nparticle", default=500,type=int,dest='Nparticle',help="number of particles to reach")
	parser.add_argument("--initial_e",type=float,nargs='+',dest='initial_e',help="inital error to tolerance")
	parser.add_argument("--nchr", default=100,type=int,dest='nchr',help="number of chromosome per simulation")
	parser.add_argument("--nhap", default=2,type=int,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", type=int,default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--G", default=30,dest='G',help="Generation time")
	parser.add_argument("--ratio", type=float,default=0.2,dest='ratio',help="ratio of r to mu when simulating using macs")
	parser.add_argument("--tprior1", type=int,default=8,dest='tprior1',help="split time prior")
	parser.add_argument("--tprior2", type=int,default=15,dest='tprior2',help="split time prior")
	parser.add_argument("--mprior1", type=int,default=0,dest='mprior1',help="migration prior")
	parser.add_argument("--mprior2", type=int,default=20,dest='mprior2',help="migration prior")
	parser.add_argument("--kernel",default="uniform",dest="kernel",help="choose from gaussian or uniform kernel")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	parser.add_argument("--pop", dest='pop',default='YRI_CEU',help="YRI_CEU")
	parser.add_argument("--type", dest='type',default='psmc',help="psmc or msmc input for file1 and file2")
	parser.add_argument("--version", dest='version',help="version")
	parser.add_argument("--checkpoint", type=int,default=1e05,dest='checkpoint',help="checkpoint in yr")
	parser.add_argument("--migration_end", type=float, default=0,dest='migration_end',help="migration end in 10kyr")
	parser.add_argument("--interval_compare", type=int,default=13,dest='interval_compare',help="number of first intervals to compare")
	parser.add_argument("--last_iter", type=int,default=0,dest='last_iter',help="last_iter")
	parser.add_argument("--log", dest='log',help="logfile")
	args = parser.parse_args()

#	sys.stdout = open(args.log,'a',0)

#	os.popen('mkdir -p version%s' %(args.version))
	ITERATION=3

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'

	TMRCA_bin_pseudo=TMRCA_BIN()   # the pseudo diploid based on real data(the one you observe)
	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,None,4)
	TMRCA_bin_pseudo.calc_covariance()
	TMRCA_bin_pseudo.print_TMRCA_bin()
	TMRCA_bin_Af=TMRCA_BIN()
	if args.type=='psmc':
		TMRCA_bin_Af.read_from_psmc_mulitple(args.file1,None,1)
	elif args.type=='msmc':
		TMRCA_bin_Af.read_from_msmc(args.file1)
	TMRCA_bin_Af.print_TMRCA_bin()
	TMRCA_bin_Eu=TMRCA_BIN()
	if args.type=='psmc':
		TMRCA_bin_Eu.read_from_psmc_mulitple(args.file2,None,1)
	elif args.type=='msmc':
		TMRCA_bin_Eu.read_from_msmc(args.file2)
	TMRCA_bin_Eu.print_TMRCA_bin()


	# Test function
	f_density,P=calc_f_density_integrate_single(6.17,7.9,2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print P1
#	print P2
#	draw_f_density_compare(6.24,987,args.migration_end,f_density1,P1,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	draw_f_density_compare_v3(6.17,7.9,2,f_density,P,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)

