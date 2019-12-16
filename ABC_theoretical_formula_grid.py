#!/usr/bin/env python
# python ABC_calc_TMRCA_v2.py
# Shiya Song
# 17th April 2015

import argparse,time,itertools,multiprocessing,math,glob,os,subprocess,pickle,sys
from NGS_utils import *
import numpy as np
from scipy import linalg
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

	def add_macs(self,file):
		f=open(file,'r')
		for line in f:
			line = line.rstrip().split('\t')
			self.genome_percentage_macs.append(float(line[4]))

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
	t_bin=np.arange(0,0.0065,0.000001)
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
	t_bin=np.arange(0,0.0065,0.000001)
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

def calc_f_density_Gravel_complex(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint):
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
#	NA=14474*2.36/1.25
#	NB=1861*2.36/1.25
	NA=40000
	NB=5000
	theta_Af=4*args.mu*NA
	theta_Eu=4*args.mu*NB
	Anc_interval = 0
	for i in range(len(TMRCA_bin_psuedo.time_bin)):
		if TMRCA_bin_psuedo.time_bin[i]*args.G>=max(checkpoint,t*1e04):
			Anc_interval=i
			break
	N_Anc=[TMRCA_bin_psuedo.pop_size[i] for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	t_change=[TMRCA_bin_psuedo.time_bin[i]*args.mu for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	theta_Anc=[4*i*args.mu for i in N_Anc]
	print 'N_Anc:',N_Anc
	print 'theta_Anc:',theta_Anc
	print 't_change:',t_change
	T=float(t*1e04)/args.G*args.mu
	print 'T',T
	t_bin=np.arange(0,0.0065,0.000001)
	Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]]);
	f_density=np.zeros(len(t_bin))
	M=linalg.expm(Q*T);
	ancestral_integral=np.zeros(len(N_Anc))
	for i in range(1,len(N_Anc)):
		if t_change[i-1]<T:
			continue
		if i==1:
			ancestral_integral[i]=-(2/theta_Anc[i-1])*(t_change[i-1]-T)
		else:
			ancestral_integral[i]=ancestral_integral[i-1]-(2/theta_Anc[i-2])*(t_change[i-1]-t_change[i-2])
	for i in range(len(t_bin)):
	    if t_bin[i]<T:
	        m=linalg.expm(Q*t_bin[i]);
	        f_density[i]=m[1,0]*2/theta_Af+m[1,2]*2/theta_Eu;
	    else:
	    	for j in range(len(N_Anc)):
	    		if t_bin[i]<t_change[j]:
	    			if j==0:
	    				f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[j]*np.exp(-(2/theta_Anc[j])*(t_bin[i]-T))
	    			else:
	    				f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[j-1]*np.exp(ancestral_integral[j]-(2/theta_Anc[j-1])*(t_bin[i]-t_change[j-1]))
	    			break
	    	if t_bin[i]>t_change[-1]:
	    		f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[-1]*np.exp(ancestral_integral[-1]-(2/theta_Anc[-1])*(t_bin[i]-t_change[-1]));
	print sum(f_density),sum(f_density)*0.000001
	return f_density

def calc_f_density_complex(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint):
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.2
	macs_cmds='macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
	interval_Af = []
	for i in range(0,len(TMRCA_bin_Af.time_bin)):
		if TMRCA_bin_Af.time_bin[i]*args.G<checkpoint:
			interval_Af.append(i)
		else:
			break
	NA=np.mean([TMRCA_bin_Af.pop_size[i] for i in interval_Af])
	theta_Af=4*args.mu*NA
	interval_Eu=[]
	for i in range(0,len(TMRCA_bin_Af.time_bin)):
		if TMRCA_bin_Eu.time_bin[i]*args.G<checkpoint:
			interval_Eu.append(i)
		else:
			break
	NB=TMRCA_bin_Eu.pop_size[interval_Eu[-2]]
#	NB=mean([TMRCA_bin_Eu.pop_size[i] for i in interval_Eu])
	theta_Eu=4*args.mu*NB
	macs_cmds+='-n 1 %f -n 2 %f ' %(NA/args.N0,NB/args.N0)
	macs_cmds+='-ej %f 2 1 ' %(float(t*1e04)/args.G/4/args.N0)
#	print 'N and theta:',NA,NB,theta_Af,theta_Eu
	Anc_interval = 0
	for i in range(len(TMRCA_bin_psuedo.time_bin)):
		if TMRCA_bin_psuedo.time_bin[i]*args.G>=max(checkpoint,t*1e04):
			Anc_interval=i
			break
	N_Anc=[TMRCA_bin_psuedo.pop_size[i] for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	t_change=[TMRCA_bin_psuedo.time_bin[i]*args.mu for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	theta_Anc=[4*i*args.mu for i in N_Anc]
	macs_cmds += '-en %f 1 %f ' %(float(t*1e04)/args.G/4/args.N0+0.000001,N_Anc[0]/args.N0)
	for i in range(len(N_Anc)):
		if i!=0:
			if N_Anc[i]==N_Anc[i-1]:
				continue
		macs_cmds+='-en %f 1 %f ' %(t_change[i]/args.mu/4/args.N0,N_Anc[i]/args.N0)
	macs_cmds += '-ema 0.000000 2 0 %f %f 0 ' %(m/1e05*4*int(args.N0),m/1e05*4*int(args.N0))   # migration starting from 15kyr
	macs_cmds += '-eM %f 0 ' %(float(t*1e04)/args.G/4/args.N0-0.000001)
#	print macs_cmds
	print 'N_Anc:',N_Anc
	print 'theta_Anc:',theta_Anc
	print 't_change:',t_change
	T=float(t*1e04)/args.G*args.mu
	print 'T',T
	t_bin=np.arange(0,0.0065,0.000001)
	Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]]);
	f_density=np.zeros(len(t_bin))
	M=linalg.expm(Q*T);
	ancestral_integral=np.zeros(len(N_Anc))
	for i in range(1,len(N_Anc)):
		if t_change[i-1]<T:
			continue
		if i==1:
			ancestral_integral[i]=-(2/theta_Anc[i-1])*(t_change[i-1]-T)
		else:
			ancestral_integral[i]=ancestral_integral[i-1]-(2/theta_Anc[i-2])*(t_change[i-1]-t_change[i-2])
	for i in range(len(t_bin)):
	    if t_bin[i]<T:
	        m=linalg.expm(Q*t_bin[i]);
	        f_density[i]=m[1,0]*2/theta_Af+m[1,2]*2/theta_Eu;
	    else:
	    	for j in range(len(N_Anc)):
	    		if t_bin[i]<t_change[j]:
	    			if j==0:
	    				f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[j]*np.exp(-(2/theta_Anc[j])*(t_bin[i]-T))
	    			else:
	    				f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[j-1]*np.exp(ancestral_integral[j]-(2/theta_Anc[j-1])*(t_bin[i]-t_change[j-1]))
	    			break
	    	if t_bin[i]>t_change[-1]:
	    		f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/theta_Anc[-1]*np.exp(ancestral_integral[-1]-(2/theta_Anc[-1])*(t_bin[i]-t_change[-1]));
	print sum(f_density),sum(f_density)*0.000001
	return f_density

def calc_f_density_complex_v2(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint):
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
	interval_Af = []
	for i in range(0,len(TMRCA_bin_Af.time_bin)):
		if TMRCA_bin_Af.time_bin[i]*args.G<checkpoint:
			interval_Af.append(i)
		else:
			break
	interval_Eu = []
	for i in range(0,len(TMRCA_bin_Af.time_bin)):
		if TMRCA_bin_Eu.time_bin[i]*args.G<checkpoint:
			interval_Eu.append(i)
		else:
			break
	NB=mean([TMRCA_bin_Eu.pop_size[i] for i in interval_Eu])
	interval_after_split=[TMRCA_bin_Eu.time_bin[i] for i in interval_Eu]+[TMRCA_bin_Af.time_bin[i] for i in interval_Af]
	interval_after_split.sort()

	theta_Af=np.zeros(len(interval_after_split))
	for i in range(len(TMRCA_bin_Af.time_bin)):
		for j in range(len(interval_after_split)-1):
			if TMRCA_bin_Af.time_bin[i]>=interval_after_split[j] and TMRCA_bin_Af.time_bin[i]<=interval_after_split[j+1]:
				theta_Af[j]=4*TMRCA_bin_Af.pop_size[i]*args.mu
				break
		if TMRCA_bin_Af.time_bin[i]>=interval_after_split[-1]:
			theta_Af[len(interval_after_split)]=4*TMRCA_bin_Af.pop_size[i]*args.mu
	theta_Eu=np.zeros(len(interval_after_split))
	for i in range(len(TMRCA_bin_Eu.time_bin)):
		for j in range(len(interval_after_split)-1):
			if TMRCA_bin_Eu.time_bin[i]>=interval_after_split[j] and TMRCA_bin_Eu.time_bin[i]<=interval_after_split[j+1]:
				theta_Eu[j]=4*TMRCA_bin_Eu.pop_size[i]*args.mu
				break
		if TMRCA_bin_Eu.time_bin[i]>=interval_after_split[-1]:
			theta_Eu[len(interval_after_split)]=4*TMRCA_bin_Eu.pop_size[i]*args.mu
	interval_after_split=[i*args.mu for i in interval_after_split]
	print 'interval after split',interval_after_split
	
	NA=mean([TMRCA_bin_Af.pop_size[i] for i in interval_Af])
	Anc_interval = 0
	for i in range(len(TMRCA_bin_Af.time_bin)):
		if TMRCA_bin_Eu.time_bin[i]*args.G>=checkpoint:
			Anc_interval=i
			break
	N_Anc=[TMRCA_bin_psuedo.pop_size[i] for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	t_change=[TMRCA_bin_psuedo.time_bin[i]*args.mu for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	theta1=4*NA*args.mu;
	theta2=4*NB*args.mu;
	theta_Anc=[4*i*args.mu for i in N_Anc]
	T=t*1e04/args.G*args.mu
	t_bin=np.arange(0,0.0065,0.000001)
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

def calc_distance_between_two_dist(A,B,C,method):
	stats=0
	if method=='chi':   # chi-square
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=pow((A[i]-B[i]),2)/A[i]
	if method=='total variation dist':
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=abs(A[i]-B[i])
	if method=='chi omit':   # chi-square
		for i in range(len(A)):
			if A[i]==0:
				continue
			if i in [10,11]:
				continue
			stats+=pow((A[i]-B[i]),2)/A[i]
	if method=='chi scaled':
		for i in range(len(A)):
			if A[i]==0:
				continue
			if abs(A[i]-C[i])>=0.005:
				scale=0.2
			elif abs(A[i]-C[i])>=0.001:
				scale=0.5
			else:
				scale=1
			stats+=pow((A[i]-B[i]),2)/A[i]*scale
	return stats

def calc_f_density_integrate(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density=calc_f_density(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,'chi')
	return stats

def calc_f_density_integrate_v2(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density=calc_f_density(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi scaled')
	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage_macs,P,TMRCA_bin_psuedo.genome_percentage,'chi scaled')
	return (stats1,stats2)

def calc_f_density_integrate_v3(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density=calc_f_density_Gravel_complex(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi scaled')
	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage_macs,P,TMRCA_bin_psuedo.genome_percentage,'chi scaled')
	stats3=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi')
	stats4=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage_macs,P,TMRCA_bin_psuedo.genome_percentage,'chi')
	return (P,stats1,stats2,stats3,stats4)

def calc_f_density_integrate_single(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density=calc_f_density_Gravel_complex(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	return f_density,P

def calc_f_density_integrate_single_v2(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density=calc_f_density_Gravel(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	return f_density,P

def calc_f_density_integrate_wrapper(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate(*args2)
	return stats

def calc_f_density_integrate_wrapper_v2(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate_v2(*args2)
	return stats

def calc_f_density_integrate_wrapper_v3(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate_v3(*args2)
	return stats

def draw_f_density(f_density):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	plt.semilogx(t_bin,f_density)
	plt.title('TMRCA density')
	plt.grid(True)
	fig.savefig("sim_macs_YRI_CEU_fosmid_f_density_v2.pdf",format='pdf')


def draw_f_density_compare(f_density1,f_density2,P1,P2,TMRCA_bin1,TMRCA_bin2):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy',label='the one from psmc from simulated data')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density1,'r',label='the one from psmc from simulated data')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density2,'b',label='the one from data input into simulation')
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax1.legend(prop={'size':10})
	ax2 = fig.add_subplot(312,title='TMRCA distribution',)
	ax2.semilogx([i*30 for i in TMRCA_bin1.time_bin],P1,'r',label='the one from psmc from simulated data')
	ax2.semilogx([i*30 for i in TMRCA_bin2.time_bin],P2,'b',label='the one from data input into simulation')
	ax2.grid(True)
	ax3 = fig.add_subplot(313,title='population size(in 10k)',label='the one from data input into simulation')
	ax3.semilogx([i*30 for i in TMRCA_bin1.time_bin],[float(i)/1e04 for i in TMRCA_bin1.pop_size],'r',label='the one from psmc from simulated data')
	ax3.semilogx([i*30 for i in TMRCA_bin2.time_bin],[float(i)/1e04 for i in TMRCA_bin1.pop_size],'b',label='the one from data input into simulation')
	ax3.grid(True)
	fig.savefig("sim_macs_Gravel_fosmid_f_density.pdf",format='pdf')

def explore_stats_in_grid(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for t in range(8,15):
		for m in range(0,51,5):
			tm_list.append([t,m])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v2,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	psmc_result=[i[0] for i in result1]
	macs_result=[i[1] for i in result1]
	print len(psmc_result),psmc_result
	print len(macs_result),macs_result
	psmc_result_matrix=np.reshape(np.array(psmc_result),(7,11))
	macs_result_matrix=np.reshape(np.array(macs_result),(7,11))
	print psmc_result_matrix
	print macs_result_matrix
	return psmc_result_matrix,macs_result_matrix

def explore_stats_in_grid_verbose(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for t in range(8,15):
		for m in range(0,51,5):
			tm_list.append([t,m])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v3,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	TMRCA_dist=[i[0] for i in result1]
	psmc_result_scaled=[i[1] for i in result1]
	macs_result_scaled=[i[2] for i in result1]
	psmc_result=[i[3] for i in result1]
	macs_result=[i[4] for i in result1]
#	print len(psmc_result),psmc_result
#	print len(macs_result),macs_result
	psmc_result_matrix_scaled=np.reshape(np.array(psmc_result_scaled),(7,11))
	macs_result_matrix_scaled=np.reshape(np.array(macs_result_scaled),(7,11))
	psmc_result_matrix=np.reshape(np.array(psmc_result),(7,11))
	macs_result_matrix=np.reshape(np.array(macs_result),(7,11))
	print 'psmc chi',psmc_result_matrix
	print 'macs chi',macs_result_matrix
	print 'psmc chi scaled',psmc_result_matrix_scaled
	print 'macs chi scaled',macs_result_matrix_scaled
	return psmc_result_matrix_scaled,macs_result_matrix_scaled,psmc_result_matrix,macs_result_matrix,TMRCA_dist,tm_list

def explore_stats_in_grid_Gravel(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for t in range(8,15):
		for m in range(0,51,5):
			tm_list.append([t,m])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v2,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	psmc_result=[i[0] for i in result1]
	macs_result=[i[1] for i in result1]
	print len(psmc_result),psmc_result
	print len(macs_result),macs_result
	psmc_result_matrix=np.reshape(np.array(psmc_result),(7,11))
	macs_result_matrix=np.reshape(np.array(macs_result),(7,11))
	print psmc_result_matrix
	print macs_result_matrix
	return psmc_result_matrix,macs_result_matrix

def draw_result_matrix(psmc_result_matrix,macs_result_matrix,title):
	fig = plt.figure()
	ax1 = fig.add_subplot(211,title='chi square stats based on PSMC')
	bins = 40
	cax=ax1.imshow(psmc_result_matrix.transpose(),interpolation="None",extent=[6,13,50,0],aspect="auto",)
	cax.set_clim(np.min(list(psmc_result_matrix)),np.max(psmc_result_matrix))
	print np.min(list(psmc_result_matrix)),np.max(psmc_result_matrix)
	ax1.legend(prop={'size':6})
	ax1.set_ylabel('migration rate')
	fig.colorbar(cax,orientation='vertical')

	ax2 = fig.add_subplot(212,title='chi square stats based on macs')
	cax=ax2.imshow(macs_result_matrix.transpose(),interpolation="None",extent=[6,13,50,0],aspect="auto",)
	cax.set_clim(np.min(macs_result_matrix),np.max(macs_result_matrix))
	print np.min(macs_result_matrix),np.max(macs_result_matrix)
	ax2.legend(prop={'size':6})
	ax2.set_xlabel('split time')
	ax2.set_ylabel('migration rate')
	fig.colorbar(cax,orientation='vertical')

	fig.savefig("%s_stats_psmc_vs_macs.pdf" %(title),format='pdf')

def plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight):
	t=[i[0] for i in tm_list_total]
	m=[i[1] for i in tm_list_total]
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior1,args.mprior2,num=(args.mprior2-args.mprior1)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t,m,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s sampled param' %(iter))
	cax=ax1.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax1.set_ylabel('migration rate')
	t_accepted=[i[0] for i in tm_list_new]
	m_accepted=[i[1] for i in tm_list_new]
	ax2 = fig.add_subplot(312,title='iter %s accepted param' %(iter))
	heatmap, xedges, yedges = np.histogram2d(t_accepted, m_accepted,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax2.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax2.set_ylabel('migration rate')
	ax3 = fig.add_subplot(313)
	ax3.set_xlim([args.tprior1,args.tprior2])
	ax3.set_ylim([args.mprior1,args.mprior2])
	ax3.set_xlabel('split time')
	ax3.plot(t,m,'ro')
	ax3.plot(t_accepted,m_accepted,'bo')
	fig.savefig("iter%s_sampled_accepted_param.pdf" %(iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s stats' %(iter))
	bins = 40
	ax1.hist(result_total,facecolor='green',alpha=0.5,histtype='stepfilled',label='genome percentage')
	ax1.legend(prop={'size':10})
	ax1.set_xlabel('chi square stats based on genome percentage')
	ax1.set_ylabel('count')
	ax1.axvline(x=e,linewidth=2,color='k')

	ax2 = fig.add_subplot(312)
	cax=ax2.scatter(t_accepted,m_accepted,c=result_new,cmap=plt.cm.coolwarm)
	ax2.set_xlabel('chi square stats based on genome percentage')
	fig.colorbar(cax,orientation='vertical')

	if iter==0:
		new_weight = [1./len(tm_list_new) for i in range(len(tm_list_new))]
	else:
		t_prev=[i[0] for i in prev_tm_list]
		m_prev=[i[1] for i in prev_tm_list]
		if args.kernel == "uniform":
			sigma_t=0.5*(max(t_prev)-min(t_prev))
			sigma_m=0.5*(max(m_prev)-min(m_prev))
			print 'sigma_t', sigma_t
			print 'sigma_m', sigma_m
			new_weight=update_weight(iter,prev_tm_list,tm_list_new,prev_weight,sigma_t,sigma_m)
	ax3 = fig.add_subplot(313)
	cax=ax3.scatter(t_accepted,m_accepted,c=new_weight,cmap=plt.cm.coolwarm)
	ax3.set_xlabel('weight')
	fig.colorbar(cax,orientation='vertical')	
	fig.savefig("iter%s_distance_stats_scatterplot.pdf" %(iter),format='pdf')

	fout = open('iter%s_stats_weight.txt' %(iter),'w')
	for i in range(len(new_weight)):
		t=tm_list_new[i][0]
		m=tm_list_new[i][1]
		print >>fout,'%i\t%s\t%s\t%f\t%f' %(i,t,m,result_new[i],new_weight[i])
	fout.close()
	return new_weight

def print_TMRCA_result(tm_list,TMRCA_dist):
	f=open('sim_macs_Gravel_grid_TMRCA_v2.txt','w')
	for i in range(len(tm_list)):
		print >>f,"%s\t%s\t%s" %(tm_list[i][0],tm_list[i][1],"\t".join(map(str,TMRCA_dist[i])))

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--psmc", dest='psmc',help="psmc file prefix")  
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--mu", type=float,default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--Ntrial", default=1000,type=int,dest='Ntrial',help="Number of trials per iteration")
	parser.add_argument("--Nparticle", default=500,type=int,dest='Nparticle',help="number of particles to reach")
	parser.add_argument("--initial_e",type=float,nargs='+',dest='initial_e',help="inital error to tolerance")
	parser.add_argument("--nchr", default=100,type=int,dest='nchr',help="number of chromosome per simulation")
	parser.add_argument("--nhap", default=2,type=int,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", type=int,default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--G", type=int,default=30,dest='G',help="Generation time")
	parser.add_argument("--ratio", type=float,default=0.2,dest='ratio',help="ratio of r to mu when simulating using macs")
	parser.add_argument("--tprior1", type=int,default=8,dest='tprior1',help="split time prior")
	parser.add_argument("--tprior2", type=int,default=15,dest='tprior2',help="split time prior")
	parser.add_argument("--mprior1", type=int,default=0,dest='mprior1',help="migration prior")
	parser.add_argument("--mprior2", type=int,default=20,dest='mprior2',help="migration prior")
	parser.add_argument("--kernel",default="uniform",dest="kernel",help="choose from gaussian or uniform kernel")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	parser.add_argument("--checkpoint", type=int,default=1e05,dest='checkpoint',help="checkpoint in yr")
	parser.add_argument("--migration_end", default=100,dest='migration_end',help="migration end in kyr")
	parser.add_argument("--interval_compare", type=int,default=13,dest='interval_compare',help="number of first intervals to compare")
	parser.add_argument("--last_iter", type=int,default=0,dest='last_iter',help="last_iter")
	parser.add_argument("--log", dest='log',help="logfile")
	args = parser.parse_args()

#	sys.stdout = open(args.log,'a',0)

	ITERATION=3

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'

#	TMRCA_bin=TMRCA_BIN()
#	TMRCA_bin.read_from_TMRCA_txt(args.file)

#	draw chi-square stats 2d plot:
#	psmc_result_matrix,macs_result_matrix=explore_stats_in_grid(TMRCA_bin,TMRCA_bin,TMRCA_bin)
#	draw_result_matrix(psmc_result_matrix,macs_result_matrix)

#	output TMRCA distribution from the formula
#	psmc_result_matrix,macs_result_matrix,TMRCA_dist,tm_list=explore_stats_in_grid_verbose(TMRCA_bin,TMRCA_bin,TMRCA_bin)
#	print_TMRCA_result(tm_list,TMRCA_dist)


	# test calc_f_density_complex function
	TMRCA_bin_pseudo=TMRCA_BIN()   # the pseudo diploid based on real data(the one you want to follow and put into macs)
	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,None,1)
	TMRCA_bin_pseudo.print_TMRCA_bin()
	TMRCA_bin_Af=TMRCA_BIN()
	TMRCA_bin_Af.read_from_psmc_mulitple(args.file1,None,1)
	TMRCA_bin_Af.print_TMRCA_bin()
	TMRCA_bin_Eu=TMRCA_BIN()
	TMRCA_bin_Eu.read_from_psmc_mulitple(args.file2,None,1)
	TMRCA_bin_Eu.print_TMRCA_bin()

	'''
	for i in range(len(TMRCA_bin_pseudo.time_bin)):
		print TMRCA_bin_pseudo.time_bin[i],TMRCA_bin_pseudo.pop_size[i]
		if TMRCA_bin_pseudo.time_bin[i]> 335300/30:
			TMRCA_bin_pseudo.pop_size[i]=7310*2.36/1.25
		elif TMRCA_bin_pseudo.time_bin[i]> 100000/30:
			TMRCA_bin_pseudo.pop_size[i]=14474*2.36/1.25
		print TMRCA_bin_pseudo.time_bin[i],TMRCA_bin_pseudo.pop_size[i]
	print TMRCA_bin_pseudo.pop_size
#	TMRCA_bin_pseudo.print_TMRCA_bin()
	'''
	TMRCA_bin_pseudo_sim=TMRCA_BIN()   # the pseudo diploid based on simulated data(the one you observe)
	TMRCA_bin_pseudo_sim.read_from_psmc_mulitple(args.psmc,None,1)
	TMRCA_bin_pseudo_sim.print_TMRCA_bin()
#	f_density_sim,P_sim=calc_f_density_integrate_single(10,15,TMRCA_bin_pseudo_sim,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density,P=calc_f_density_integrate_single_v2(10,15,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	draw_f_density_compare(f_density_sim,f_density,P_sim,P,TMRCA_bin_pseudo_sim,TMRCA_bin_pseudo)

#	f_density_sim,P_sim=calc_f_density_integrate_single(8,20,TMRCA_bin_pseudo_sim,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density,P=calc_f_density_integrate_single(8,20,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)

#	draw_f_density_compare(f_density_sim,f_density,P_sim,P,TMRCA_bin_pseudo_sim,TMRCA_bin_pseudo)
	
	TMRCA_bin_pseudo_sim.add_macs(args.file)
	psmc_result_matrix_scaled,macs_result_matrix_scaled,psmc_result_matrix,macs_result_matrix,TMRCA_dist,tm_list=explore_stats_in_grid_verbose(TMRCA_bin_pseudo_sim,TMRCA_bin_Af,TMRCA_bin_Eu)
	print_TMRCA_result(tm_list,TMRCA_dist)
#	draw_result_matrix(psmc_result_matrix,macs_result_matrix,'chi')
#	draw_result_matrix(psmc_result_matrix_scaled,macs_result_matrix_scaled,'chi_scaled')


