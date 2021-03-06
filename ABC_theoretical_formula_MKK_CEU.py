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
		print self.time_bin
		print self.genome_percentage
		print self.genome_percentage_macs

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

class particle:
	def __init__(self):
		self.t=0
		self.m=0
		self.NAf=0
		self.NEu=0

def sample_parameters(Ntrial,tprior,mprior1,mprior2):
	tprior_list=np.random.uniform(tprior[0],tprior[1],Ntrial)
	mprior_list1=np.random.uniform(mprior1[0],mprior1[1],Ntrial)
	mprior_list2=np.random.uniform(mprior2[0],mprior2[1],Ntrial)
	tm_list=zip(tprior_list,mprior_list1,mprior_list2)
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

# various population size after split; various ancestral population size
# hasn't got it right
'''
def calc_f_density_complex_v2(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint):  
	m1=m*1e-05/args.mu
	m2=m*1e-05/args.mu
	interval_Af = [0,4,6]
	NA=[TMRCA_bin_Af.pop_size[i] for i in interval_Af]
	interval_Eu=[0,4,8]
	NB=[TMRCA_bin_Eu.pop_size[i] for i in interval_Eu]
	t_change_recent=[TMRCA_bin_Eu.time_bin[i]*args.mu for i in [0,4,6]]
	Anc_interval = 0
	for i in range(len(TMRCA_bin_psuedo.time_bin)):
		if TMRCA_bin_psuedo.time_bin[i]*args.G>=max(checkpoint,t*1e04):
			Anc_interval=i
			break
	print 't_change_recent',t_change_recent
	print 'NA',NA,'NB',NB
	N_Anc=[TMRCA_bin_psuedo.pop_size[i] for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	t_change=[TMRCA_bin_psuedo.time_bin[i]*args.mu for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	theta_Anc=[4*i*args.mu for i in N_Anc]
	T=float(t*1e04)/args.G*args.mu
#	migration_end=float(migration_end*1e04)/args.G*args.mu
	print 'T',T,'migration end',migration_end
	t_bin=np.arange(0,0.0065,0.000001)
	f_density=np.zeros(len(t_bin))

	Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]]);
	M=linalg.expm(Q*T);
	ancestral_integral=np.zeros(len(N_Anc))
	recent_integral=[]
	for i in range(1,len(N_Anc)):
		if t_change[i-1]<T:
			continue
		if i==1:
			ancestral_integral[i]=-(2/theta_Anc[i-1])*(t_change[i-1]-T)
		else:
			ancestral_integral[i]=ancestral_integral[i-1]-(2/theta_Anc[i-2])*(t_change[i-1]-t_change[i-2])
	for i in range(len(t_change_recent)-1):
		if i==0:
			theta_Af=4*args.mu*NA[i]
	    	theta_Eu=4*args.mu*NB[i]
	    	Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]])
	    	m=linalg.expm(Q*t_change_recent[i+1])
	    	recent_integral.append(m)
	    else:
	    	theta_Af=4*args.mu*NA[i]
	    	theta_Eu=4*args.mu*NB[i]
	    	Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]])
			m=linalg.expm(Q*t_change_recent[i+1])
			recent_integral.append(recent_integral[-1]*m)
	for i in range(len(t_bin)):
	    if t_bin[i]<T:
	    	for j in range(1,len(t_change_recent)):
	    		if t_bin[i]<t_change_recent[j]:
	    			theta_Af=4*args.mu*NA[j]
	    			theta_Eu=4*args.mu*NB[j]
	    			Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]])
					m=linalg.expm(Q*t_bin[i])
					if j==1:
						m_final=m
					else:
						m_final=recent_integral[j-2]*m
					break
			if t_bin[i]>=t_change_recent[-1]:
				theta_Af=4*args.mu*NA[-1]
	    		theta_Eu=4*args.mu*NB[-1]
				Q=np.array([[-2*m1-2/theta_Af,2*m1,0,2/theta_Af,0],[m2,-m2-m1,m1,0,0],[0,2*m2,-2*m2-2/theta_Eu,0,2/theta_Eu],[0,0,0,-m1,m1],[0,0,0,m2,-m2]])
				m=linalg.expm(Q*t_bin[i])
				m_final=recent_integral[-1]*m
	        f_density[i]=m_final[1,0]*2/theta_Af+m_final[1,2]*2/theta_Eu;
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
#	print sum(f_density),sum(f_density)*0.000001
	return f_density
'''
# various population size after split; various ancestral population size
# with migration end time
def calc_f_density_complex_v3(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	interval_Af = [0,4,6]
	migration_done=False
	if migration_end==0:
		migration_done=True
	for i in range(len(interval_Af)):
		if TMRCA_bin_Eu.time_bin[interval_Af[i]]>=float(migration_end*1e04)/args.G:
			if migration_done is False:
				time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
				migration_done=True
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
		else:
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,label='after split'))
	if TMRCA_bin_psuedo.time_bin[checkpoint]*args.G>=t*1e04:
		Anc_interval=checkpoint
	else:
		for i in range(len(TMRCA_bin_psuedo.time_bin)):
			if TMRCA_bin_psuedo.time_bin[i]*args.G<=t*1e04 and TMRCA_bin_psuedo.time_bin[i+1]*args.G>=t*1e04:
				Anc_interval=i
				break
	theta_Anc=4*TMRCA_bin_psuedo.pop_size[Anc_interval]*args.mu
	time_bin_list.append(time_bin(time=float(t*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[-1]]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,theta_anc=theta_Anc,label='split time'))
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
#	for i in range(len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
	T=float(t*1e04)/args.G*args.mu
	print len(t_bin),T
	for i in range(len(t_bin)):
		if t_bin[i]<T:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time:
					m_final=time_bin_list[j].recent_integral*np.mat(linalg.expm(time_bin_list[j].Q*(t_bin[i]-time_bin_list[j].time)))
#					print i,t_bin[i],f_density[i]
					f_density[i]=m_final[1,0]*2./time_bin_list[j].theta1+m_final[1,2]*2./time_bin_list[j].theta2
					break
		else:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time and time_bin_list[j].label in ['split time','before split']:
#					if time_bin_list[j].theta_anc==0:
#						print t_bin[i],j,time_bin_list[j].theta_anc
					f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					print i,t_bin[i],f_density[i]
					break
#		print i,t_bin[i],f_density[i]
#	print sum(f_density),sum(f_density)*0.000001
	return f_density

# two periods of migration event, m1 is the recent one,small, m2 is the second, huge
def calc_f_density_complex_v4(t,m1,m2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_transit):
#	print t,m1,m2,migration_transit
	time_bin_list=[]
	interval_Af = [0,4,6]
	migration_done=False
	if migration_transit==0:
		migration_done=True
	for i in range(len(interval_Af)):
		if TMRCA_bin_Eu.time_bin[interval_Af[i]]>=float(migration_transit*1e04)/args.G:
			if migration_done is False:
#				if m1!=m2:
				time_bin_list.append(time_bin(time=float(migration_transit*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=math.pow(10,-m2)/args.mu,m2=math.pow(10,-m2)/args.mu,label='migration_end'))
				migration_done=True
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=math.pow(10,-m2)/args.mu,m2=math.pow(10,-m2)/args.mu,label='after split'))
		else:
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=math.pow(10,-m1)/args.mu,m2=math.pow(10,-m1)/args.mu,label='after split'))
#	print 'after split'
#	for i in range(len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
	if TMRCA_bin_psuedo.time_bin[checkpoint]*args.G>=t*1e04:
		Anc_interval=checkpoint
	else:
		for i in range(len(TMRCA_bin_psuedo.time_bin)):
			if TMRCA_bin_psuedo.time_bin[i]*args.G<=t*1e04 and TMRCA_bin_psuedo.time_bin[i+1]*args.G>=t*1e04:
				Anc_interval=i
				break
	theta_Anc=4*TMRCA_bin_psuedo.pop_size[Anc_interval]*args.mu
	time_bin_list.append(time_bin(time=float(t*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[-1]]*args.mu,theta_anc=theta_Anc,label='split time'))
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
#		time_bin_list[-1].print_time_bin()
	t_bin=np.arange(0,0.0065,0.000001)
	f_density=np.zeros(len(t_bin))
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
#	for i in range(len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
	T=float(t*1e04)/args.G*args.mu
#	print len(t_bin),T
	for i in range(len(t_bin)):
		if t_bin[i]<T:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time:
					m_final=time_bin_list[j].recent_integral*np.mat(linalg.expm(time_bin_list[j].Q*(t_bin[i]-time_bin_list[j].time)))
#					print i,t_bin[i],f_density[i]
					f_density[i]=m_final[1,0]*2./time_bin_list[j].theta1+m_final[1,2]*2./time_bin_list[j].theta2
					break
		else:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time and time_bin_list[j].label in ['split time','before split']:
#					if time_bin_list[j].theta_anc==0:
#						print t_bin[i],j,time_bin_list[j].theta_anc
					f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					print i,t_bin[i],f_density[i]
					break
#		print i,t_bin[i],f_density[i]
#	print sum(f_density),sum(f_density)*0.000001
	return f_density

def calc_f_density_complex(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint): # constant population size after split; various ancestral population size
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
#	print 'N and theta:',NA,NB,theta_Af,theta_Eu
	Anc_interval = 0
	for i in range(len(TMRCA_bin_psuedo.time_bin)):
		if TMRCA_bin_psuedo.time_bin[i]*args.G>=max(checkpoint,t*1e04):
			Anc_interval=i
			break
	N_Anc=[TMRCA_bin_psuedo.pop_size[i] for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	t_change=[TMRCA_bin_psuedo.time_bin[i]*args.mu for i in range(Anc_interval,len(TMRCA_bin_psuedo.time_bin))]
	theta_Anc=[4*i*args.mu for i in N_Anc]
#	print 'N_Anc:',N_Anc
#	print 'theta_Anc:',theta_Anc
#	print 't_change:',t_change
	T=float(t*1e04)/args.G*args.mu
#	print 'T',T
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
#	print sum(f_density),sum(f_density)*0.000001
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
	t_bin=np.arange(0,0.004,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi scaled')
	return stats

def calc_f_density_integrate_v2(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
#	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi')
	return (P,stats)

def calc_f_density_integrate_v3(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
#	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi')
	return (P,stats)

def calc_f_density_integrate_v4(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
#	f_density=calc_f_density_complex_v4(tm_list[index][0],tm_list[index][1],tm_list[index][2],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][3])
	f_density=calc_f_density_complex_v4(tm_list[index][0],tm_list[index][1],tm_list[index][2],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_transit)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi')
	return (P,stats)

def calc_f_density_integrate_single(t,m,migration_end,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density=calc_f_density_complex_v3(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	return f_density,P

def calc_f_density_integrate_single_v2(t,m1,m2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,migration_end):
	f_density=calc_f_density_complex_v4(t,m1,m2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	return f_density,P

def draw_f_density_compare(t,m,migration_end,f_density,P,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density,'r')
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax2 = fig.add_subplot(312,title='TMRCA distribution')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P,'r',label='simulation')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'b',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)
	ax3 = fig.add_subplot(313,title='population size(in 10k)')
	ax3.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r')
	ax3.grid(True)
	ax3.set_ylim([0,10])
	fig.savefig("sim_YRI_CEU_%s_%s_%s_f_density_v2.pdf" %(t,m,migration_end),format='pdf')

def draw_f_density_compare_v2(t,m,migration_end,f_density1,P1,f_density2,P2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density1,'r')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density2,'g')
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax2 = fig.add_subplot(312,title='TMRCA distribution')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P1,'r',label='simulation1')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P2,'g',label='simulation2')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'b',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)
	ax3 = fig.add_subplot(313,title='population size(in 10k)')
	ax3.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r')
	ax3.grid(True)
	ax3.set_ylim([0,10])
	fig.savefig("sim_MKK_CEU_%s_%s_%s_f_density_compare.pdf" %(t,m,migration_end),format='pdf')

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

def calc_f_density_integrate_wrapper_v4(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate_v4(*args2)
	return stats

def explore_stats_in_grid_verbose(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for migration_end in range(3):
		for t in range(args.tprior1,args.tprior2+1):
			for m in range(args.mprior1,args.mprior2+1,5):
				tm_list.append([t,m,migration_end])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v2,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	TMRCA_dist=[i[0] for i in result1]
	psmc_result=[i[1] for i in result1]
	psmc_result_matrix=psmc_result
	psmc_result_matrix=np.reshape(np.array(psmc_result),(3,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
	print 'psmc chi',psmc_result_matrix
	return psmc_result_matrix,TMRCA_dist,tm_list

def explore_stats_in_grid_verbose_v2(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for t in range(args.tprior1,args.tprior2+1):
		for m1 in range(2,5):
			for m2 in range(2,5):
				for migration_end in range(1,6):
					tm_list.append([t,m1,m2,migration_end])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v4,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	TMRCA_dist=[i[0] for i in result1]
	psmc_result=[i[1] for i in result1]
#	psmc_result_matrix=np.reshape(np.array(psmc_result),(3,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
#	print 'psmc chi',psmc_result_matrix
	return psmc_result,TMRCA_dist,tm_list

def draw_f_density(f_density):
	t_bin=np.arange(0,0.004,0.000001)
	fig = plt.figure()
	plt.semilogx(t_bin,f_density)
	plt.title('TMRCA density')
	plt.grid(True)
	fig.savefig("f_density",format='pdf')

def ABC_SMC_iter(iter,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,error_value,prev_weight,prev_tm_list):
	tm_list_new=[]
	done = 0
	tm_list_total=[]
	result_total=[]
	result_new=[]
	TMRCA_total =[]
	while done<args.Nparticle:
		if iter==0:
			tm_list=sample_parameters(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior11,args.mprior12),(args.mprior21,args.mprior22))
		else:
			tm_list=sample_parameters_weight(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior11,args.mprior12),(args.mprior21,args.mprior22),prev_tm_list,prev_weight)
		pool = multiprocessing.Pool(processes=30)
		result1=pool.map(calc_f_density_integrate_wrapper_v4,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(args.Nparticle)))
		TMRCA_dist=[i[0] for i in result1]
		result=[i[1] for i in result1]
		pool.close()
		pool.join()
		tm_list_total+=tm_list
		result_total+=result
#		print 'result:',result
		if len(error_value)<=iter:
			e = np.percentile(result,50)
			print 'iter',iter,'error value',e
			error_value.append(e)
		else:
			e=error_value[iter]
		for i in range(len(result)):
			if result[i]<=e:
				tm_list_new.append(tm_list[i])
				result_new.append(result[i])
				TMRCA_total.append(TMRCA_dist[i])
				done+=1
		print iter,subprocess.check_output('date').strip(),done
	tm_list_new=tm_list_new[:args.Nparticle]
	result_new=result_new[:args.Nparticle]
	TMRCA_total=TMRCA_total[:args.Nparticle]
	new_weight=plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight,TMRCA_total)
	return error_value,new_weight,tm_list_new

def update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_t,sigma_m1,sigma_m2):

	new_weight = np.zeros(len(new_tm_list))
	for i in range(len(new_weight)):
		for j in range(len(prev_weight)):
			if abs(prev_tm_list[j][0]-new_tm_list[i][0])<=sigma_t and abs(prev_tm_list[j][1]-new_tm_list[i][1])<=sigma_m1 and abs(prev_tm_list[j][2]-new_tm_list[i][2])<=sigma_m2:
				new_weight[i]+=prev_weight[j]*1/(2*sigma_t*2*sigma_m1*2*sigma_m2)
	new_weight = [1/float(i) for i in new_weight]
	new_weight = [float(i)/sum(new_weight) for i in new_weight]
	return new_weight

def plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight,TMRCA_dist):
	t=[i[0] for i in tm_list_total]
	m1=[i[1] for i in tm_list_total]
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior11,args.mprior12,num=(args.mprior12-args.mprior11)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t,m1,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s sampled param' %(iter))
	cax=ax1.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax1.set_ylabel('migration rate')
	t_accepted=[i[0] for i in tm_list_new]
	m1_accepted=[i[1] for i in tm_list_new]
	ax2 = fig.add_subplot(312,title='iter %s accepted param' %(iter))
	heatmap, xedges, yedges = np.histogram2d(t_accepted, m1_accepted,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax2.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax2.set_ylabel('migration rate')
	ax3 = fig.add_subplot(313)
	ax3.set_xlim([args.tprior1,args.tprior2])
	ax3.set_ylim([args.mprior11,args.mprior12])
	ax3.set_xlabel('split time')
	ax3.plot(t,m1,'ro')
	ax3.plot(t_accepted,m1_accepted,'bo')
	fig.savefig("version%s/iter%s_sampled_accepted_param_t_m1.pdf" %(args.version,iter),format='pdf')

	t=[i[0] for i in tm_list_total]
	m2=[i[2] for i in tm_list_total]
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior21,args.mprior22,num=(args.mprior22-args.mprior21)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t,m2,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s sampled param' %(iter))
	cax=ax1.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax1.set_ylabel('migration rate')
	t_accepted=[i[0] for i in tm_list_new]
	m2_accepted=[i[2] for i in tm_list_new]
	ax2 = fig.add_subplot(312,title='iter %s accepted param' %(iter))
	heatmap, xedges, yedges = np.histogram2d(t_accepted, m2_accepted,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax2.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax2.set_ylabel('migration rate')
	ax3 = fig.add_subplot(313)
	ax3.set_xlim([args.tprior1,args.tprior2])
	ax3.set_ylim([args.mprior21,args.mprior22])
	ax3.set_xlabel('split time')
	ax3.plot(t,m2,'ro')
	ax3.plot(t_accepted,m2_accepted,'bo')
	fig.savefig("version%s/iter%s_sampled_accepted_param_t_m2.pdf" %(args.version,iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s stats' %(iter))
	bins = 40
	ax1.hist(result_total,facecolor='green',alpha=0.5,histtype='stepfilled',label='genome percentage')
	ax1.legend(prop={'size':10})
	ax1.set_xlabel('chi square stats based on genome percentage')
	ax1.set_ylabel('count')
	ax1.axvline(x=e,linewidth=2,color='k')

	ax2 = fig.add_subplot(312)
	cax=ax2.scatter(t_accepted,m1_accepted,c=result_new,cmap=plt.cm.coolwarm)
	ax2.set_xlabel('chi square stats based on genome percentage')
	fig.colorbar(cax,orientation='vertical')

	ax3 = fig.add_subplot(313)
	cax=ax3.scatter(t_accepted,m2_accepted,c=result_new,cmap=plt.cm.coolwarm)
	ax3.set_xlabel('chi square stats based on genome percentage')
	fig.colorbar(cax,orientation='vertical')

	if iter==0:
		new_weight = [1./len(tm_list_new) for i in range(len(tm_list_new))]
	else:
		t_prev=[i[0] for i in prev_tm_list]
		m_prev1=[i[1] for i in prev_tm_list]
		m_prev2=[i[2] for i in prev_tm_list]
		if args.kernel == "uniform":
			sigma_t=0.5*(max(t_prev)-min(t_prev))
			sigma_m1=0.5*(max(m_prev1)-min(m_prev1))
			sigma_m2=0.5*(max(m_prev2)-min(m_prev2))
			print 'sigma_t', sigma_t
			print 'sigma_m1', sigma_m1
			print 'sigma_m2', sigma_m2
			new_weight=update_weight(iter,prev_tm_list,tm_list_new,prev_weight,sigma_t,sigma_m1,sigma_m2)
#	ax3 = fig.add_subplot(313)
#	cax=ax3.scatter(t_accepted,m_accepted,c=new_weight,cmap=plt.cm.coolwarm)
#	ax3.set_xlabel('weight')
#	fig.colorbar(cax,orientation='vertical')	
	fig.savefig("version%s/iter%s_distance_stats_scatterplot.pdf" %(args.version,iter),format='pdf')

	fout = open('version%s/iter%s_stats_weight.txt' %(args.version,iter),'w')
	for i in range(len(new_weight)):
		print >>fout,'%i\t%s\t%f\t%f\t%s' %(i,"\t".join(map(str,tm_list_new[i])),result_new[i],new_weight[i],"\t".join(map(str,TMRCA_dist[i])))
	fout.close()
	return new_weight

def get_id_from_sample(weight,random_number):
	weight_cumsum=np.cumsum(weight)
	random_id = 0
	assert random_number >=0 and random_number<1
	for i in range(len(weight_cumsum)):
		if weight_cumsum[i]>random_number:
			random_id=i
			break
	return random_id

def sample_parameters_weight(Ntrial,tprior,mprior1,mprior2,tm_list,weight):
	print 'sum of weight',sum(weight)
	t_list=np.zeros(Ntrial)
	m_list1=np.zeros(Ntrial)
	m_list2=np.zeros(Ntrial)
	t=[i[0] for i in tm_list]
	m1=[i[1] for i in tm_list]
	m2=[i[2] for i in tm_list]
	index = 0
	if args.kernel == "uniform":
		sigma_t=0.5*(max(t)-min(t))
		sigma_m1=0.5*(max(m1)-min(m1))
		sigma_m2=0.5*(max(m2)-min(m2))
		print '# passed:',len(t),len(m)
		print 'sigma_t', sigma_t
		print 'sigma_m', sigma_m
	while index<Ntrial:
		random_number=np.random.uniform(0,1)
		random_id = get_id_from_sample(weight,random_number)
		assert weight[random_id]>0
		rand_a=np.random.uniform(0,1)
		rand_direction=np.random.randint(2)
		if rand_direction:
			new_t=tm_list[random_id][0]+rand_a*sigma_t
		else:	
			new_t=tm_list[random_id][0]-rand_a*sigma_t
		if new_t<tprior[0] or new_t >tprior[1]:
			continue
		rand_b=np.random.uniform(0,1)
		rand_direction=np.random.randint(2)
		if rand_direction:
			new_m1=tm_list[random_id][1]+rand_b*sigma_m1
		else:	
			new_m1=tm_list[random_id][1]-rand_b*sigma_m1
		if new_m1<mprior1[0] or new_m1 >mprior1[1]:
			continue
		rand_c=np.random.uniform(0,1)
		rand_direction=np.random.randint(2)
		if rand_direction:
			new_m2=tm_list[random_id][2]+rand_c*sigma_m2
		else:	
			new_m2=tm_list[random_id][2]-rand_c*sigma_m2
		if new_m2<mprior2[0] or new_m2 >mprior2[1]:
			continue
		t_list[index]=new_t
		m_list1[index]=new_m1
		m_list2[index]=new_m2
		index+=1
		if index%100==0:
			print index,'sample done'
	return zip(t_list,m_list1,m_list2)


def read_previous(iter):
	f=open('iter%s_stats_weight.txt' %(iter),'r')
	tm_list = []
	weight=[]
	for line in f:
		line=line.rstrip().split('\t')
		t=float(line[1])
		m=float(line[2])
		w=float(line[4])
		tm_list.append([t,m])
		weight.append(w)
	return weight,tm_list

def print_TMRCA_result(tm_list,TMRCA_dist):
	f=open('sim_macs_MKK_CEU_grid_2mig_TMRCA_new2.txt','w')
	for i in range(len(tm_list)):
		print >>f,"%s\t%s" %("\t".join(map(str,tm_list[i])),"\t".join(map(str,TMRCA_dist[i])))

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
	parser.add_argument("--mprior11", type=int,default=0,dest='mprior11',help="migration prior")
	parser.add_argument("--mprior12", type=int,default=20,dest='mprior12',help="migration prior")
	parser.add_argument("--mprior21", type=int,default=0,dest='mprior21',help="migration prior")
	parser.add_argument("--mprior22", type=int,default=20,dest='mprior22',help="migration prior")
	parser.add_argument("--kernel",default="uniform",dest="kernel",help="choose from gaussian or uniform kernel")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	parser.add_argument("--version", dest='version',help="version")
	parser.add_argument("--checkpoint", type=int,default=1e05,dest='checkpoint',help="checkpoint in yr")
	parser.add_argument("--migration_transit", type=float, default=0,dest='migration_transit',help="migration end in 10kyr")
	parser.add_argument("--interval_compare", type=int,default=13,dest='interval_compare',help="number of first intervals to compare")
	parser.add_argument("--last_iter", type=int,default=0,dest='last_iter',help="last_iter")
	parser.add_argument("--log", dest='log',help="logfile")
	args = parser.parse_args()

#	sys.stdout = open(args.log,'a',0)

	os.popen('mkdir -p version%s' %(args.version))
	ITERATION=3

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'

	TMRCA_bin_pseudo=TMRCA_BIN()   # the pseudo diploid based on real data(the one you observe)
	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,None,4)
	TMRCA_bin_pseudo.print_TMRCA_bin()
	TMRCA_bin_Af=TMRCA_BIN()
	TMRCA_bin_Af.read_from_psmc_mulitple(args.file1,None,1)
	TMRCA_bin_Af.print_TMRCA_bin()
	TMRCA_bin_Eu=TMRCA_BIN()
	TMRCA_bin_Eu.read_from_psmc_mulitple(args.file2,None,1)
	TMRCA_bin_Eu.print_TMRCA_bin()


	# Test function
#	f_density1,P1=calc_f_density_integrate_single(6,1000,0,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density2,P2=calc_f_density_integrate_single(6,100,0,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print P1
#	print P2
#	draw_f_density_compare(6,40,args.migration_end,f_density,P,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	draw_f_density_compare_v2(6,1000,100,f_density1,P1,f_density2,P2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)

#	f_density1,P1=calc_f_density_integrate_single_v2(6,2,2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,1)
#	f_density2,P2=calc_f_density_integrate_single_v2(6,3,3,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,1)
#	print P1
#	print P2
#	draw_f_density_compare_v2(6,2,2,f_density1,P1,f_density2,P2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	## PERFORM GRID SEARCH

#	psmc_result_matrix,TMRCA_dist,tm_list=explore_stats_in_grid_verbose_v2(TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print_TMRCA_result(tm_list,TMRCA_dist)



	## PERFORM ABC
	if args.initial_e is None:
		error_value = []
		for i in range(args.last_iter):
			error_value.append(0)
	else:
		error_value = args.initial_e
	prev_weight = [1 for i in range(args.Nparticle)]
	prev_tm_list = [particle() for i in range(args.Nparticle)]
	prev_weight = [float(i)/sum(prev_weight) for i in prev_weight]
	for iter in range(ITERATION):
		if iter<args.last_iter-1:
			continue
		print 'iter%s begin' %(iter)
		print subprocess.check_output('date').strip()
		if iter==args.last_iter-1:
			print 'read output of this iteration'
			new_weight,new_tm_list=read_previous(iter)
		elif iter==args.last_iter:
			if iter>=1:
				prev_weight=new_weight
				prev_tm_list=new_tm_list
			error_value,new_weight,new_tm_list=ABC_SMC_iter(iter,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,error_value,prev_weight,prev_tm_list)
		else:
			prev_weight=new_weight
			prev_tm_list=new_tm_list
			error_value,new_weight,new_tm_list=ABC_SMC_iter(iter,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,error_value,prev_weight,prev_tm_list)

