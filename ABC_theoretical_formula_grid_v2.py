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

# various population size after split; various ancestral population size
# with migration end time
def calc_f_density_complex_v3(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12]
	elif args.type=='msmc':
		interval_Af = [0,1,2]
	migration_done=False
	if migration_end==0:
		migration_done=True
	for i in range(len(interval_Af)):
		if TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.G>=t*1e04:
			break
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
	for i in range(1,len(time_bin_list)):
		if time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',i,t,m
			time_bin_list[i].print_time_bin()
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

# two periods of migration event, m1 is the recent one,small, m2 is the second, huge
def calc_f_density_complex_v4(t,m1,m2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_transit):
	print t,m1,m2,migration_transit
	time_bin_list=[]
	interval_Af = [0,4,6]
	migration_done=False
	if migration_transit==0:
		migration_done=True
	for i in range(len(interval_Af)):
		if TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.G>=t*1e04:
			break
		if TMRCA_bin_Eu.time_bin[interval_Af[i]]>=float(migration_transit*1e04)/args.G:
			if migration_done is False:
#				if m1!=m2:
				time_bin_list.append(time_bin(time=float(migration_transit*1e04)/args.G*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i-1]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i-1]]*args.mu,m1=math.pow(10,-m2)/args.mu,m2=math.pow(10,-m2)/args.mu,label='migration_end'))
				migration_done=True
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=math.pow(10,-m2)/args.mu,m2=math.pow(10,-m2)/args.mu,label='after split'))
		else:
			time_bin_list.append(time_bin(time=TMRCA_bin_Eu.time_bin[interval_Af[i]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[i]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[i]]*args.mu,m1=math.pow(10,-m1)/args.mu,m2=math.pow(10,-m1)/args.mu,label='after split'))

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
	for i in range(len(time_bin_list)):
		time_bin_list[i].print_time_bin()
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

def calc_f_density_complex_v5(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12,14]
	elif args.type=='msmc':
		interval_Af = [0,1,2,3,4,5,6]
	migration_done=False
	if migration_end==0:
		migration_done=True
#	print t,m,migration_end
	Af_index=0
	Eu_index=0
	current_t=0
	while True:
#		print Af_index,Eu_index,current_t
		if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]==TMRCA_bin_Eu.time_bin[interval_Af[Eu_index]]:
			if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=float(migration_end*1e04):
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=time_bin_list[-1].theta1,theta2=time_bin_list[-1].theta2,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
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
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=time_bin_list[-1].theta1,theta2=time_bin_list[-1].theta2,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',time_bin_list[i].print_time_bin()

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

def calc_f_density_complex_v6_simple(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12]
	elif args.type=='msmc':
		interval_Af = [0,1,2,3,4,5,6]
	migration_done=False
	if migration_end==0:
		migration_done=True
#	print t,m
	Af_index=0
	Eu_index=0
	current_t=0
	time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[0]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[0]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
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
#	f_density11=np.zeros(len(t_bin))
#	f_density22=np.zeros(len(t_bin))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,i,time_bin_list[i].print_time_bin()

	T=float(t*1e04)/args.G*args.mu
#	print len(t_bin),T
	for i in range(len(t_bin)):
		if t_bin[i]<T:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time:
					m_final=time_bin_list[j].recent_integral*np.mat(linalg.expm(time_bin_list[j].Q*(t_bin[i]-time_bin_list[j].time)))
#					print i,t_bin[i],f_density[i]
					f_density[i]=m_final[1,0]*2./time_bin_list[j].theta1+m_final[1,2]*2./time_bin_list[j].theta2
#					f_density11[i]=m_final[0,0]*2./time_bin_list[j].theta1+m_final[0,2]*2./time_bin_list[j].theta2
#					f_density22[i]=m_final[2,0]*2./time_bin_list[j].theta1+m_final[2,2]*2./time_bin_list[j].theta2
					break
		else:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time and time_bin_list[j].label in ['split time','before split']:
#					if time_bin_list[j].theta_anc==0:
#						print t_bin[i],j,time_bin_list[j].theta_anc
					f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					f_density11[i]=(M[0,0]+M[0,1]+M[0,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					f_density22[i]=(M[2,0]+M[2,1]+M[2,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))

#					print i,t_bin[i],f_density[i]
					break
#		print i,t_bin[i],f_density[i]
#	print sum(f_density),sum(f_density)*0.000001
	return f_density

def calc_f_density_complex_v7_simple(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12]
	elif args.type=='msmc':
		interval_Af = [0,1,2,3,4,5,6]
	migration_done=False
	growth_begin=4
	NAf=TMRCA_bin_Af.pop_size[0]
	NEu_growth=TMRCA_bin_Eu.pop_size[0]
	NEu=min(TMRCA_bin_Eu.pop_size)
	print t,m,migration_end,NAf,NEu_growth,NEu
	if migration_end==0:
		migration_done=True
	if migration_end==0:
		time_bin_list.append(time_bin(time=TMRCA_bin_psuedo.time_bin[0]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[0]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[0]*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
	else:
		time_bin_list.append(time_bin(time=TMRCA_bin_psuedo.time_bin[0]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[0]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[0]*args.mu,m1=0,m2=0,label='after split'))
		if growth_begin>float(migration_end):
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
			time_bin_list.append(time_bin(time=float(growth_begin*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
		elif growth_begin<float(migration_end):
			time_bin_list.append(time_bin(time=float(growth_begin*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=0,m2=0,label='after split'))
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
		else:
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
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
#	f_density11=np.zeros(len(t_bin))
#	f_density22=np.zeros(len(t_bin))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,i,time_bin_list[i].print_time_bin()

	T=float(t*1e04)/args.G*args.mu
#	print len(t_bin),T
	for i in range(len(t_bin)):
		if t_bin[i]<T:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time:
					m_final=time_bin_list[j].recent_integral*np.mat(linalg.expm(time_bin_list[j].Q*(t_bin[i]-time_bin_list[j].time)))
#					print i,t_bin[i],f_density[i]
					f_density[i]=m_final[1,0]*2./time_bin_list[j].theta1+m_final[1,2]*2./time_bin_list[j].theta2
#					f_density11[i]=m_final[0,0]*2./time_bin_list[j].theta1+m_final[0,2]*2./time_bin_list[j].theta2
#					f_density22[i]=m_final[2,0]*2./time_bin_list[j].theta1+m_final[2,2]*2./time_bin_list[j].theta2
					break
		else:
			for j in range(len(time_bin_list)):
				if t_bin[i]>=time_bin_list[j].time and t_bin[i]<=time_bin_list[j+1].time and time_bin_list[j].label in ['split time','before split']:
#					if time_bin_list[j].theta_anc==0:
#						print t_bin[i],j,time_bin_list[j].theta_anc
					f_density[i]=(M[1,0]+M[1,1]+M[1,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					f_density11[i]=(M[0,0]+M[0,1]+M[0,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))
#					f_density22[i]=(M[2,0]+M[2,1]+M[2,2])*2/time_bin_list[j].theta_anc*np.exp(time_bin_list[j].ancestral_integral-(2/time_bin_list[j].theta_anc)*(t_bin[i]-time_bin_list[j].time))

#					print i,t_bin[i],f_density[i]
					break
#		print i,t_bin[i],f_density[i]
#	print sum(f_density),sum(f_density)*0.000001
	return f_density
def calc_distance_between_two_dist(A,B,C,cov,method):
	stats=0
	if method=='chi':   # chi-square
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=pow((A[i]-B[i]),2)/A[i]
	if method=='chi covariance':   # chi-square
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=pow((A[i]-B[i]),2)/cov[i,i]
#			print i,pow((A[i]-B[i]),2)/cov[i,i]
	if method=='total variation dist':
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=abs(A[i]-B[i])
	if method=='euclidean':
		for i in range(len(A)):
			if A[i]==0:
				continue
			stats+=(A[i]-B[i])*(A[i]-B[i])
		stats=math.sqrt(stats)
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
	f_density=calc_f_density_complex_v5(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
#	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi covariance')
	return (P,stats1,stats2)

def calc_f_density_integrate_v3(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
#	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	f_density=calc_f_density_complex_v5(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	return (P,stats)

def calc_f_density_integrate_v4(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density=calc_f_density_complex_v4(tm_list[index][0],tm_list[index][1],tm_list[index][2],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][3])
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,'chi')
	return (P,stats)

def calc_f_density_integrate_v5(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
#	f_density=calc_f_density_complex_v3(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	f_density=calc_f_density_complex_v7_simple(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'dist')
	return (P,stats1,stats2)

def calc_f_density_integrate_single(t,m,migration_end,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density=calc_f_density_complex_v5(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
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
	fig.savefig("sim_MKK_CEU_%s_%s_%s_f_density.pdf" %(t,m,migration_end),format='pdf')

def draw_f_density_compare_v2(t,m,migration_end,f_density1,P1,f_density2,P2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density1,'r')
	ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density2,'g')
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax2 = fig.add_subplot(312,title='TMRCA distribution')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P1,'r',label='70k,20')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P2,'g',label='100k,20')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'b',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)
	ax3 = fig.add_subplot(313,title='population size(in 10k)')
	ax3.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r')
	ax3.grid(True)
	ax3.set_ylim([0,10])
	fig.savefig("sim_San_CEU_%s_%s_%s_f_density_compare.pdf" %(t,m,migration_end),format='pdf')

def draw_f_density_compare_v3(t,m,migration_end,f_density,P,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA denstiy')
	color=['r','g','b','y','c']
	label=['70k,20','80k,20','90k,20','100k,20']
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
	fig.savefig("sim_San_CEU_%s_%s_%s_f_density_compare.pdf" %(t,m,migration_end),format='pdf')

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

def calc_f_density_integrate_wrapper_v5(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate_v5(*args2)
	return stats

def explore_stats_in_grid_verbose(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for migration_end in range(5):
		for t in range(args.tprior1,args.tprior2+1):
			for m in range(args.mprior1,args.mprior2+1,5):
				tm_list.append([t,m,migration_end])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	psmc_result1=[]
#	for i in [len(tm_list)-1]:
#		a,b=calc_f_density_integrate_v3(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,i)
#		psmc_result1.append(b)
	print 'Use calc_f_density_complex_v7_simple'
	result1=pool.map(calc_f_density_integrate_wrapper_v5,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	TMRCA_dist=[i[0] for i in result1]
	psmc_result1=[i[1] for i in result1]
	psmc_result2=[i[2] for i in result1]
	psmc_result_matrix=psmc_result1
	print len(psmc_result_matrix),1*(args.tprior2-args.tprior1+1)*((args.mprior2-args.mprior1)/5+1)
	psmc_result_matrix1=np.reshape(np.array(psmc_result1),(5,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
	psmc_result_matrix2=np.reshape(np.array(psmc_result2),(5,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
	print 'psmc chi',psmc_result_matrix1
	print 'psmc dist',psmc_result_matrix2
	return result1,tm_list

def explore_stats_in_grid_verbose_v2(TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	tm_list=[]
	for t in range(args.tprior1,args.tprior2+1):
		for m1 in range(2,6):
			for m2 in range(2,6):
				for migration_end in range(1,4):
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
			tm_list=sample_parameters(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2))
		else:
			tm_list=sample_parameters_weight(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2),prev_tm_list,prev_weight)
		pool = multiprocessing.Pool(processes=30)
		result1=pool.map(calc_f_density_integrate_wrapper_v3,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(args.Nparticle)))
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

def update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_t,sigma_m):

	new_weight = np.zeros(len(new_tm_list))
	for i in range(len(new_weight)):
		for j in range(len(prev_weight)):
			if abs(prev_tm_list[j][0]-new_tm_list[i][0])<=sigma_t and abs(prev_tm_list[j][1]-new_tm_list[i][1])<=sigma_m:
				new_weight[i]+=prev_weight[j]*1/(2*sigma_t*2*sigma_m)
	new_weight = [1/float(i) for i in new_weight]
	new_weight = [float(i)/sum(new_weight) for i in new_weight]
	return new_weight

def plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight,TMRCA_dist):
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
	fig.savefig("version%s/iter%s_sampled_accepted_param.pdf" %(args.version,iter),format='pdf')

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
	fig.savefig("version%s/iter%s_distance_stats_scatterplot.pdf" %(args.version,iter),format='pdf')

	fout = open('version%s/iter%s_stats_weight.txt' %(args.version,iter),'w')
	for i in range(len(new_weight)):
		t=tm_list_new[i][0]
		m=tm_list_new[i][1]
		print >>fout,'%i\t%s\t%s\t%f\t%f\t%s' %(i,t,m,result_new[i],new_weight[i],"\t".join(map(str,TMRCA_dist[i])))
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

def sample_parameters_weight(Ntrial,tprior,mprior,tm_list,weight):
	print 'sum of weight',sum(weight)
	t_list=np.zeros(Ntrial)
	m_list=np.zeros(Ntrial)
	t=[i[0] for i in tm_list]
	m=[i[1] for i in tm_list]
	index = 0
	if args.kernel == "uniform":
		sigma_t=0.5*(max(t)-min(t))
		sigma_m=0.5*(max(m)-min(m))
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
			new_m=tm_list[random_id][1]+rand_b*sigma_m
		else:	
			new_m=tm_list[random_id][1]-rand_b*sigma_m
		if new_m<mprior[0] or new_m >mprior[1]:
			continue
		t_list[index]=new_t
		m_list[index]=new_m
		index+=1
		if index%100==0:
			print index,'sample done'
	return zip(t_list,m_list)


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
	f=open('sim_macs_GIH_CEU_grid_TMRCA_use_v5.txt','w')
	for i in range(len(tm_list)):
		print >>f,"%s\t%s" %("\t".join(map(str,tm_list[i])),"\t".join(map(str,TMRCA_dist[i])))

def print_TMRCA_result_v2(tm_list,result):
	f=open('%s_grid_TMRCA_psmc.txt' %(args.file3),'w')
	for i in range(len(tm_list)):
		print >>f,"%s\t%s\t%s" %("\t".join(map(str,tm_list[i])),"\t".join(map(str,result[i][0])),"\t".join(map(str,result[i][1:])))

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--psmc", dest='psmc',help="psmc file prefix")  
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--Ntrial", default=1500,type=int,dest='Ntrial',help="Number of trials per iteration")
	parser.add_argument("--Nparticle", default=1000,type=int,dest='Nparticle',help="number of particles to reach")
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

	os.popen('mkdir -p version%s' %(args.version))
	ITERATION=3

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'

	TMRCA_bin_pseudo=TMRCA_BIN()   # the pseudo diploid based on real data(the one you observe)
	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,['0','1','2','3'],4)
#	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,['0','1','2','3'],4)
#	TMRCA_bin_pseudo.calc_covariance()
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
#	f_density1,P1=calc_f_density_integrate_single(7,50,4,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density2,P2=calc_f_density_integrate_single(8,50,4,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density3,P3=calc_f_density_integrate_single(9,50,4,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	f_density4,P4=calc_f_density_integrate_single(10,50,4,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print P1
#	print P2
#	draw_f_density_compare(6.24,987,args.migration_end,f_density1,P1,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	draw_f_density_compare_v3(7,50,4,[f_density1,f_density2,f_density3,f_density4],[P1,P2,P3,P4],TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)

#	f_density1,P1=calc_f_density_integrate_single_v2(6,2,2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,1)
#	f_density2,P2=calc_f_density_integrate_single_v2(6,3,3,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,1)
#	print P1
#	print P2
#	draw_f_density_compare_v2(6,2,2,f_density1,P1,f_density2,P2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	## PERFORM GRID SEARCH

#	f_density1,P1=calc_f_density_integrate_single(8,100,0,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	psmc_result_matrix,tm_list=explore_stats_in_grid_verbose(TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	print_TMRCA_result_v2(tm_list,psmc_result_matrix)
#	calc_f_density_integrate_v2([[4,100,0]],TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu,0)

	'''
	## PERFORM ABC
	if args.initial_e is None:
		error_value = []
		for i in range(args.last_iter):
			error_value.append(0)
	else:
		error_value = args.initial_e
	prev_weight = [1 for i in range(args.Nparticle)]
	prev_tm_list = [[0,0] for i in range(args.Nparticle)]
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
	'''
