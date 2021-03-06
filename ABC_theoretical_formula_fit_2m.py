#!/usr/bin/env python
# python ABC_calc_TMRCA_v2.py
# Shiya Song
# 17th April 2015

import argparse,time,itertools,multiprocessing,math,glob,os,subprocess,pickle,sys
import numpy as np
from scipy import linalg
from scipy.optimize import fsolve
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


class TMRCA_BIN:
	def __init__(self,label='',time_bin=[],genome_percentage=[],genome_percentage_formula=[],pop_size=[],pop_size_formula=[]):
		self.time_bin=time_bin   # in generations
		self.genome_percentage=genome_percentage
		self.genome_percentage_macs=[]
		self.genome_percentage_formula=genome_percentage_formula
		self.pop_size=pop_size
		self.pop_size_formula=pop_size_formula
		self.rho=0
		self.theta=0
		self.covariance=0
		self.label=label

	def calc_TMRCA_dist_formula_v2(self):   # given pop_size, calculate the expected TMRCA, should match the TMRCA from a diploid genome
		TMRCA=np.zeros(len(self.time_bin))
		TMRCA_v2=np.zeros(len(self.time_bin))
		summation=np.zeros(len(self.time_bin))
		for i in range(1,len(self.time_bin)):
			theta=4*self.pop_size[i-1]*args.mu
			summation[i]=-2/theta*args.mu*(self.time_bin[i]-self.time_bin[i-1])+summation[i-1]
		for j in range(len(self.time_bin)-1):
			theta=4*self.pop_size[j]*args.mu
			TMRCA_v2[j]=math.exp(summation[j])*(1-math.exp(-2/theta*args.mu*(self.time_bin[j+1]-self.time_bin[j])))
		TMRCA_v2[-1]=math.exp(summation[-1])*(1-math.exp(-2/theta*(0.0065-self.time_bin[-1]*args.mu)))
		# numerical integration as below
		'''
		t_bin=np.arange(0,0.0065,0.000001)
		f_density=np.zeros(len(t_bin))
		for i in range(len(t_bin)):
			for j in range(len(self.time_bin)-1):
				if t_bin[i]>=self.time_bin[j]*args.mu and t_bin[i]<=self.time_bin[j+1]*args.mu:
					theta=4*self.pop_size[j]*args.mu
					f_density[i]=2/theta*np.exp(-2/theta*(t_bin[i]-self.time_bin[j]*args.mu)+summation[j])
					break
			if t_bin[i]>=self.time_bin[-1]*args.mu:
				theta=4*self.pop_size[-1]*args.mu
				f_density[i]=2/theta*np.exp(-2/theta*args.mu*(t_bin[i]-self.time_bin[j]*args.mu)+summation[-1])
		print np.sum(f_density),sum(f_density)*0.000001
		for i in range(len(t_bin)):
			for j in range(len(self.time_bin)-1):
				if t_bin[i]>=self.time_bin[j]*args.mu and t_bin[i]<=self.time_bin[j+1]*args.mu:
					TMRCA[j]=TMRCA[j]+f_density[i]/np.sum(f_density)
					break
			if t_bin[i]>=self.time_bin[-1]*args.mu:
				TMRCA[-1]=TMRCA[-1]+f_density[i]/np.sum(f_density)
		'''
		print TMRCA_v2,np.sum(TMRCA_v2)
		print self.genome_percentage,sum(self.genome_percentage)
		self.genome_percentage_formula=TMRCA_v2

	def calc_theta_t_from_TMRCA_v2(self):
		index = 0
		lambda_array=[]  # lambda_t=1/2/N_t
		prev_left=0
		print self.genome_percentage
		zero=False
		for i in range(1,len(self.time_bin)):
			if self.genome_percentage[i-1]==0:
				lambda_array.append(0)
				zero=True
				continue
			if zero:
				lambda_t=0-math.log(1-self.genome_percentage[i-1]/math.exp(prev_left))/(self.time_bin[i]-self.time_bin[0])
				prev_left+=0-lambda_t*(self.time_bin[i]-self.time_bin[0])
				zero=False
			else:
				try:
					lambda_t=0-math.log(1-self.genome_percentage[i-1]/math.exp(prev_left))/(self.time_bin[i]-self.time_bin[i-1])
				except ValueError:
					print 'Wrong',i,self.genome_percentage[i-1],math.exp(prev_left),self.time_bin[i]-self.time_bin[i-1]
				prev_left+=0-lambda_t*(self.time_bin[i]-self.time_bin[i-1])
			lambda_array.append(lambda_t)
#		lambda_t=0-math.log(1-self.genome_percentage[-1]/math.exp(prev_left))/(0.0065/args.mu-self.time_bin[-1])
		lambda_array.append(lambda_t)
		print lambda_array
		self.pop_size_formula=np.zeros(len(lambda_array))
		a=range(len(lambda_array))
		a.reverse()
		for i in a:
			if lambda_array[i]!=0:
				self.pop_size_formula[i]=1./2/lambda_array[i]
			else:
				self.pop_size_formula[i]=self.pop_size_formula[i+1]
		print self.pop_size_formula
		print len(self.pop_size_formula)
#		print self.pop_size

	def calc_theta_t_from_TMRCA(self):
		prev_size=self.pop_size[0]
		index = 0
		lambda_array=[]  # lambda_t=1/2/N_t
		prev_left=0
		for i in range(1,len(self.pop_size)):
			if self.pop_size[i]==prev_size:
				continue
			else:
				if index==0:
					lambda_t=0-math.log(1-sum(self.genome_percentage[index:i]))/(self.time_bin[i]-self.time_bin[index])
					for j in range(index,i):
						lambda_array.append(lambda_t)
					print i,index,lambda_array
					index=i
					prev_size=self.pop_size[i]
					prev_left+=0-lambda_t*(self.time_bin[i]-self.time_bin[index])
				else:
					c=sum(self.genome_percentage[index:i])/math.exp(prev_left)
#					lambda_t=fsolve(func,x0=lambda_array[-1],args=(self.time_bin[index],self.time_bin[i],c))
					lambda_t=0-math.log(1-c)/(self.time_bin[i]-self.time_bin[index])
					for j in range(index,i):
						lambda_array.append(lambda_t)
					print i,index,lambda_array
					index=i
					prev_size=self.pop_size[i]
					prev_left+=0-lambda_t*(self.time_bin[i]-self.time_bin[index])
		c=sum(self.genome_percentage[index:])/math.exp(prev_left)
#		lambda_t=fsolve(func,x0=lambda_array[-1],args=(self.time_bin[index],self.time_bin[-1],c))
		lambda_t=0-math.log(1-c)/(self.time_bin[i]-self.time_bin[index])
		for j in range(index,len(self.pop_size)):
			lambda_array.append(lambda_t)
		print i,index,lambda_array
		self.pop_size_formula=[1./2/i for i in lambda_array]
		print self.pop_size_formula
		print self.pop_size

	def draw_TMRCA(self):
		t_bin=np.arange(0,0.0065,0.000001)
		fig = plt.figure()
		ax1 = fig.add_subplot(211,title='TMRCA distribution')
		ax1.semilogx([i*30 for i in self.time_bin],self.genome_percentage_formula,'r',label='formula')
#		ax1.semilogx([i*30 for i in self.time_bin],self.genome_percentage_formula_v2,'b',label='formula integration')
		ax1.semilogx([i*30 for i in self.time_bin],self.genome_percentage,'k',label='observed')
		ax1.legend(prop={'size':10})
		ax1.grid(True)
		ax1.set_xlim([1e03,1e08])
		ax2 = fig.add_subplot(212,title='population size')
		ax2.semilogx([i*30 for i in self.time_bin],[float(i)/1e04 for i in self.pop_size_formula],'r',label='formula')
		ax2.semilogx([i*30 for i in self.time_bin],[float(i)/1e04 for i in self.pop_size],'k',label='observed')
		ax2.grid(True)
		ax2.set_ylim([0,10])
		fig.savefig("%s_pop_size_vs_TMRCA_dist.pdf" %(self.label),format='pdf')

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

	def read_from_msmc(self,file,verbose=False):
		f = open(file,'r')
		line = f.readline()
		for line in f:
			line=line.rstrip().split('\t')
			self.time_bin.append(abs(float(line[1]))/args.G)
			self.pop_size.append(abs(float(line[3]))*10000)
			self.genome_percentage.append(0)
		if verbose==True:
			print self.time_bin
			print self.pop_size

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

def calc_f_density_complex_v6(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
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
	time_bin_list[0].print_time_bin()
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

def calc_f_density_complex_v6_simple(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
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
	time_bin_list[0].print_time_bin()
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

def calc_f_density_complex_v6_explore(t,m,NAf,NEu,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index   # needs to fix code when migration_end!=0
	time_bin_list=[]
	migration_done=False
#	print t,m,NAf,NEu,migration_end
	time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
	if migration_end!=0:
		time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=0,m2=0,label='migration_end'))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,NAf,NEu,migration_end
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

def calc_f_density_complex_v6_explore_YRI_CEU(t,m,NAf,NEu,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	migration_done=False
#	print t,m,NAf,NEu,migration_end
	NEu_growth=5e04
	growth_begin=3e04
	if migration_end==0:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
	else:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=0,m2=0,label='after split'))
		if growth_begin>float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
		elif growth_begin<float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,NAf,NEu,migration_end
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

def calc_f_density_complex_v6_explore_GIH_CEU(t,m,NAf,NEu,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	migration_done=False
#	print t,m,NAf,NEu,migration_end
	NEu_growth=5e04
	growth_begin=3e04
	if migration_end==0:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NEu_growth*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
	else:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NEu_growth*args.mu,theta2=4*NEu_growth*args.mu,m1=0,m2=0,label='after split'))
		if growth_begin>float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
		elif growth_begin<float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,NAf,NEu,migration_end
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

def calc_f_density_complex_v6_explore_YRI_MKK(t,m,NAf,NEu,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	migration_done=False
#	print t,m,NAf,NEu,migration_end
	NEu_growth=5e04
	growth_begin=3e04
	if migration_end==0:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
	else:
		time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=0,m2=0,label='after split'))
		if growth_begin>float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu_growth*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='migration_end'))
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
		elif growth_begin<float(migration_end*1e04):
			time_bin_list.append(time_bin(time=float(growth_begin)/args.G*args.mu,theta1=4*NAf*args.mu,theta2=4*NEu*args.mu,m1=m*1e-05/args.mu,m2=m*1e-05/args.mu,label='after split'))
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
#	time_bin_list[0].print_time_bin()
	for i in range(1,len(time_bin_list)):
#		time_bin_list[i].print_time_bin()
		if  time_bin_list[i].time < time_bin_list[i-1].time:
			print 'wrong',t,m,NAf,NEu,migration_end
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

def calc_f_density_complex_v7(t,m1,m2,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,checkpoint,migration_end):  # checkpoint is an index
	time_bin_list=[]
	if args.type=='psmc':
		interval_Af = [0,4,6,8,10,12,14,16]
	elif args.type=='msmc':
		interval_Af = [0,1,2,3,4,5,6]
	migration_done=False
	if migration_end==0:
		migration_done=True
#	print t,m1,m2
	Af_index=0
	Eu_index=0
	current_t=0
#	time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[0]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[0]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[0]*args.mu,m1=m1*1e-05/args.mu,m2=m2*1e-05/args.mu,label='after split'))

	while True:
		print Af_index,Eu_index,current_t
		if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]==TMRCA_bin_Eu.time_bin[interval_Af[Eu_index]]:
			if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=float(migration_end*1e04):
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=time_bin_list[-1].theta1,theta2=time_bin_list[-1].theta2,m1=m1*1e-05/args.mu,m2=m2*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,m1=m1*1e-05/args.mu,m2=m2*1e-05/args.mu,label='after split'))
			else:
				if TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=TMRCA_bin_Af.time_bin[interval_Af[Af_index]]*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,label='after split'))
		else:
			if current_t*args.G>=float(migration_end*1e04):
				if migration_done is False:
					time_bin_list.append(time_bin(time=float(migration_end*1e04)/args.G*args.mu,theta1=time_bin_list[-1].theta1,theta2=time_bin_list[-1].theta2,m1=m1*1e-05/args.mu,m2=m2*1e-05/args.mu,label='migration_end'))
					migration_done=True
				if current_t*args.G>=t*1e04:
					break
				time_bin_list.append(time_bin(time=current_t*args.mu,theta1=4*TMRCA_bin_Af.pop_size[interval_Af[Af_index]]*args.mu,theta2=4*TMRCA_bin_Eu.pop_size[interval_Af[Eu_index]]*args.mu,m1=m1*1e-05/args.mu,m2=m2*1e-05/args.mu,label='after split'))
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
			print 'wrong',t,m1,m2,i,time_bin_list[i].print_time_bin()

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

def calc_f_density_integrate_single(t,m,migration_end,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density,f_density11,f_density22=calc_f_density_complex_v5(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P11=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P22=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density)
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Af.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Af.time_bin[j+1]*args.mu:
				P11[j]=P11[j]+f_density11[i]/sum(f_density11);
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Eu.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Eu.time_bin[j+1]*args.mu:
				P22[j]=P22[j]+f_density22[i]/sum(f_density22);
				break
	return [f_density,f_density11,f_density22],[P,P11,P22]

def calc_relative_cross_coalescent_curve(t,m,migration_end,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	f_density,f_density11,f_density22=calc_f_density_complex_v6(t,m,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P11=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P22=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density)
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Af.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Af.time_bin[j+1]*args.mu:
				P11[j]=P11[j]+f_density11[i]/sum(f_density11);
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Eu.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Eu.time_bin[j+1]*args.mu:
				P22[j]=P22[j]+f_density22[i]/sum(f_density22);
				break

	TMRCA_bin_11=TMRCA_BIN(label='YRI',genome_percentage=list(P11),time_bin=TMRCA_bin_Af.time_bin)
	TMRCA_bin_12=TMRCA_BIN(label='YRI-CEU',genome_percentage=list(P),time_bin=TMRCA_bin_psuedo.time_bin)
	TMRCA_bin_22=TMRCA_BIN(label='CEU',genome_percentage=list(P22),time_bin=TMRCA_bin_Eu.time_bin)
	TMRCA_bin_11.calc_theta_t_from_TMRCA_v2()
	TMRCA_bin_12.calc_theta_t_from_TMRCA_v2()
	TMRCA_bin_22.calc_theta_t_from_TMRCA_v2()
	return f_density,f_density11,f_density22,TMRCA_bin_11,TMRCA_bin_12,TMRCA_bin_22

def draw_TMRCA_and_coalescence_rate(t,m,migration_end,TMRCA_bin_11,TMRCA_bin_12,TMRCA_bin_22):
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='TMRCA distribution')
	color=['r','g','b','y','c']
	label=['pop 1-2','pop 1','pop 2']
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax1.semilogx([i*30 for i in TMRCA_bin_12.time_bin],TMRCA_bin_12.genome_percentage,color[0],label=label[0])
	ax1.semilogx([i*30 for i in TMRCA_bin_11.time_bin],TMRCA_bin_11.genome_percentage,color[1],label=label[1])
	ax1.semilogx([i*30 for i in TMRCA_bin_22.time_bin],TMRCA_bin_22.genome_percentage,color[2],label=label[2])
	ax1.legend(prop={'size':10})
	ax2 = fig.add_subplot(312,title='population size')
	ax2.semilogx([i*30 for i in TMRCA_bin_12.time_bin],[float(i)/1e04 for i in TMRCA_bin_12.pop_size_formula],color[0],label=label[0])
	ax2.semilogx([i*30 for i in TMRCA_bin_11.time_bin],[float(i)/1e04 for i in TMRCA_bin_11.pop_size_formula],color[1],label=label[1])
	ax2.semilogx([i*30 for i in TMRCA_bin_22.time_bin],[float(i)/1e04 for i in TMRCA_bin_22.pop_size_formula],color[2],label=label[2])
	ax2.grid(True)
	ax2.set_xlim([1e03,1e08])
	ax2.set_ylim([0,10])
	relative_cross_coalescent_rate=np.zeros(len(TMRCA_bin_11.time_bin))
	for i in range(len(TMRCA_bin_11.time_bin)):
		relative_cross_coalescent_rate[i]=2./TMRCA_bin_12.pop_size_formula[i]/(1./TMRCA_bin_11.pop_size_formula[i]+1./TMRCA_bin_22.pop_size_formula[i])
	ax3 = fig.add_subplot(313,title='relative cross coalescence rate')
	ax3.semilogx([i*30 for i in TMRCA_bin_11.time_bin],relative_cross_coalescent_rate)
	ax3.set_xlim([1e03,1e08])
	ax3.grid(True)
	fig.savefig("sim_%s_%s_%s_%s_%s_relative_cross_coalecence_rate_v2.pdf" %(args.pop,t,m,migration_end,args.type),format='pdf')

def draw_f_density_compare_v4(t,m,migration_end,f_density,P,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(511,title='TMRCA denstiy')
	color=['r','g','b','y','c']
#	label=['70k,20','80k,20','90k,20','100k,20']
	label=['pop 1-2','pop 1','pop 2']
	for index in range(len(f_density)):
		ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density[index],color[index])
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])

	ax2 = fig.add_subplot(512,title='TMRCA distribution')
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P[0],color[0],label=label[0])
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'k',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)

	ax3 = fig.add_subplot(513,title='TMRCA distribution')
	ax3.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],P[1],color[1],label=label[1])
	ax3.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],TMRCA_bin_Af.genome_percentage,'k',label='observed')
	ax3.legend(prop={'size':10})
	ax3.grid(True)
	ax4 = fig.add_subplot(514,title='TMRCA distribution')
	ax4.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],P[2],color[2],label=label[2])
	ax4.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],TMRCA_bin_Eu.genome_percentage,'k',label='observed')
	ax4.legend(prop={'size':10})
	ax4.grid(True)
	ax5 = fig.add_subplot(515,title='population size(in 10k)')
	ax5.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r',label='pop 1-2')
	ax5.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],[float(i)/1e04 for i in TMRCA_bin_Af.pop_size],'g',label='pop 1')
	ax5.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],[float(i)/1e04 for i in TMRCA_bin_Eu.pop_size],'b',label='pop 2')
	ax5.legend(prop={'size':10})
	ax5.grid(True)
	ax5.set_ylim([0,10])

	fig.savefig("sim_%s_%s_%s_%s_%s_TMRCA.pdf" %(args.pop,t,m,migration_end,args.type),format='pdf')

def draw_f_density_compare_v3(t,m,migration_end,f_density,P,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu):
	t_bin=np.arange(0,0.0065,0.000001)
	fig = plt.figure()
	ax1 = fig.add_subplot(411,title='TMRCA denstiy')
	color=['r','g','b','y','c']
	label=['70k,20','80k,20','90k,20','100k,20']
#	label=['50k,50','60k,50','70k,50','80k,50','90k,50']
	for index in range(len(f_density)):
		ax1.semilogx([i/1.25e-08*30 for i in t_bin],f_density[index],color[index])
	ax1.grid(True)
	ax1.set_xlim([1e03,1e08])
	ax2 = fig.add_subplot(412,title='TMRCA distribution')
	for index in range(len(P)):
		ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],P[index],color[index],label=label[index])
	ax2.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],TMRCA_bin_psuedo.genome_percentage,'k',label='observed')
	ax2.legend(prop={'size':10})
	ax2.grid(True)
	ax3 = fig.add_subplot(413,title='population size(in 10k)')
	ax3.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.pop_size],'r',label='pop 1-2')
	ax3.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],[float(i)/1e04 for i in TMRCA_bin_Af.pop_size],'b',label='pop 1')
	ax3.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],[float(i)/1e04 for i in TMRCA_bin_Eu.pop_size],'g',label='pop 2')
	ax3.legend(prop={'size':10})
	ax3.grid(True)
	ax3.set_ylim([0,10])
	ax4 = fig.add_subplot(414,title='TMRCA distribution')
	ax4.semilogx([i*30 for i in TMRCA_bin_psuedo.time_bin],[float(i)/1e04 for i in TMRCA_bin_psuedo.genome_percentage],'r',label='pop 1-2')
	ax4.semilogx([i*30 for i in TMRCA_bin_Af.time_bin],[float(i)/1e04 for i in TMRCA_bin_Af.genome_percentage],'b',label='pop 1')
	ax4.semilogx([i*30 for i in TMRCA_bin_Eu.time_bin],[float(i)/1e04 for i in TMRCA_bin_Eu.genome_percentage],'g',label='pop 2')
	ax4.grid(True)
	fig.savefig("sim_%s_%s_%s_%s_%s_f_density_compare_v2.pdf" %(args.pop,t,m,migration_end,args.type),format='pdf')

def calc_f_density_integrate_v2(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):
	f_density,f_density11,f_density22=calc_f_density_complex_v6_simple(tm_list[index][0],tm_list[index][1],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,tm_list[index][2])
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P11=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P22=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density)
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Af.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Af.time_bin[j+1]*args.mu:
				P11[j]=P11[j]+f_density11[i]/sum(f_density11);
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Eu.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Eu.time_bin[j+1]*args.mu:
				P22[j]=P22[j]+f_density22[i]/sum(f_density22);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	stats2=calc_distance_between_two_dist(TMRCA_bin_Af.genome_percentage,P11,TMRCA_bin_Af.genome_percentage_macs,TMRCA_bin_Af.covariance,'chi')
	stats3=calc_distance_between_two_dist(TMRCA_bin_Eu.genome_percentage,P22,TMRCA_bin_Eu.genome_percentage_macs,TMRCA_bin_Eu.covariance,'chi')
#	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi covariance')
	return (P,stats1,stats2,stats3)

def calc_f_density_integrate_v3(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):  # fit 1 stats
	f_density=calc_f_density_complex_v7(tm_list[index][0],tm_list[index][1],tm_list[index][2],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density)
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	info=list(P)
	info1="\t".join(map(str,info))
	del info
	info2=stats1
	return (info1,info2)

def calc_f_density_integrate_v4(tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,index):  # fit 3 stats
	f_density,f_density11,f_density22=calc_f_density_complex_v6_explore(tm_list[index][0],tm_list[index][1],tm_list[index][2],tm_list[index][3],TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,args.checkpoint,args.migration_end)
	P=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P11=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	P22=np.zeros(len(TMRCA_bin_psuedo.time_bin))
	t_bin=np.arange(0,0.0065,0.000001)
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_psuedo.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_psuedo.time_bin[j+1]*args.mu:
				P[j]=P[j]+f_density[i]/sum(f_density)
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Af.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Af.time_bin[j+1]*args.mu:
				P11[j]=P11[j]+f_density11[i]/sum(f_density11);
				break
	for i in range(len(t_bin)):
		for j in range(len(P)-1):
			if t_bin[i]>=TMRCA_bin_Eu.time_bin[j]*args.mu and t_bin[i]<=TMRCA_bin_Eu.time_bin[j+1]*args.mu:
				P22[j]=P22[j]+f_density22[i]/sum(f_density22);
				break
	stats1=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi')
	stats2=calc_distance_between_two_dist(TMRCA_bin_Af.genome_percentage,P11,TMRCA_bin_Af.genome_percentage_macs,TMRCA_bin_Af.covariance,'chi')
	stats3=calc_distance_between_two_dist(TMRCA_bin_Eu.genome_percentage,P22,TMRCA_bin_Eu.genome_percentage_macs,TMRCA_bin_Eu.covariance,'chi')
#	stats2=calc_distance_between_two_dist(TMRCA_bin_psuedo.genome_percentage,P,TMRCA_bin_psuedo.genome_percentage_macs,TMRCA_bin_psuedo.covariance,'chi covariance')
	info=list(P)+list(P11)+list(P22)+[stats1]+[stats2]+[stats3]
	info1="\t".join(map(str,info))
	del info
	info2=stats1+stats2+stats3
	return (info1,info2)

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
	for migration_end in [0]:
		for t in range(args.tprior1,args.tprior2+1):
			for m in range(args.mprior1,args.mprior2+1,5):
				tm_list.append([t,m,migration_end])
	print len(tm_list)
	pool = multiprocessing.Pool(processes=30)
	result1=pool.map(calc_f_density_integrate_wrapper_v2,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(len(tm_list))))
	pool.close()
	pool.join()
	TMRCA_dist=[i[0] for i in result1]
	psmc_result1=[i[1] for i in result1]
	psmc_result2=[i[2] for i in result1]
	psmc_result3=[i[3] for i in result1]
	psmc_result4=[i[1]+i[2]+i[3] for i in result1]
	psmc_result5=[i[1]+i[2] for i in result1]
	psmc_result_matrix=psmc_result1
	print len(psmc_result_matrix),1*(args.tprior2-args.tprior1+1)*((args.mprior2-args.mprior1)/5+1)
	psmc_result_matrix1=np.reshape(np.array(psmc_result1),(1,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
	psmc_result_matrix4=np.reshape(np.array(psmc_result4),(1,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
	psmc_result_matrix5=np.reshape(np.array(psmc_result5),(1,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/5+1))
#	psmc_result_matrix2=np.reshape(np.array(psmc_result2),(1,args.tprior2-args.tprior1+1,(args.mprior2-args.mprior1)/50+1))
	print 'psmc chi',psmc_result_matrix1
	print 'psmc chi',psmc_result_matrix4
	print 'psmc chi',psmc_result_matrix5
#	print 'psmc chi covariance',psmc_result_matrix2
	return result1,tm_list

def ABC_SMC_iter(iter,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu,error_value,prev_weight,prev_tm_list):
	tm_list_new=[]
	done = 0
	tm_list_total=[]
	result_total=[]
	result_new=[]
	TMRCA_total =[]
	prior_list=[(args.tprior1,args.tprior2),(args.mprior11,args.mprior12),(args.mprior21,args.mprior22)]
	while done<args.Nparticle:
		if iter==0:
			tm_list=sample_parameters(args.Ntrial,prior_list)
		else:
			tm_list=sample_parameters_weight(args.Ntrial,prior_list,prev_tm_list,prev_weight)
		pool = multiprocessing.Pool(processes=30)
		result1=pool.map(calc_f_density_integrate_wrapper_v3,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(args.Nparticle)))
#		TMRCA_dist=[i[0]+i[1]+i[2]+[i[3]]+[i[4]]+[i[5]] for i in result1]  # 3 TMRCA
		TMRCA_dist=[i[0] for i in result1]
#		print 'TMRCA_dist',len(TMRCA_dist[0])
#		result=[i[3]+i[4]+i[5] for i in result1]	# 3 stats
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
				assert len(tm_list[i])==len(prior_list)
				tm_list_new.append(tm_list[i])
				result_new.append(result[i])
				TMRCA_total.append(TMRCA_dist[i])
				done+=1
		print iter,subprocess.check_output('date').strip(),done
	tm_list_new=tm_list_new[:args.Nparticle]
	result_new=result_new[:args.Nparticle]
	TMRCA_total=TMRCA_total[:args.Nparticle]
	print len(tm_list_new),len(result_new),len(TMRCA_total)
#	print len(tm_list_new[0]),len(TMRCA_total[0])
	new_weight=plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight,TMRCA_total)
	return error_value,new_weight,tm_list_new

def sample_parameters(Ntrial,prior_list):
	sampled_parameters=[]
#	print 'prior_list',len(prior_list)
	for i in range(len(prior_list)):
		sampled_parameters.append(np.random.uniform(prior_list[i][0],prior_list[i][1],Ntrial))
	tm_list=[[sampled_parameters[i][k] for i in range(len(prior_list))] for k in range(Ntrial)]
#	for i in tm_list:
#		print i
#		assert len(i)==len(prior_list)
	return tm_list

def sample_parameters_weight(Ntrial,prior_list,tm_list,weight):
	print 'sum of weight',sum(weight)
	sampled_parameters=[]
	sigma_list=[]
	index = 0
	if args.kernel == "uniform":
		for i in range(len(prior_list)):
			print len(tm_list)
			param=[tm_list[j][i] for j in range(len(tm_list))]
			sigma=0.2*(max(param)-min(param))
			print 'sigma',i,sigma
			sigma_list.append(sigma)
	while index<Ntrial:
		random_number=np.random.uniform(0,1)
		random_id = get_id_from_sample(weight,random_number)
		assert weight[random_id]>0
		param=[]
		for i in range(len(prior_list)):
			while True:
				rand_a=np.random.uniform(0,1)
				rand_direction=np.random.randint(2)
				if rand_direction:
					new_t=tm_list[random_id][i]+rand_a*sigma_list[i]
				else:	
					new_t=tm_list[random_id][i]-rand_a*sigma_list[i]
				if new_t>=prior_list[i][0] and new_t<=prior_list[i][1]:
					param.append(new_t)
					break
		sampled_parameters.append(param)
		index+=1
		if index%100==0:
			print index,'sample done'
	return sampled_parameters

def plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight,TMRCA_dist):
	t=[i[0] for i in tm_list_total]
	m1=[i[1] for i in tm_list_total]
	m2=[i[2] for i in tm_list_total]
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior11,args.mprior12,num=(args.mprior12-args.mprior11)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t,m1,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	fig = plt.figure()
	ax1 = fig.add_subplot(411,title='iter %s sampled param' %(iter))
	cax=ax1.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax1.set_ylabel('migration rate 1')
	t_accepted=[i[0] for i in tm_list_new]
	m1_accepted=[i[1] for i in tm_list_new]
	m2_accepted=[i[2] for i in tm_list_new]
	print 't,m',np.mean(t_accepted),np.mean(m1_accepted),np.mean(m2_accepted)
	ax2 = fig.add_subplot(412,title='iter %s accepted param' %(iter))
	heatmap, xedges, yedges = np.histogram2d(t_accepted, m1_accepted,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax2.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax3 = fig.add_subplot(413)
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior11,args.mprior12,num=(args.mprior22-args.mprior21)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t,m2,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax3.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax3.set_ylabel('migration rate 2')
	ax4 = fig.add_subplot(414,title='iter %s accepted param' %(iter))
	heatmap, xedges, yedges = np.histogram2d(t_accepted, m2_accepted,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	cax=ax4.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax4.set_xlabel('t')
	fig.savefig("version%s/iter%s_sampled_accepted_param.pdf" %(args.version,iter),format='pdf')

	fig = plt.figure()
	ax1 = fig.add_subplot(411,title='iter %s stats' %(iter))
	bins = 40
	ax1.hist(result_total,facecolor='green',alpha=0.5,histtype='stepfilled',label='genome percentage')
	ax1.legend(prop={'size':10})
	ax1.set_ylabel('count')
	ax1.axvline(x=e,linewidth=2,color='k')

	ax2 = fig.add_subplot(412)
	cax=ax2.scatter(t_accepted,m1_accepted,c=result_new,cmap=plt.cm.coolwarm)
	fig.colorbar(cax,orientation='vertical')
	ax2.set_ylabel('migration1')
	ax3 = fig.add_subplot(413)
	cax=ax3.scatter(t_accepted,m2_accepted,c=result_new,cmap=plt.cm.coolwarm)
	fig.colorbar(cax,orientation='vertical')
	ax3.set_ylabel('migration2')
	if iter==0:
		new_weight = [1./len(tm_list_new) for i in range(len(tm_list_new))]
	else:
		prior_list=[(args.tprior1,args.tprior2),(args.mprior11,args.mprior12),(args.mprior21,args.mprior22)]
		sigma_list=[]
		if args.kernel == "uniform":
			for i in range(len(prior_list)):
				param=[prev_tm_list[j][i] for j in range(len(prev_tm_list))]
				sigma=0.2*(max(param)-min(param))
				print 'sigma',i,sigma
				sigma_list.append(sigma)
			new_weight=update_weight(iter,prev_tm_list,tm_list_new,prev_weight,sigma_list)
	ax4 = fig.add_subplot(414)
	cax=ax4.scatter(m1_accepted,m2_accepted,c=result_new,cmap=plt.cm.coolwarm)
	ax4.set_xlabel('migration1')
	ax4.set_ylabel('migration2')
	fig.colorbar(cax,orientation='vertical')	
	fig.savefig("version%s/iter%s_distance_stats_scatterplot.pdf" %(args.version,iter),format='pdf')

	fout = open('version%s/iter%s_stats_weight.txt' %(args.version,iter),'w')
	for i in range(len(new_weight)):
#		print >>fout,'%i\t%s\t%f\t%f\t%s' %(i,"\t".join(map(str,tm_list_new[i])),result_new[i],new_weight[i],"\t".join(map(str,TMRCA_dist[i])))
		print >>fout,'%i\t%s\t%f\t%f\t%s' %(i,"\t".join(map(str,tm_list_new[i])),result_new[i],new_weight[i],TMRCA_dist[i])	
	fout.close()
	return new_weight

def update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_list):
	new_weight = np.zeros(len(new_tm_list))
	for i in range(len(new_weight)):
		for j in range(len(prev_weight)):
			within=True
			for k in range(len(sigma_list)):
				if abs(prev_tm_list[j][k]-new_tm_list[i][k])>sigma_list[k]:
					within=False
					break
			if within:
#				new_weight[i]+=prev_weight[j]*1/(2*sigma_t*2*sigma_m)
				new_weight[i]+=prev_weight[j]
	new_weight = [1/float(i) for i in new_weight]
	new_weight = [float(i)/sum(new_weight) for i in new_weight]
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

def print_TMRCA_result_v2(tm_list,result):
	f=open('sim_macs_YRI_CEU_8_20_grid_TMRCA_3_stats.txt','w')
	for i in range(len(tm_list)):
		print >>f,"%s\t%s\t%s" %("\t".join(map(str,tm_list[i])),"\t".join(map(str,result[i][0])),"\t".join(map(str,result[i][1:])))

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

	TMRCA_bin_pseudo=TMRCA_BIN(label='YRI_CEU',time_bin=[],genome_percentage=[],genome_percentage_formula=[],pop_size=[],pop_size_formula=[])  # the pseudo diploid based on real data(the one you observe)
	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,None,4)
#	TMRCA_bin_pseudo.read_from_psmc_mulitple(args.file3,['0','1','2','3'],4)
#	TMRCA_bin_pseudo.calc_covariance()
	TMRCA_bin_pseudo.print_TMRCA_bin()
	TMRCA_bin_Af=TMRCA_BIN(label='YRI',time_bin=[],genome_percentage=[],genome_percentage_formula=[],pop_size=[],pop_size_formula=[])
	if args.type=='psmc':
		TMRCA_bin_Af.read_from_psmc_mulitple(args.file1,None,1)
	elif args.type=='msmc':
		TMRCA_bin_Af.read_from_msmc(args.file1)
	TMRCA_bin_Af.print_TMRCA_bin()
#	TMRCA_bin_Af.calc_TMRCA_dist_formula_v2()
#	TMRCA_bin_Af.calc_theta_t_from_TMRCA_v2()
#	TMRCA_bin_Af.draw_TMRCA()
	TMRCA_bin_Eu=TMRCA_BIN(label='CEU',time_bin=[],genome_percentage=[],genome_percentage_formula=[],pop_size=[],pop_size_formula=[])
	print 'start'
	TMRCA_bin_Eu.print_TMRCA_bin()
	if args.type=='psmc':
		TMRCA_bin_Eu.read_from_psmc_mulitple(args.file2,None,1)
	elif args.type=='msmc':
		print args.file2
		TMRCA_bin_Eu.read_from_msmc(args.file2,verbose=True)
	TMRCA_bin_Eu.print_TMRCA_bin()
#	TMRCA_bin_Eu.calc_TMRCA_dist_formula_v2()
#	TMRCA_bin_Eu.calc_theta_t_from_TMRCA_v2()
#	TMRCA_bin_Eu.draw_TMRCA()

	'''
	f_density_array=[]
	P_array=[]
	for t in range(7,11):
		f_density,f_density11,f_density22,TMRCA_bin_11,TMRCA_bin_12,TMRCA_bin_22=calc_relative_cross_coalescent_curve(t,15,4,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
		draw_TMRCA_and_coalescence_rate(t,15,4,TMRCA_bin_11,TMRCA_bin_12,TMRCA_bin_22)
		draw_f_density_compare_v4(t,15,4,[f_density,f_density11,f_density22],[TMRCA_bin_12.genome_percentage,TMRCA_bin_11.genome_percentage,TMRCA_bin_22.genome_percentage],TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
		f_density_array.append(f_density)
		P_array.append(TMRCA_bin_12.genome_percentage)

	draw_f_density_compare_v3(7,15,4,f_density_array,P_array,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
	'''
#	psmc_result_matrix,tm_list=explore_stats_in_grid_verbose(TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print_TMRCA_result_v2(tm_list,psmc_result_matrix)
	# Test function
#	f_density,P=calc_f_density_integrate_single(10,20,2,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	print P1
#	print P2
#	draw_f_density_compare(6.24,987,args.migration_end,f_density1,P1,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)
#	draw_f_density_compare_v3(10,20,2,f_density,P,TMRCA_bin_pseudo,TMRCA_bin_Af,TMRCA_bin_Eu)


	## PERFORM ABC
	if args.initial_e is None:
		error_value = []
		for i in range(args.last_iter):
			error_value.append(0)
	else:
		error_value = args.initial_e
	prev_weight = [1 for i in range(args.Nparticle)]
	prev_tm_list = [[0,0,0] for i in range(args.Nparticle)]
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

