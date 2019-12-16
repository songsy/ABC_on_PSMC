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

def calc_f_density_integrate_wrapper(args):
	args2 = args[0] + (args[1],)
	stats=calc_f_density_integrate(*args2)
	return stats

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
	while done<args.Nparticle:
		if iter==0:
			tm_list=sample_parameters(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2))
		else:
			tm_list=sample_parameters_weight(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2),prev_tm_list,prev_weight)
		pool = multiprocessing.Pool(processes=30)
		result=pool.map(calc_f_density_integrate_wrapper,itertools.izip(itertools.repeat((tm_list,TMRCA_bin_psuedo,TMRCA_bin_Af,TMRCA_bin_Eu)),range(args.Nparticle)))
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
				done+=1
		print iter,subprocess.check_output('date').strip(),done
	tm_list_new=tm_list_new[:args.Nparticle]
	result_new=result_new[:args.Nparticle]
	new_weight=plot_stats(iter,tm_list_total,tm_list_new,prev_tm_list,result_total,result_new,e,prev_weight)
	return error_value,new_weight,tm_list_new

def update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_t,sigma_m):
	'''
	new_weight = np.zeros(len(new_tm_list))
	for i in range(len(new_weight)):
		for j in range(len(prev_weight)):
			if abs(prev_tm_list[j][0]-new_tm_list[i][0])<=sigma_t and abs(prev_tm_list[j][1]-new_tm_list[i][1])<=sigma_m:
				new_weight[i]+=prev_weight[j]*1/(2*sigma_t*2*sigma_m)
#	new_weight = [float(i)/sum(new_weight) for i in new_weight]
	new_weight = [1/float(i) for i in new_weight]
	new_weight = [float(i)/sum(new_weight) for i in new_weight]
	'''
	new_weight = np.zeros(len(new_tm_list))
	new_weight=[1./len(new_weight) for i in new_weight]
	return new_weight

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
	#	print '# passed:',len(t),len(m)
	#	print 'sigma_t', sigma_t
	#	print 'sigma_m', sigma_m
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

def calc_stats_from_psmc(iter,genome_percentage,chunk_percentage,interval_compare):
	f=open('macs_iter%s_param.txt' %(iter),'r')
#	fout = open('macs_iter%s_stats.txt' %(iter),'w')
	tm_list = []
	diff_genome_list=[]
	diff_chunk_list=[]
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		t=line[0]
		m=line[1]
		diff_genome = [0 for x in range(COPY)]
		diff_chunk = [0 for x in range(COPY)]
		file='iter%s/sim_macs_%i_%i_%s_%s.psmc' %(iter,iter,index,t,m)
		print file
		TMRCA_bin1,genome_percentage1,rho1=get_TMRCA_bin_single_psmc(file)
		chunk_percentage1,scale_factor1=calc_chunk_percentage(genome_percentage1,TMRCA_bin1,rho1,1)
		f_stats=open(file+'.txt','w')
		for i in range(len(TMRCA_bin1[0])):
			print >>f_stats,'%f\t%f\t%f\t%f' %(TMRCA_bin1[0][i],genome_percentage1[0][i],scale_factor1[0][i],chunk_percentage1[0][i])
			if i < interval_compare:
				for j in range(COPY):
					diff_genome[j]+=abs(genome_percentage[j][i]-genome_percentage1[0][i])
					diff_chunk[j]+=abs(chunk_percentage[j][i]-chunk_percentage1[0][i])
		mean_diff_genome=float(np.mean(diff_genome))/2
		mean_diff_chunk=float(np.mean(diff_chunk))/2
		f_stats.close()
#		print >>fout,'%i\t%s\t%s\t%f\t%f' %(index,t,m,mean_diff_genome,mean_diff_chunk)
		if index%100==0:
			print index,'done'
		tm_list.append([float(t),float(m)])
		diff_genome_list.append(mean_diff_genome)
		diff_chunk_list.append(mean_diff_chunk)
		f_stats.close()
#	fout.close()
	return tm_list,diff_genome_list,diff_chunk_list

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
	parser.add_argument("--checkpoint", default=100,dest='checkpoint',help="checkpoint in yr")
	parser.add_argument("--migration_end", default=100,dest='migration_end',help="migration end in kyr")
	parser.add_argument("--interval_compare", type=int,default=13,dest='interval_compare',help="number of first intervals to compare")
	parser.add_argument("--last_iter", type=int,default=0,dest='last_iter',help="last_iter")
	parser.add_argument("--log", dest='log',help="logfile")
	args = parser.parse_args()

#	sys.stdout = open(args.log,'a',0)

	ITERATION=4

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'


#	TMRCA_bin=TMRCA_BIN()
#	TMRCA_bin.read_from_psmc_mulitple(args.psmc,None,1)
#	TMRCA_bin.print_TMRCA_bin()
#	f_density=calc_f_density(10,10,TMRCA_bin,TMRCA_bin,TMRCA_bin)
#	draw_f_density(f_density)

	TMRCA_bin=TMRCA_BIN()
	TMRCA_bin.read_from_TMRCA_txt(args.file)

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
			error_value,new_weight,new_tm_list=ABC_SMC_iter(iter,TMRCA_bin,TMRCA_bin,TMRCA_bin,error_value,prev_weight,prev_tm_list)
		else:
			prev_weight=new_weight
			prev_tm_list=new_tm_list
			error_value,new_weight,new_tm_list=ABC_SMC_iter(iter,TMRCA_bin,TMRCA_bin,TMRCA_bin,error_value,prev_weight,prev_tm_list)
