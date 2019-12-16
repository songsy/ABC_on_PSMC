#!/usr/bin/env python
# python ABC_calc_TMRCA_v2.py
# Shiya Song
# 17th April 2015

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

def runCMD(cmd):
	val = subprocess.Popen(cmd, shell=True,bufsize=0).wait()
	if val == 0:
		pass
	else:
		print 'command failed'
		exit(1)

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

def sample_paramters(Ntrial,tprior,mprior):
	tprior_list=np.random.uniform(tprior[0],tprior[1],Ntrial)
	mprior_list=np.random.uniform(mprior[0],mprior[1],Ntrial)
	return tprior_list,mprior_list

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
	t=[tm_list[i][0] for i in range(len(tm_list)) if weight[i]]
	m=[tm_list[i][1] for i in range(len(tm_list)) if weight[i]]
	index = 0
	if args.kernel == "uniform":
		sigma_t=0.5*(max(t)-min(t))
		sigma_m=0.5*(max(m)-min(m))
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
	return t_list,m_list

def msmc_to_macs_v3(file1,file2,file3,t,m,migration_end,CHECKPOINT): # no exponential growth, all psmc result
	theta = 4*int(args.simN0)*float(args.mu)*int(args.L)
	pho = theta*0.13
	split_time = float(t*10000)/4/int(args.simN0)/G
	migration_rate = m/1e05*4*int(args.simN0)
	print t,m,split_time,migration_rate
	check_point = float(CHECKPOINT)/4/int(args.simN0)/G
#	print check_point
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.simN0)
		size = float(i[1])*10000/int(args.simN0)
		start_size1 = size
		break
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.simN0)
		size = float(i[1])*10000/int(args.simN0)
		start_size2 = size
		break	
	cmds += '-n 1 %f -n 2 %f ' %(start_size1,start_size2)
	interval1 = []
	interval2 = []
	share_interval = []
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.simN0)
		size = float(i[1])*10000/int(args.simN0)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
#	print interval1
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.simN0)
		size = float(i[1])*10000/int(args.simN0)
		if time<=split_time and time>0.005:
			interval2.append([time,size])
#	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.simN0)
		size = float(i[1])*10000/int(args.simN0)
		if i[0][0]=='t':
			continue
		if float(i[0])<CHECKPOINT:
			continue
		if t*10000> CHECKPOINT:
			if float(i[0])< t*10000:
				continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
#		if i<len(interval1):
#			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
#	cmds += '-en %s 1 %f ' %(split_time+0.000001,start_size1)          ## change it! ##
	cmds += '-en %s 1 %f ' %(split_time+0.000001,share_interval[0][1])
#	for i in range(len(interval1_inbetween)):
#		cmds += '-en %f 1 %f ' %(interval1_inbetween[i][0],interval1_inbetween[i][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.simN0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time+0.000002)
	return cmds

def print_macs_cmds(iter,Ntrial,tprior,mprior):
	f_out=open('macs_iter%s.cmds' %(iter),'w')
	f_summary=open('macs_iter%s_param.txt' %(iter),'w')
	for i in range(Ntrial):
		t=round(tprior[i],1)
		m=round(mprior[i],1)
		print >>f_summary,'%.1f\t%.1f' %(t,m)
		cmds=msmc_to_macs_v3(args.file1,args.file2,args.file3,t,m,int(args.migration_end),int(args.checkpoint))
		for j in range(args.nchr):
			new_cmds="%s-T 2>iter%i/sim_macs_%i_%i_%.1f_%.1f_%i_trees.txt 1>iter%i/sim_macs_%i_%i_%.1f_%.1f_%i_sites.txt" %(cmds,iter,iter,i,t,m,j,iter,iter,i,t,m,j)
			new_cmds = new_cmds.replace("SEED",str(142084502+iter*19+i*13+j*17))
			print >>f_out, new_cmds
		print i,'done'

def rewrite_macs_cmds(iter):
	f_summary=open('macs_iter%s_param.txt' %(iter),'r')
	f_out=open('macs_iter%s_new.cmds' %(iter),'w')
	for i,line in enumerate(f_summary):
		col=line.rstrip().split('\t')
		t=float(col[0])
		m=float(col[1])
		cmds=msmc_to_macs_v3(args.file1,args.file2,args.file3,t,m,int(args.migration_end),int(args.checkpoint))
		for j in range(args.nchr):
			new_cmds="%s-T 2>iter%i/sim_macs_%i_%i_%.1f_%.1f_%i_trees.txt|python /home/songsy/script/pop-gen/ABC_TMRCA_from_macs_v1.py --file iter%i/sim_macs_%i_%i_%.1f_%.1f --psmc ../../psmc/simulation/macs/sim_macs_Gravel_2pop_no_mig --index %i" %(cmds,iter,iter,i,t,m,j,iter,iter,i,t,m,j)
			new_cmds = new_cmds.replace("SEED",str(142084502+iter*19+i*13+j*17))
			print >>f_out, new_cmds
#		print i,'done'

def rerun_macs_cmds(iter):
	f_summary=open('macs_iter%s_param.txt' %(iter),'r')
	for i,line in enumerate(f_summary):
		col=line.rstrip().split('\t')
		t=float(col[0])
		m=float(col[1])
#		cmds=msmc_to_macs_v3(args.file1,args.file2,args.file3,t,m,int(args.migration_end),int(args.checkpoint))
		for j in range(args.nchr):
			if os.path.isfile('iter%i/sim_macs_%i_%i_%.1f_%.1f_TMRCA_%i.txt' %(iter,iter,i,t,m,j)) is False:
				os.popen("grep sim_macs_%i_%i_%.1f_%.1f macs_iter%s_new.cmds|grep 'index %i' >> macs_iter%s_redo.cmds" %(iter,i,t,m,iter,j,iter))
		print i,'done'

def run_pbs_macs(iter,Ntrial):
	os.popen('cp run-macs-iter.pbs run-macs-iter%s.pbs' %(iter))
	os.popen("sed -i 's/ITER/iter%s/g' run-macs-iter%s.pbs" %(iter,iter))
	os.popen("sed -i 's/NUM/%s/g' run-macs-iter%s.pbs" %(Ntrial/10,iter))
	output=subprocess.check_output("qsub -N macs_iter%s run-macs-iter%s.pbs" %(iter,iter),shell=True)
	job_id=output.rstrip().split('[')[0]
	print 'job submitted',job_id
	job_name='macs_iter%s' %(iter)
	return job_name

def run_pbs_macs_to_psmcfa(iter,Ntrial):
	os.popen('cp run-macs-to-psmcfa-iter.pbs run-macs-to-psmcfa-iter%s.pbs' %(iter))
	os.popen("sed -i 's/ITER/iter%s/g' run-macs-to-psmcfa-iter%s.pbs" %(iter,iter))
	os.popen("sed -i 's/NUM/%s/g' run-macs-to-psmcfa-iter%s.pbs" %(Ntrial/10,iter))
	output=subprocess.check_output("qsub -N macs_fa_iter%s run-macs-to-psmcfa-iter%s.pbs" %(iter,iter),shell=True)
	job_id=output.rstrip().split('[')[0]
	print 'job submitted',job_id
	job_name='macs_fa_iter%s' %(iter)
	return job_name

def run_pbs_psmc(iter,Ntrial):
	os.popen('cp run-psmc-iter-0.pbs run-psmc-iter%s-0.pbs' %(iter))
	os.popen('cp run-psmc-iter-1.pbs run-psmc-iter%s-1.pbs' %(iter))
	os.popen("sed -i 's/ITER/iter%s/g' run-psmc-iter%s-0.pbs" %(iter,iter))
	os.popen("sed -i 's/ITER/iter%s/g' run-psmc-iter%s-1.pbs" %(iter,iter))
	output=subprocess.check_output("qsub -N psmc_iter%s run-psmc-iter%s-0.pbs" %(iter,iter),shell=True)
	job_id=output.rstrip().split('[')[0]
	print 'job submitted',job_id
	time.sleep(3600)
	output=subprocess.check_output("qsub -N psmc_iter%s run-psmc-iter%s-1.pbs" %(iter,iter),shell=True)
	job_id=output.rstrip().split('[')[0]
	print 'job submitted',job_id
	job_name='psmc_iter%s' %(iter)
	time.sleep(3600*4)
	return job_name

def check_job(iter,job_name,sleep_time):
	while True:
		time.sleep(sleep_time)
		process=subprocess.Popen('qstat -t -u songsy|grep %s' %(job_name),shell=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
		output,err=process.communicate()
		status=output.split('\n')
		finish=True
		for line in status:
			if line=='':
				continue
			col=line.split()[-2]
			print line,col
			if col!='C':
				finish=False
				break
		print finish
		if finish:
			break

def get_TMRCA_from_macs_v3(file,TMRCA_bin,start):
	TMRCA=[[] for x in range(COPY)]
	prev_index_i=[0 for x in range(COPY)]
	prev_TMRCA=[0 for x in range(COPY)]
	size_chunk=[0 for x in range(COPY)]
	chunk_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)] 
	genome_percentage_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	average_length_bin=[[0 for x in range(len(TMRCA_bin[0]))] for x in range(COPY)]
	print len(genome_percentage_bin),len(genome_percentage_bin[0])
	for i in range(args.nchr):
		first = [True for x in range(COPY)]
		done = [True for x in range(COPY)]
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
#				print size,newick
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
						index +=1
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

def check_macs_result(iter):
	f=open('macs_iter%s_param.txt' %(iter),'r')
	fout = open('calc_TMRCA_from_macs_iter%s.cmds' %(iter),'w')
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		t=line[0]
		m=line[1]
		file='sim_macs_%i_%i_%s_%s' %(iter,index,t,m)
		cmds = 'python /home/songsy/script/pop-gen/ABC_TMRCA_from_macs.py --file iter%s/%s --psmc %s' %(iter,file,args.psmc)
		print >>fout,cmds
	fout.close()

def print_macs_to_psmcfa_cmds(iter):
	f=open('macs_iter%s_param.txt' %(iter),'r')
	fout1 = open('macs_to_psmcfa_iter%s.cmds' %(iter),'w')
	fout2 = open('psmc_iter%s.cmds' %(iter),'w')
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		t=line[0]
		m=line[1]
		file='sim_macs_%i_%i_%s_%s' %(iter,index,t,m)
		cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-macs.py --file iter%s/%s --nrep 100 --nchr 100' %(iter,file)
		print >>fout1,cmds
		cmds = 'psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o iter%s/%s.psmc iter%s/%s_0_hetFQ.0' %(iter,file,iter,file)
		print >>fout2,cmds
	fout1.close()
	fout2.close()

def run_pbs_calc(iter,Ntrial):
	os.popen('cp run-calc-iter.pbs run-calc-iter%s.pbs' %(iter))
	os.popen("sed -i 's/ITER/iter%s/g' run-calc-iter%s.pbs" %(iter,iter))
	os.popen("sed -i 's/NUM/%s/g' run-calc-iter%s.pbs" %(Ntrial,iter))
	output=subprocess.check_output("qsub -N calc_iter%s run-calc-iter%s.pbs" %(iter,iter),shell=True)
	job_id=output.rstrip().split('[')[0]
	print 'job submitted',job_id
	return job_id

def get_stats(iter,Ntrial,genome_percentage,chunk_percentage,interval_compare):
	f=open('macs_iter%s_param.txt' %(iter),'r')
	fout = open('macs_iter%s_stats.txt' %(iter),'w')
	tm_list = []
	diff_genome_list=[]
	diff_chunk_list=[]
	for index,line in enumerate(f):
		diff_genome = [0 for x in range(COPY)]
		diff_chunk = [0 for x in range(COPY)]
		line=line.rstrip().split('\t')
		t=line[0]
		m=line[1]
		file='iter%s/sim_macs_%i_%i_%s_%s_TMRCA.txt' %(iter,iter,index,t,m)
		f_stats=open(file,'r')
		for line_index,each_line in enumerate(f_stats):
			each_line=each_line.rstrip()
			col = each_line.split('\t')
			if line_index > interval_compare-1:
				break
			for j in range(4):
				diff_genome[j]+=abs(genome_percentage[j][line_index]-float(col[j*2+1]))
				diff_chunk[j]+=abs(chunk_percentage[j][line_index]-float(col[j*2+2]))
		mean_diff_genome=float(np.mean(diff_genome))/2
		mean_diff_chunk=float(np.mean(diff_chunk))/2
		print >>fout,'%i\t%s\t%s\t%f\t%f' %(index,t,m,mean_diff_genome,mean_diff_chunk)
		if index%100==0:
			print index,'done'
		tm_list.append([float(t),float(m)])
		diff_genome_list.append(mean_diff_genome)
		diff_chunk_list.append(mean_diff_chunk)
		f_stats.close()
	fout.close()
	return tm_list,diff_genome_list,diff_chunk_list

def plot_stats(iter,new_tm_list,prev_tm_list,diff_genome_list,diff_chunk_list,error_value,prev_weight):
	fig = plt.figure()
	ax1 = fig.add_subplot(211,title='iter %s stats' %(iter))
	ax2 = fig.add_subplot(212)
	bins = 40
	ax1.hist(diff_genome_list,facecolor='green',alpha=0.5,histtype='stepfilled',label='genome percentage')
	ax1.legend(prop={'size':10})
	ax1.set_xlabel('variation distance based on genome percentage')
	ax1.set_ylabel('count')
#	fig.suptitle(sample_name, fontsize=10)
	ax2.hist(diff_chunk_list,facecolor='blue',alpha=0.5,histtype='stepfilled',label='chunk percentage')
	ax2.legend(prop={'size':10})
	ax2.set_xlabel('variation distance based on chunk percentage')
	ax2.set_ylabel('count')
	'''
	if iter<len(error_value):
		ax2.axvline(x=error_value[iter],linewidth=2,color='k')
		accept_indicator=[1 if i<=error_value[iter] else 0 for i in diff_chunk_list]
		print 'iter%s error value:%f' %(iter,error_value[iter])
	else:
		e = np.percentile(diff_chunk_list,50)
		error_value.append(e)
		ax2.axvline(x=error_value[iter],linewidth=2,color='k')
		accept_indicator=[1 if i<=error_value[iter] else 0 for i in diff_chunk_list]
		print 'iter%s error value:%f' %(iter,error_value[iter])
	fig.savefig("iter%s_distance_stats.pdf" %(iter),format='pdf')
	'''
	if iter<=1:
		e = np.percentile(diff_chunk_list,50)
	else:
		e = np.percentile(diff_chunk_list,20)
	error_value.append(e)
	ax2.axvline(x=error_value[iter],linewidth=2,color='k')
	accept_indicator=[1 if i<=error_value[iter] else 0 for i in diff_chunk_list]
	print 'iter%s error value:%f' %(iter,error_value[iter])
	fig.savefig("iter%s_distance_stats.pdf" %(iter),format='pdf')

	t=[i[0] for i in new_tm_list]
	m=[i[1] for i in new_tm_list]
	print accept_indicator
	print accept_indicator.count(1),accept_indicator.count(0)
	xedges=list(np.linspace(args.tprior1,args.tprior2,num=(args.tprior2-args.tprior1)/0.5+1))
	yedges=list(np.linspace(args.mprior1,args.mprior2,num=(args.mprior2-args.mprior1)/0.5+1))
	heatmap, xedges, yedges = np.histogram2d(t, m,bins=(xedges,yedges))
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s sampled param' %(iter))
	cax=ax1.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	fig.colorbar(cax,orientation='vertical')
	ax1.set_ylabel('migration rate')
	t_accepted=[new_tm_list[i][0] for i in range(len(new_tm_list)) if accept_indicator[i]]
	m_accepted=[new_tm_list[i][1] for i in range(len(new_tm_list)) if accept_indicator[i]]
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
#	plt.tight_layout()
	fig.savefig("iter%s_sampled_accepted_param.pdf" %(iter),format='pdf')
	
	fig = plt.figure()
	ax1 = fig.add_subplot(311,title='iter %s stats' %(iter))
	cax=ax1.scatter(t,m,c=diff_chunk_list,cmap=plt.cm.coolwarm)
	ax1.set_xlabel('variation distance based on genome percentage')
	fig.colorbar(cax,orientation='vertical')
	ax2 = fig.add_subplot(312)
	cax=ax2.scatter(t,m,c=diff_genome_list,cmap=plt.cm.coolwarm)
	ax2.set_xlabel('variation distance based on chunk percentage')
	fig.colorbar(cax,orientation='vertical')
	print 'figure saved'
	if iter==0:
		new_weight = accept_indicator
		new_weight = [float(i)/sum(new_weight) for i in new_weight]
	else:
		t_prev=[prev_tm_list[i][0] for i in range(len(prev_tm_list)) if prev_weight[i]]
		m_prev=[prev_tm_list[i][1] for i in range(len(prev_tm_list)) if prev_weight[i]]
		print 't_prev',t_prev
		print 'm_prev',m_prev
		if args.kernel == "uniform":
			sigma_t=0.5*(max(t_prev)-min(t_prev))
			sigma_m=0.5*(max(m_prev)-min(m_prev))
			print 'sigma_t', sigma_t
			print 'sigma_m', sigma_m
			new_weight=update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_t,sigma_m,accept_indicator)
	ax3 = fig.add_subplot(313)
	cax=ax3.scatter(t,m,c=new_weight,cmap=plt.cm.coolwarm)
	ax3.set_xlabel('weight')
	fig.colorbar(cax,orientation='vertical')	
	fig.savefig("iter%s_distance_stats_scatterplot.pdf" %(iter),format='pdf')	
	return new_weight,error_value

def update_weight(iter,prev_tm_list,new_tm_list,prev_weight,sigma_t,sigma_m,accept_indicator):
	new_weight = accept_indicator
	for i in range(len(new_weight)):
		if new_weight[i]==1:
			new_weight[i] =0 
			for j in range(len(prev_weight)):
				if prev_weight[j]:
					if abs(prev_tm_list[j][0]-new_tm_list[i][0])<=sigma_t and abs(prev_tm_list[j][1]-new_tm_list[i][1])<=sigma_m:
						new_weight[i]+=prev_weight[j]*1/(2*sigma_t*2*sigma_m)
	new_weight = [float(i)/sum(new_weight) for i in new_weight]
	return new_weight

def write_output_v2():
	for index in range(COPY):
		fout=open(args.psmc+'_TMRCA_ave_length_%s_test.txt' %(index),'w')
		for i in range(len(TMRCA_bin[index])):
			print >>fout,'%f\t%f\t%f\t%f' %(TMRCA_bin[index][i],genome_percentage[index][i],scale_factor[index][i],float(chunk_percentage[index][i])/sum(chunk_percentage[index]))
	fout.close()

def calc_stats_from_psmc(iter,genome_percentage,chunk_percentage,interval_compare):
	f=open('macs_iter%s_param.txt' %(iter),'r')
	fout = open('macs_iter%s_stats.txt' %(iter),'w')
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
		print >>fout,'%i\t%s\t%s\t%f\t%f' %(index,t,m,mean_diff_genome,mean_diff_chunk)
		if index%100==0:
			print index,'done'
		tm_list.append([float(t),float(m)])
		diff_genome_list.append(mean_diff_genome)
		diff_chunk_list.append(mean_diff_chunk)
		f_stats.close()
	fout.close()
	return tm_list,diff_genome_list,diff_chunk_list

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--psmc", dest='psmc',help="psmc file prefix")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--N0", default=5e04*2.36/1.25,dest='N0',help="N0")
	parser.add_argument("--simN0", default=5e04,dest='simN0',help="N0")
	parser.add_argument("--Ntrial", default=500,type=int,dest='Ntrial',help="Number of trials per iteration")
	parser.add_argument("--Nparticle", default=200,type=int,dest='Nparticle',help="number of particles to reach")
	parser.add_argument("--initial_e",type=float,nargs='+',dest='initial_e',help="inital error to tolerance")
	parser.add_argument("--nchr", default=100,type=int,dest='nchr',help="number of chromosome per simulation")
	parser.add_argument("--nhap", default=2,type=int,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", type=int,default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
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
	parser.add_argument("--log", dest='log',help="logfile")
	args = parser.parse_args()

	sys.stdout = open(args.log,'a',0)
	COPY = 4
	G=30
	ITERATION=5

	print subprocess.check_output('date').strip()
	print 'get TMRCA bin'
	TMRCA_bin,genome_percentage,rho=get_TMRCA_bin(args.psmc)
	print 'get chunk percentage'
	chunk_percentage,scale_factor=calc_chunk_percentage(genome_percentage,TMRCA_bin,rho,COPY)
	print 'get TMRCA from macs'

	write_output_v2()

	print subprocess.check_output('date').strip()
	if args.initial_e is None:
		error_value = []
	else:
		error_value = args.initial_e
	prev_weight = [1 for i in range(args.Ntrial)]
	prev_weight = [float(i)/sum(prev_weight) for i in prev_weight]
	prev_tm_list = [[0,0] for i in range(args.Ntrial)]
	for iter in range(ITERATION):
		print 'iter%s begin' %(iter)
		print subprocess.check_output('date').strip()
		if iter==0:
			'''
			if os.path.isdir('iter%i' %(iter)) is False:
				os.popen('mkdir iter%i' %(iter))
			tprior_list,mprior_list=sample_paramters(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2))
			print 'sampling iter%i done' %(iter)
			print_macs_cmds(iter,args.Ntrial,tprior_list,mprior_list)
			job_name=run_pbs_macs(iter,args.Ntrial)
			print 'running macs'
			print subprocess.check_output('date').strip()
			sys.stdout.flush()
			check_job(iter,job_name,60)
			print subprocess.check_output('date').strip()
			print_macs_to_psmcfa_cmds(iter)
			job_name=run_pbs_macs_to_psmcfa(iter,args.Ntrial)
			sys.stdout.flush()
			print 'running macs to psmcfa'
			print subprocess.check_output('date').strip()
			check_job(iter,job_name,60)
			print subprocess.check_output('date').strip()
			job_name=run_pbs_psmc(iter,args.Ntrial)
			print 'running psmc'
			print subprocess.check_output('date').strip()
			check_job(iter,job_name,60)
			'''
			if iter>=1:
				prev_weight=new_weight
				prev_tm_list=new_tm_list
			new_tm_list,diff_genome_list,diff_chunk_list=calc_stats_from_psmc(iter,genome_percentage,chunk_percentage,args.interval_compare)
			new_weight,error_value=plot_stats(iter,new_tm_list,prev_tm_list,diff_genome_list,diff_chunk_list,error_value,prev_weight)
		else:
			if os.path.isdir('iter%i' %(iter)) is False:
				os.popen('mkdir iter%i' %(iter))
			tprior_list,mprior_list=sample_parameters_weight(args.Ntrial,(args.tprior1,args.tprior2),(args.mprior1,args.mprior2),new_tm_list,new_weight)
			print 'sampling iter%i done' %(iter)
			print subprocess.check_output('date').strip()
			print_macs_cmds(iter,args.Ntrial,tprior_list,mprior_list)
			job_name=run_pbs_macs(iter,args.Ntrial)
			print 'running macs'
			print subprocess.check_output('date').strip()
			sys.stdout.flush()
			check_job(iter,job_name,60)
			print subprocess.check_output('date').strip()
			print_macs_to_psmcfa_cmds(iter)
			job_name=run_pbs_macs_to_psmcfa(iter,args.Ntrial)
			sys.stdout.flush()
			print 'running macs to psmcfa'
			print subprocess.check_output('date').strip()
			check_job(iter,job_name,60)
			print subprocess.check_output('date').strip()
			job_name=run_pbs_psmc(iter,args.Ntrial)
			print 'running psmc'
			print subprocess.check_output('date').strip()
			check_job(iter,job_name,60)

			prev_weight=new_weight
			prev_tm_list=new_tm_list
#			print 'prev_weight',prev_weight
#			print 'prev_tm_list',prev_tm_list
			new_tm_list,diff_genome_list,diff_chunk_list=calc_stats_from_psmc(iter,genome_percentage,chunk_percentage,args.interval_compare)
			new_weight,error_value=plot_stats(iter,new_tm_list,prev_tm_list,diff_genome_list,diff_chunk_list,error_value,prev_weight)
		print 'iter%s ends' %(iter)
		print subprocess.check_output('date').strip()


