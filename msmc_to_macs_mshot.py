#!/usr/bin/env python
# python msmc_to_macs.py
# Shiya Song
# 9 Jan 2015

#python /home/songsy/script/pop-gen/msmc_to_macs_YRI_MKK.py --file1 ../YRI/NA19240_gvcf_msmc.final.txt --file2 ../MKK/NA21302_gvcf_msmc.final.txt
# python /home/songsy/script/pop-gen/msmc_to_macs.py --file1 ../../MKK/NA21302_gvcf_msmc.final.txt --file2 ../../CEU/NA12878_1KG_msmc.final.txt --file3 ../../African_CEU/psmcfa/MKK-CEU-true-fosmid-psmc.3.txt --pop MKK --migration_end 50 --checkpoint 105660 -L 1e08 -nrep 3000
import argparse
from NGS_utils import *
import math

G=30

def msmc_to_macs(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*1
	split_time = float(t*10000)/4/int(args.N0)/G
	check_point = float(CHECKPOINT)/4/int(args.N0)/G
	print CHECKPOINT,check_point
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	growth = 0.0086848
	if growth>split_time:
		growth=split_time-0.0001
	print file1
	interval = []
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth = interval[-1][1]
	growth_factor = -math.log(size_before_growth*int(args.N0)/33813.99)/growth
	first = True
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False
			continue
		if time<growth:
			interval1_beginning.append([time,size])
		if t*10000<CHECKPOINT:
			if time > split_time and time<check_point:
				interval1_inbetween.append([time,size])
	cmds += '-n 1 %f -n 2 0.67628 -g 2 %f -eg %f 2 0 ' %(start_size,growth_factor,growth)
	cmds += '-en %f 2 %f ' %(growth+0.000001,size_before_growth)
	for i in range(len(interval1_beginning)):
		cmds += '-en %f 1 %f ' %(interval1_beginning[i][0],interval1_beginning[i][1])
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
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
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	for i in range(len(interval1_inbetween)):
		cmds += '-en %f 1 %f ' %(interval1_inbetween[i][0],interval1_inbetween[i][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time+0.000001)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*10))
		print >>f_cmds, new_cmds

def msmc_to_macs_no_checkpoint(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.20
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	growth = 0.0086848
	if growth>split_time:
		growth=split_time-0.0001
#	print file1
	interval = []
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth = interval[-1][1]
	growth_factor = -math.log(size_before_growth*int(args.N0)/33813.99)/growth
	first = True
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False
			continue
		if time<growth:
			interval1_beginning.append([time,size])
	cmds += '-n 1 %f -n 2 0.67628 -g 2 %f -eg %f 2 0 ' %(start_size,growth_factor,growth)
	cmds += '-en %f 2 %f ' %(growth+0.000001,size_before_growth)
	for i in range(len(interval1_beginning)):
		cmds += '-en %f 1 %f ' %(interval1_beginning[i][0],interval1_beginning[i][1])
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
#		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if float(i[0])<t*10000:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	cmds += '-en %f 1 %f ' %(split_time+0.000002,share_interval[0][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time+0.000001)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*107))
		print >>f_cmds, new_cmds

def msmc_to_macs_no_checkpoint_v2(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.20
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	growth = 0.0086848
	if growth>split_time:
		growth=split_time-0.0001
#	print file1
	interval = []
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth = interval[-1][1]
	growth_factor = -math.log(size_before_growth*int(args.N0)/(33813.99*2.36/1.25))/growth
	first = True
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False
			continue
		if time<growth:
			interval1_beginning.append([time,size])
	cmds += '-n 1 %f -n 2 1.27 -g 2 %f -eg %f 2 0 ' %(start_size,growth_factor,growth)
	cmds += '-en %f 2 %f ' %(growth+0.000001,size_before_growth)
	for i in range(len(interval1_beginning)):
		cmds += '-en %f 1 %f ' %(interval1_beginning[i][0],interval1_beginning[i][1])
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
#		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint (constrain it to 100kyrs ago)
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if float(i[0])<10*10000:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-en %f 1 %f ' %(split_time-0.000001,share_interval[0][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time-0.000002)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*107))
		print >>f_cmds, new_cmds

def msmc_to_msHOT_no_checkpoint(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	args.L = 1e07
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.25
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'msHOT 4 100 -t %i -r %i %i -I 2 2 2 ' %(int(theta),int(pho),args.L)
	growth = 0.0086848
	if growth>split_time:
		growth=split_time
	print file1
	interval = []
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth = interval[-1][1]
	growth_factor = -math.log(size_before_growth*int(args.N0)/33813.99)/growth
	first = True
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False
			continue
		if time<growth:
			interval1_beginning.append([time,size])
	cmds += '-n 1 %f -n 2 0.67628 -g 2 %f -eg %f 2 0 ' %(start_size,growth_factor,growth)
	cmds += '-en %f 2 %f ' %(growth,size_before_growth)
	for i in range(len(interval1_beginning)):
		cmds += '-en %f 1 %f ' %(interval1_beginning[i][0],interval1_beginning[i][1])
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if float(i[0])<t*10000:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	cmds += '-en %f 1 %f ' %(split_time,share_interval[0][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time)

	cmds += '>sim_msHOT_%s_%i_%i.txt' %(args.pop,t,m)
	print >>f_cmds,cmds
	cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-ms-multiple.py sim_%s_%i_%i.txt 2 %i both > sim_%s_%i_%i_psmcfa.log' %(args.pop,t,m,args.L,args.pop,t,m)
	print >>f_psmcfa,cmds


def msmc_to_macs_no_checkpoint_YRI_MKK(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.2
	split_time = float(t*10000)/4/int(args.N0)/G
	check_point = float(CHECKPOINT)/4/int(args.N0)/G
	print check_point
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		start_size1 = size
		break
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		start_size2 = size
		break	
	cmds += '-n 1 %f -n 2 %f ' %(start_size1,start_size2)
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<=split_time and time>0.005:
			interval2.append([time,size])
		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if float(i[0])<t*10000:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	cmds += '-en %f 1 %f ' %(split_time+0.000002,share_interval[0][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time+0.000001)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*10))
		print >>f_cmds, new_cmds

def msmc_to_macs_no_checkpoint_GIH(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # have exponential growth 
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.14
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	growth = 0.0086848
	if growth>split_time:
		growth=split_time-0.0001
	print file1
	interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth1 = interval[-1][1]
	growth_factor1 = -math.log(size_before_growth1*int(args.N0)/45042)/growth
	interval = []
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth:
			break
		interval.append([time,size])
	size_before_growth2 = interval[-1][1]
	growth_factor2 = -math.log(size_before_growth2*int(args.N0)/33813.99)/growth
	first = True
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False
			continue
		if t*10000<CHECKPOINT:
			if time > split_time and time<check_point:
				interval1_inbetween.append([time,size])
	cmds += '-n 1 0.90084 -n 2 0.67628 -g 1 %f -g 2 %f -eg %f 1 0 -eg %f 2 0 ' %(growth_factor1,growth_factor2,growth,growth+0.000001)
	cmds += '-en %f 1 %f -en %f 2 %f ' %(growth+0.000002,size_before_growth1,growth+0.000003,size_before_growth2)
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval1.append([time,size])
		print time,size
	print interval1
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
		print time,size
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if float(i[0])<t*10000:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][1]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	cmds += '-en %f 1 %f ' %(split_time+0.000002,share_interval[0][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-eM %f 0 ' %(split_time+0.000001)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*107))
		print >>f_cmds, new_cmds

def msmc_to_macs_v3(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # no exponential growth, all psmc result
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.14
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	check_point = float(CHECKPOINT)/4/int(args.N0)/G
	print check_point
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	print file1
	interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		start_size1 = size
		break
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		start_size2 = size
		break	
	cmds += '-n 1 %f -n 2 %f ' %(start_size1,start_size2)
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
	print interval1
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if time<=split_time and time>0.005:
			interval2.append([time,size])
	print interval2
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
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
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
#		cmds += '-ema 0.0025 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 15kyr
#		cmds += '-ema 0.004166667 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 25kyr
#		cmds += '-ema 0.008333333 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 50kyr
		cmds += '-eM %f 0 ' %(split_time+0.000002)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*10))
		print >>f_cmds, new_cmds
#		msmc_cmds += 'sim_macs_MKK_CEU_%i_%i_%i_sites_2.txt ' %(t,m,i)
#	print >>f_psmc,msmc_cmds

def msmc_to_macs_v4(file1,file2,file3,t,m,f_cmds,f_psmcfa,migration_end,CHECKPOINT): # no exponential growth, all psmc result, simulate a single bottleneck for CEU, don't match the curve for CEU after the split, test
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.14
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	check_point = float(CHECKPOINT)/4/int(args.N0)/G
	print check_point
	cmds = 'macs 4 %i -s SEED -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
	print file1
	interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		start_size1 = size
		break
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		start_size2 = size
		break	
	cmds += '-n 1 %f -n 2 %f ' %(start_size1,start_size2)
	interval1 = []
	interval2 = []
	share_interval = []
#	print file1
	interval1_beginning = []
	interval1_inbetween = []   # time interval after split time but before checkpoint 
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
	print interval1
	small_index = 0
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
		if time<=split_time and time>0.005:
			interval2.append([time,size])
			if len(interval2)>1:
				if interval2[-1][1]<interval2[small_index][1]:
					small_index=len(interval2)-1
	print interval2
	print 'small',interval2[small_index]
	for i in BedIterator(file3):   # only before the split, when two pop have shared history, also must before the checkpoint
		time = abs(float(i[0]))/G/4/int(args.N0)
		size = float(i[1])*10000/int(args.N0)
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
#	for i in range(len(interval2)):
#		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-en %f 2 %f ' %(interval2[small_index][0],interval2[small_index][1])
	cmds += '-ej %s 2 1 ' %(split_time)
#	cmds += '-en %s 1 %f ' %(split_time+0.000001,start_size1)          ## change it! ##
	cmds += '-en %s 1 %f ' %(split_time+0.000001,share_interval[0][1])
#	for i in range(len(interval1_inbetween)):
#		cmds += '-en %f 1 %f ' %(interval1_inbetween[i][0],interval1_inbetween[i][1])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		migration_end = float(migration_end)*1000/G/4/args.N0
		cmds += '-ema %f 2 0 %f %f 0 ' %(migration_end,migration_rate,migration_rate)   # migration starting from 15kyr
#		cmds += '-ema 0.0025 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 15kyr
#		cmds += '-ema 0.004166667 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 25kyr
#		cmds += '-ema 0.008333333 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 50kyr
		cmds += '-eM %f 0 ' %(split_time+0.000002)

	for i in range(args.nrep):
		new_cmds="%s-T 2>sim_macs_%s_%i_%i_%i_trees.txt 1>sim_macs_%s_%i_%i_%i_sites.txt" %(cmds,args.pop,t,m,i,args.pop,t,m,i)
		new_cmds = new_cmds.replace("SEED",str(1420845033+i*10))
		print >>f_cmds, new_cmds
#		msmc_cmds += 'sim_macs_MKK_CEU_%i_%i_%i_sites_2.txt ' %(t,m,i)
#	print >>f_psmc,msmc_cmds
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--nrep", default=300,dest='nrep',help="length of fragments")
	parser.add_argument("--pop", default=100,dest='pop',help="population name")
	parser.add_argument("--migration_end", default=100,dest='migration_end',help="migration end in kyr")
	parser.add_argument("--checkpoint", default=100,dest='checkpoint',help="checkpoint in yr")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	parser.add_argument("--version", dest='version',help="psmc input or not")
	args = parser.parse_args()
	f_psmcfa = open('macs_to_psmc_grid_search.cmds','a')
	args.L=int(args.L)
	args.nrep = int(args.nrep)
	f_cmds = open('macs_grid_search.cmds','a')
#	f_cmds_mshot = open('msHOT_grid_search_changer.cmds','a')
#	f_psmcfa_cmds = open('psmcfa.cmds','a')
	for i in range(6,16,2):
#		f_cmds = open('macs_grid_search_%s.cmds' %(i),'w')
		for j in [10]:
			print i,j
			if args.file3 is not None:
				if args.version is not None:
					msmc_to_macs_v3(args.file1,args.file2,args.file3,i,j,f_cmds,f_psmcfa,int(args.migration_end),int(args.checkpoint))
				else:
#					msmc_to_macs_no_checkpoint(args.file1,args.file2,args.file3,i,j,f_cmds,f_psmcfa,int(args.migration_end),int(args.checkpoint))
					msmc_to_macs_no_checkpoint_v2(args.file1,args.file2,args.file3,i,j,f_cmds,f_psmcfa,int(args.migration_end),int(args.checkpoint))
#					msmc_to_msHOT_no_checkpoint(args.file1,args.file2,args.file3,i,j,f_cmds_mshot,f_psmcfa_cmds,int(args.migration_end),int(args.checkpoint))
				print >>f_psmcfa,"python /home/songsy/script/pop-gen/make-psmc-from-macs.py --file sim_macs_%s_%i_%i --nrep %s --nchr 100" %(args.pop,i,j,args.nrep)

	
