#90000

#python /home/songsy/script/pop-gen/msmc_to_msHOT_GIH_CEU.py --file1 ../NA20847/NA20847_gvcf_v2_msmc.final.txt --file2 ../CEU/NA12878_1KG_msmc.final.txt --file3 ../African_CEU/psmcfa/CEU-GIH-psmc.0.txt
#python /home/songsy/script/pop-gen/msmc_to_msHOT_GIH_CEU.py --file1 ../NA20847/NA20847_gvcf_v2_msmc.final.txt --file2 ../CEU/NA12878_1KG_msmc.final.txt --file3 ../African_CEU/psmcfa/CEU-GIH-true-fosmid-psmc.0.txt

import argparse
from NGS_utils import *
import math

G=30
def msmc_to_msHOT(file1,file2,t,m,f_cmds,f_psmcfa,f_psmc):
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.13
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'msHOT 4 %i -t %i -r %i %i -I 2 2 2 ' %(args.nrep,int(theta),int(pho),args.L)
	growth = 0.0086848
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
	growth_factor = -math.log(size_before_growth*int(args.N0)/33813.99)/0.0086848
	first = True
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if first == True:
			start_size = size
			first = False		
	cmds += '-n 1 %f -n 2 0.67628 -g 2 %f -eg 0.0086848 2 0 ' %(start_size,growth_factor)
	cmds += '-en 0.0086848 2 %f ' %(size_before_growth)
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
		if time > split_time:
			if len(share_interval)==0:
				share_interval.append([time,size])
			else:
				if size==share_interval[-1][0]:
					continue
				share_interval.append([time,size])
#		print time,size
#	print file2
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time > growth and time<=split_time:
			interval2.append([time,size])
#		print time,size
	for i in range(len(interval1)):
		cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
		cmds += '-ema 0.008333333 2 0 %f %f 0 ' %(migration_rate,migration_rate)
		cmds += '-eM %f 0 ' %(split_time)
	cmds += '>sim_GIH_CEU_%i_%i.txt' %(t,m)
	print >>f_cmds,cmds
	print cmds
	cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-ms-multiple.py sim_GIH_CEU_%i_%i.txt 2 %i both > sim_GIH_CEU_%i_%i_psmcfa.log' %(t,m,args.L,t,m)
	print >>f_psmcfa,cmds
	print cmds
	for i in range(0,3):
		cmds = 'psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o sim_GIH_CEU_%i_%i.%i.psmc sim_GIH_CEU_%i_%i.hetFQ.%i' %(t,m,i,t,m,i)
		print >>f_psmc,cmds
		print cmds

def msmc_to_msHOT_v2(file1,file2,file3,t,m,f_cmds,f_psmcfa,f_psmc):
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.13
	split_time = float(t*10000)/4/int(args.N0)/G
	check_point = float(90000)/4/int(args.N0)/G
	print check_point
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'msHOT 4 %i -t %i -r %i %i -I 2 2 2 ' %(args.nrep,int(theta),int(pho),args.L)
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
	growth_factor1 = -math.log(size_before_growth1*int(args.N0)/33813.99)/growth
	growth_factor2 = -math.log(size_before_growth2*int(args.N0)/45042.13)/growth
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
		if t*10000<90000:
			if time > split_time and time<check_point:
				interval1_inbetween.append([time,size])
#	cmds += '-n 1 1.69032 -n 2 0.67628 -g 1 %f -g 2 %f -eg %f 2 0 -eg %f 1 0 ' %(growth_factor1,growth_factor2,growth,growth)
	cmds += '-n 1 0.90084 -n 2 0.67628 -g 1 %f -g 2 %f -eg %f 2 0 -eg %f 1 0 ' %(growth_factor1,growth_factor2,growth,growth)
	cmds += '-en %f 1 %f ' %(growth,size_before_growth1)
	cmds += '-en %f 2 %f ' %(growth,size_before_growth2)
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
		if float(i[0])<90000:
			continue
		if t*10000> 90000:
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
#		cmds += '-ema 0.0025 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-ema 0.004166667 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 25kyr
		cmds += '-eM %f 0 ' %(split_time)
	cmds += '>sim_GIH_CEU_%i_%i.txt' %(t,m)
	print >>f_cmds,cmds
	print cmds
	cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-ms-multiple.py sim_GIH_CEU_%i_%i.txt 2 %i both > sim_GIH_CEU_%i_%i_psmcfa.log' %(t,m,args.L,t,m)
	print >>f_psmcfa,cmds
	print cmds
	for i in range(0,3):
		cmds = 'psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o sim_GIH_CEU_%i_%i.%i.psmc sim_GIH_CEU_%i_%i.hetFQ.%i' %(t,m,i,t,m,i)
		print >>f_psmc,cmds
		print cmds
		
def msmc_to_msHOT_v3(file1,file2,file3,t,m,f_cmds,f_psmcfa,f_psmc):   # simulate  no exponential growth
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.13
	split_time = float(t*10000)/4/int(args.N0)/G
	CHECKPOINT=90000
	check_point = float(CHECKPOINT)/4/int(args.N0)/G
	print check_point
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'msHOT 4 %i -t %i -r %i %i -I 2 2 2 ' %(args.nrep,int(theta),int(pho),args.L)
	print file1
	interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if size<0.0000001:
			continue
		start_size1 = size
		break
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if size<0.0000001:
			continue
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
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
		if t*10000<CHECKPOINT:
			if time > split_time and time<check_point:
				interval1_inbetween.append([time,size])
	print interval1
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
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
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	for i in range(len(interval1_inbetween)):
		cmds += '-en %f 1 %f ' %(interval1_inbetween[i][0],interval1_inbetween[i][1])
#		cmds += '-en %f 1 0.3 ' %(interval1_inbetween[i][0])
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
	if m!=0:
#		cmds += '-ema 0.0025 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 15kyr
		cmds += '-ema 0.004166667 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 25kyr
#		cmds += '-ema 0.008333333 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 50kyr
		cmds += '-eM %f 0 ' %(split_time)
	cmds += '>sim_GIH_CEU_%i_%i.txt' %(t,m)
	print >>f_cmds,cmds
	print cmds
	cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-ms-multiple.py sim_GIH_CEU_%i_%i.txt 2 %i both > sim_GIH_CEU_%i_%i_psmcfa.log' %(t,m,args.L,t,m)
	print >>f_psmcfa,cmds
	print cmds
	for i in range(0,3):
		cmds = 'psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o sim_GIH_CEU_%i_%i.%i.psmc sim_GIH_CEU_%i_%i.hetFQ.%i' %(t,m,i,t,m,i)
		print >>f_psmc,cmds
		print cmds
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--nrep", default=100,dest='nrep',help="length of fragments")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	args = parser.parse_args()
	f_cmds = open('msHOT_grid_search.cmds','a')
	f_psmcfa = open('psmcfa_grid_search.cmds','a')
	f_psmc = open('psmc_grid_search.cmds','a')
	'''
	if args.file3 is not None:
		msmc_to_msHOT_v2(args.file1,args.file2,args.file3,10,0,f_cmds,f_psmcfa,f_psmc)
	else:
		msmc_to_msHOT(args.file1,args.file2,10,0,f_cmds,f_psmcfa,f_psmc)
	'''
	for i in range(4,11):
		for j in range(0,22,2):
			print i,j
			if args.file3 is not None:
				msmc_to_msHOT_v2(args.file1,args.file2,args.file3,i,j,f_cmds,f_psmcfa,f_psmc)
			else:
				msmc_to_msHOT(args.file1,args.file2,i,j,f_cmds,f_psmcfa,f_psmc)


	