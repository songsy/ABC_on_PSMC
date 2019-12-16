#137465.484312

#python /home/songsy/script/pop-gen/msmc_to_macs_YRI_MKK.py --file1 ../YRI/NA19240_gvcf_msmc.final.txt --file2 ../MKK/NA21302_gvcf_msmc.final.txt
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
	cmds += '>sim_MKK_CEU_%i_%i.txt' %(t,m)
	print >>f_cmds,cmds
	print cmds
	cmds = 'python /home/songsy/script/pop-gen/make-psmc-from-ms-multiple.py sim_MKK_CEU_%i_%i.txt 2 %i both > sim_MKK_CEU_%i_%i_psmcfa.log' %(t,m,args.L,t,m)
	print >>f_psmcfa,cmds
	print cmds
	for i in range(0,3):
		cmds = 'psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o sim_MKK_CEU_%i_%i.%i.psmc sim_MKK_CEU_%i_%i.hetFQ.%i' %(t,m,i,t,m,i)
		print >>f_psmc,cmds
		print cmds

def msmc_to_msHOT_v2(file1,file2,t,m,f_cmds,f_psmcfa,f_psmc):
	theta = 4*int(args.N0)*float(args.mu)*int(args.L)
	pho = theta*0.25
	split_time = float(t*10000)/4/int(args.N0)/G
	migration_rate = m/1e05*4*int(args.N0)
	print split_time,migration_rate
	cmds = 'macs 4 %i -t %f -r %f -I 2 2 2 ' %(args.L,theta/int(args.L),pho/int(args.L))
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
	interval1 = []
	interval2 = []
	share_interval = []
	for i in BedIterator(file1):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<=split_time and time>0.005:
			interval1.append([time,size])
	for i in BedIterator(file2):
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<=split_time and time>0.005:
			interval2.append([time,size])
	for i in BedIterator(file1):   # only before the split, when two pop have shared history, also must before the checkpoint
		if i[0][0]=='t':
			continue
		time = abs(float(i[1]))/theta*int(args.L)
		size = 1/float(i[3])/theta*2*int(args.L)
		if time<split_time or time>0.6:
			continue
		if len(share_interval)==0:
			share_interval.append([time,size])
		else:
			if size==share_interval[-1][0]:
				continue
			share_interval.append([time,size])
	for i in range(len(interval2)):
		if i<len(interval1):
			cmds += '-en %f 1 %f ' %(interval1[i][0],interval1[i][1])
		cmds += '-en %f 2 %f ' %(interval2[i][0],interval2[i][1])
	cmds += '-ej %s 2 1 ' %(split_time)
	print share_interval
	for i in range(len(share_interval)):
		cmds += '-en %f 1 %f ' %(share_interval[i][0],share_interval[i][1])
#	if m!=0:
#		cmds += '-ema 0.0025 2 0 %f %f 0 ' %(migration_rate,migration_rate)   # migration starting from 25kyr
#		cmds += '-eM %f 0 ' %(split_time-0.0001)
	msmc_cmds = 'msmc --fixedRecombination --skipAmbiguous -P 0,0,1,1 -o sim_macs_MKK_CEU_%i_%i_bundle_2 ' %(t,m)
	for i in range(args.nrep):
		print >>f_cmds,"%s-T 2>sim_macs_MKK_CEU_%i_%i_%i_trees.txt 1>sim_macs_MKK_CEU_%i_%i_%i_sites.txt" %(cmds,t,m,i,t,m,i)
		print >>f_psmcfa,"python /home/songsy/script/pop-gen/make-msmc-from-macs.py sim_macs_MKK_CEU_%i_%i_%i_sites.txt %i 2 2" %(t,m,i,args.L)
		msmc_cmds += 'sim_macs_MKK_CEU_%i_%i_%i_sites_2.txt ' %(t,m,i)
	print >>f_psmc,msmc_cmds

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--N0", default=5e04,dest='N0',help="N0")
	parser.add_argument("--mu", default=1.25e-08,dest='mu',help="mutation rate")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--nrep", default=100,dest='nrep',help="length of fragments")
	parser.add_argument("--file1", dest='file1',help="African")
	parser.add_argument("--file2", dest='file2',help="CEU")
	parser.add_argument("--file3", dest='file3',help="African-vs-CEU")
	parser.add_argument("--i",dest='i',help='split time')
	parser.add_argument("--j",dest='j',help='migration time')
	args = parser.parse_args()
	f_cmds = open('macs_%s_%s.cmds' %(args.i,args.j),'a')
	f_psmcfa = open('macs_to_msmc_%s_%s.cmds' %(args.i,args.j),'a')
	f_psmc = open('msmc_grid_search.cmds','a')

	i=int(args.i)
	j=int(args.j)
	msmc_to_msHOT_v2(args.file1,args.file2,i,j,f_cmds,f_psmcfa,f_psmc)

	