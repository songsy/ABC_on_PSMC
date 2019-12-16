#!/usr/bin/env python
# python make-psmc-from-macs.py
# Shiya Song
# 11 Jan 2015

import sys
import argparse
import random
import os
def code_to_string(code):
	string = ''
	for i in range(len(code)):
		if code[i]=='0':
			string+='A'
		else:
			string+='T'
	return string

def have_het(string):
	a=list(string)
	zero = a.count('0')
	one = a.count('1')
	if zero==0 or one ==0:
		return False
	else:
		return True

def macs_to_record(f):
	pos = [[],[],[],[]]
	prev_pos = [0,0,0,0]
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			for j in range(int(args.nsample)):
				for hap in range(int(args.nhap)):
#					print j,hap+int(args.nhap)
					substring = col[4][j]+col[4][hap+int(args.nhap)]
					if have_het(substring) is True:
						if int(int(args.L)*float(col[2]))<=prev_pos[index]:
							pos[index].append(prev_pos[index]+1)
							prev_pos[index] +=1
						else:
							pos[index].append(int(int(args.L)*float(col[2])))
							prev_pos[index] = int(int(args.L)*float(col[2]))
					index +=1
	return pos

def macs_to_record_v2(f):
	pos = []
	prev_pos = 0
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			substring = col[4][0]+col[4][2]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos:
					pos.append(prev_pos+1)
					prev_pos +=1
				else:
					pos.append(int(int(args.L)*float(col[2])))
					prev_pos = int(int(args.L)*float(col[2]))
	return pos

def macs_to_record_v3(f):    # first pseudo-diploid, then diploi
	pos = [[],[],[],[],[],[]]
	prev_pos = [0,0,0,0,0,0]
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			for j in range(int(args.nsample)):
				for hap in range(int(args.nhap)):
#					print j,hap+int(args.nhap)
					substring = col[4][j]+col[4][hap+int(args.nhap)]
					if have_het(substring) is True:
						if int(int(args.L)*float(col[2]))<=prev_pos[index]:
							pos[index].append(prev_pos[index]+1)
							prev_pos[index] +=1
						else:
							pos[index].append(int(int(args.L)*float(col[2])))
							prev_pos[index] = int(int(args.L)*float(col[2]))
					index +=1
			substring = col[4][0]+col[4][1]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[4]:
					pos[4].append(prev_pos[4]+1)
					prev_pos[4] +=1
				else:
					pos[4].append(int(int(args.L)*float(col[2])))
					prev_pos[4] = int(int(args.L)*float(col[2]))
			substring = col[4][2]+col[4][3]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[5]:
					pos[5].append(prev_pos[5]+1)
					prev_pos[5] +=1
				else:
					pos[5].append(int(int(args.L)*float(col[2])))
					prev_pos[5] = int(int(args.L)*float(col[2]))
	return pos

def macs_to_record_with_msmc(f):
	map_code={'1':'A','0':'T'}
	pos = [[],[],[]]
	prev_pos = [0,0,0]
	STRING = []
	POS = []
	sample = args.nsample
	a=random.randint(0,sample/2-1)
	b=random.randint(0,sample/2-1)
	hap=args.nhap
	print a,b
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			substring = col[4][a*hap]+col[4][a*hap+1]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[0]:
					pos[0].append(prev_pos[0]+1)
					prev_pos[0] +=1
				else:
					pos[0].append(int(int(args.L)*float(col[2])))
					prev_pos[0] = int(int(args.L)*float(col[2]))
			substring = col[4][sample+b*hap]+col[4][sample+b*hap+1]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[1]:
					pos[1].append(prev_pos[1]+1)
					prev_pos[1] +=1
				else:
					pos[1].append(int(int(args.L)*float(col[2])))
					prev_pos[1] = int(int(args.L)*float(col[2]))
			substring = col[4][a*hap]+col[4][sample+b*hap]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[2]:
					pos[2].append(prev_pos[2]+1)
					prev_pos[2] +=1
				else:
					pos[2].append(int(int(args.L)*float(col[2])))
					prev_pos[2] = int(int(args.L)*float(col[2]))
			wholestring=col[4][a*hap]+col[4][a*hap+1]+col[4][sample+b*hap]+col[4][sample+b*hap+1]
			if have_het(wholestring) is True:
				if len(POS)==0:
					POS.append(int(int(args.L)*float(col[2])))
				elif int(int(args.L)*float(col[2]))<=POS[-1]:
					POS.append(POS[-1]+1)
				else:
					POS.append(int(int(args.L)*float(col[2])))
				wholestring=wholestring.replace('1','A')
				wholestring=wholestring.replace('0','T')
				STRING.append(wholestring)
	return pos,POS,STRING

def write_psmcfa(pos,f_out,index):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	for comparison in range(len(pos)):
		outFile = f_out[comparison]
		outFile.write('>%s\n' % (index))
		windowStart = 0
		windowEnd = windowStart + winLen - 1
		windows = []
		while windowStart < int(args.L):
			windows.append(homChar)
			windowStart = windowEnd + 1
			windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
		for p in pos[comparison]:
			wID = int(p/winLen)
			try:
				windows[wID] = hetChar
			except IndexError:
				print "Error:",wID,p,winLen,len(windows)
				windows[wID-1] = hetChar
		numHet = 0
		numHom = 0
		for i in range(len(windows)):
			outFile.write('%s' % (windows[i]))
			if windows[i] == hetChar:
				numHet += 1
			else:
				numHom += 1
			if (i+1) % 50 == 0:
				outFile.write('\n')
		outFile.write('\n')
		print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i' % (index,comparison,len(pos[comparison]),len(windows),numHom,numHet)

def write_psmcfa_v2(pos,f_out,index):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	if True:
		outFile = f_out
		outFile.write('>%s\n' % (index))
		windowStart = 0
		windowEnd = windowStart + winLen - 1
		windows = []
		while windowStart < int(args.L):
			windows.append(homChar)
			windowStart = windowEnd + 1
			windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
		for p in pos:
			wID = int(p/winLen)
			try:
				windows[wID] = hetChar
			except IndexError:
				print "Error:",wID,p,winLen,len(windows)
				windows[wID-1] = hetChar
		numHet = 0
		numHom = 0
		for i in range(len(windows)):
			outFile.write('%s' % (windows[i]))
			if windows[i] == hetChar:
				numHet += 1
			else:
				numHom += 1
			if (i+1) % 50 == 0:
				outFile.write('\n')
		outFile.write('\n')
		print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i' % (index,1,len(pos),len(windows),numHom,numHet)

def read_macs(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ." 
			f_out = []
			index +=1
			for j in range(1,5):
				f_out.append(open(f_out_prefix+str(j),"w"))
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))

def read_macs_v2(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ.0" 
			index +=1
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			f_out = open(f_out_prefix,"w")
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))

def read_macs_v3(macs_file):  #diploid+4*pseudo-diploid
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ." 
			f_out = []
			index +=1
			for j in range(6):
				f_out.append(open(f_out_prefix+str(j),"w"))
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v3(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v3(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))

def write_msmcfa_multiple(POS,STRING,f_msmc_out,index):
	for i in range(len(POS)):
		if i==0:
			print >>f_msmc_out,'%s\t%i\t%i\t%s' %(index,POS[i],POS[i],STRING[i])
		else:
			print >>f_msmc_out,'%s\t%i\t%i\t%s' %(index,POS[i],POS[i]-POS[i-1],STRING[i])

def read_macs_with_msmc(macs_file):  #diploid+4*pseudo-diploid
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ." 
			f_out = []
			index +=1
			for j in range(3):
				f_out.append(open(f_out_prefix+str(j),"w"))
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,POS,STRING = macs_to_record_with_msmc(f)
			print 'len:',len(pos),len(POS),len(STRING)
			write_psmcfa(pos,f_out,i % int(args.nchr))
			f_msmc_out = open(macs_file+"_"+str(i)+"_sites_0.txt" ,'w')
			write_msmcfa_multiple(POS,STRING,f_msmc_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,POS,STRING = macs_to_record_with_msmc(f)
			print 'len:',len(pos),len(POS),len(STRING)
			write_psmcfa(pos,f_out,i % int(args.nchr))
			f_msmc_out = open(macs_file+"_"+str(i)+"_sites_0.txt" ,'w')
			write_msmcfa_multiple(POS,STRING,f_msmc_out,i % int(args.nchr))

if __name__=="__main__":

	parser = argparse.ArgumentParser(description='macs to psmcfa file')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--nrep", default=100,dest='nrep',help="number of repetition")
	parser.add_argument("--nchr", default=100,dest='nchr',help="number of chromosome")
	parser.add_argument("--nhap", default=2,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", type=int,default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--version",default=None,dest='version')
	parser.add_argument("--mode",default=None,dest='mode')

	args = parser.parse_args()
	if args.mode=='all':
		read_macs_v3(args.file)
	elif args.mode=='msmc':
		read_macs_with_msmc(args.file)
	else:
		read_macs_v2(args.file)
#	os.popen('rm %s_*sites.txt' %(args.file))
#	os.popen('rm %s_*trees.txt' %(args.file))
	