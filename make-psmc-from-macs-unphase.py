#!/usr/bin/env python
# python make-psmc-from-macs.py
# Shiya Song
# 11 Jan 2015

import sys
import argparse
import random

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

def macs_to_record_single(f):
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

def macs_to_record_unphased(f):
	pos = [[],[],[],[],[]]	 # normal, random,all_het, skip_mask, skip
	prev_pos = [0,0,0,0,0]
	mask_pos = []
	mask_pos_v2 = []
	unknown_status = False
	unphase = 0
	tot = 0
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			index = 0
			tot +=1
			col = line.split()
			substring = col[4][0]+col[4][2]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[index]:
					pos[index].append(prev_pos[index]+1)
					prev_pos[index] +=1
				else:
					pos[index].append(int(int(args.L)*float(col[2])))
					prev_pos[index] = int(int(args.L)*float(col[2]))
			s=random.uniform(0, 1)
			if s<=args.freq:
				unphase +=1
				index +=1
				# random mode 
				a = random.randint(0,1)
				b = random.randint(0,1)
				substring = col[4][0+a]+col[4][2+b]
				if have_het(substring) is True:
					if int(int(args.L)*float(col[2]))<=prev_pos[index]:
						pos[index].append(prev_pos[index]+1)
						prev_pos[index] +=1
					else:
						pos[index].append(int(int(args.L)*float(col[2])))
						prev_pos[index] = int(int(args.L)*float(col[2]))
				index +=1  #all het mode
				het = False
				for i in range(2):
					for j in range(2):
						substring = col[4][0+i]+col[4][2+j]
						if have_het(substring) is True:
							het = True
				if het:
					if int(int(args.L)*float(col[2]))<=prev_pos[index]:
						pos[index].append(prev_pos[index]+1)
						prev_pos[index] +=1
					else:
						pos[index].append(int(int(args.L)*float(col[2])))
						prev_pos[index] = int(int(args.L)*float(col[2]))
				index +=1 # skip-mask mode
				if int(int(args.L)*float(col[2]))<=prev_pos[index]:
					pos[index].append(prev_pos[index]+1)
					mask_pos.append(1)
					prev_pos[index] +=1
				else:
					pos[index].append(int(int(args.L)*float(col[2])))
					prev_pos[index] = int(int(args.L)*float(col[2]))
					mask_pos.append(1)
				index +=1  # skip no mask mode
				if int(int(args.L)*float(col[2]))<=prev_pos[index]:
					pos[index].append(prev_pos[index]+1)
					mask_pos_v2.append(1)
					prev_pos[index] +=1
				else:
					pos[index].append(int(int(args.L)*float(col[2])))
					prev_pos[index] = int(int(args.L)*float(col[2]))
					mask_pos_v2.append(1)
			else:
				substring = col[4][0]+col[4][2]
				if have_het(substring) is True:
					for index in range(1,5):
						if int(int(args.L)*float(col[2]))<=prev_pos[index]:
							pos[index].append(prev_pos[index]+1)
							prev_pos[index] +=1
							if index==3:
								mask_pos.append(0)
								mask_pos_v2.append(0)
						else:
							pos[index].append(int(int(args.L)*float(col[2])))
							prev_pos[index] = int(int(args.L)*float(col[2]))
							if index==3:
								mask_pos.append(0)
								mask_pos_v2.append(0)
	assert len(mask_pos)==len(pos[3])
	print 'Length:',len(pos[0]),len(pos[1]),len(pos[2]),len(pos[3]),len(pos[4]),unphase,tot
	return pos,mask_pos,mask_pos_v2

def macs_to_record_all(f):
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

def macs_to_record_error(f):
	map_code={'1':'A','0':'T'}
	rate = 1.0/args.length
	switch_state = [False,False]
	next_switch = [0,0]
	pos = [[],[],[],[]]
	prev_pos = [0,0,0,0]
	first_pos = None
	STRING = []
	POS = []
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			if first_pos is None:
				first_pos = int(int(args.L)*float(col[2]))
				for j in range(2):
					next_switch[j]=int(random.expovariate(rate))
			substring = ["","","",""]  # 11,12,21,22
			for j in range(2):  # number of individual
				if int(int(args.L)*float(col[2]))-first_pos>=next_switch[j]:
					switch_state[j]= not switch_state[j]
					next_switch[j]+=int(random.expovariate(rate))
				if j ==0:
					if switch_state[j]:
						substring[0]+=col[4][j*2+1]
						substring[1]+=col[4][j*2+1]
						substring[2]+=col[4][j*2]
						substring[3]+=col[4][j*2]
					else:
						substring[0]+=col[4][j*2]
						substring[1]+=col[4][j*2]
						substring[2]+=col[4][j*2+1]
						substring[3]+=col[4][j*2+1]
				elif j ==1:
					if switch_state[j]:
						substring[0]+=col[4][j*2+1]
						substring[1]+=col[4][j*2]
						substring[2]+=col[4][j*2+1]
						substring[3]+=col[4][j*2]
					else:
						substring[0]+=col[4][j*2]
						substring[1]+=col[4][j*2+1]
						substring[2]+=col[4][j*2]
						substring[3]+=col[4][j*2+1]
			wholestring=substring[0][0]+substring[3][0]+substring[0][1]+substring[3][1]
			wholestring=wholestring.replace('1','A')
			wholestring=wholestring.replace('0','T')
			STRING.append(wholestring)
			if len(POS)==0:
				POS.append(int(int(args.L)*float(col[2])))
			elif int(int(args.L)*float(col[2]))<=POS[-1]:
				POS.append(POS[-1]+1)
			else:
				POS.append(int(int(args.L)*float(col[2])))
			for index in range(4):
				if have_het(substring[index]) is True:
					if int(int(args.L)*float(col[2]))<=prev_pos[index]:
						pos[index].append(prev_pos[index]+1)
						prev_pos[index] +=1
					else:
						pos[index].append(int(int(args.L)*float(col[2])))
						prev_pos[index] = int(int(args.L)*float(col[2]))
	return pos,POS,STRING

def write_psmcfa_unphased(pos,f_out,index,mask_pos,mask_pos_v2):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	for comparison in range(len(pos)-2):
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
	# last comparison: with skip-mask mode
	outFile = f_out[-2]
	outFile.write('>%s\n' % (index))
	windowStart = 0
	windowEnd = windowStart + winLen - 1
	windows = []
	while windowStart < int(args.L):
		windows.append(homChar)
		windowStart = windowEnd + 1
		windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
	pID = 0
	start = True
	while pID < len(pos[-2]):
		if mask_pos[pID]==1:
			if start:
				start_ID = pID
				end_ID = pID
				start = False
			end_ID = pID
		else:
			if start == False:
				start = True
#				print start_ID,end_ID,pos[-1][start_ID],pos[-1][end_ID]
				if start_ID!=0:
					s_wID=int(pos[-2][start_ID-1]/winLen)
#					print 'start',s_wID
					leftover = pos[-2][start_ID-1]%winLen
#					if leftover>90:
					if leftover>10:
						try:
							windows[s_wID] = hetChar
						except IndexError:
							print "Error:",s_wID,pID,winLen,len(windows)
							windows[s_wID-1] = hetChar
					else:
#						print 'Fail',s_wID
						windows[s_wID] = failChar
				else:
					s_wID=int(pos[-2][start_ID]/winLen)
				if end_ID != len(pos[-2])-1:
					e_wID=int(pos[-2][end_ID+1]/winLen)
					leftover = winLen-pos[-2][start_ID-1]%winLen
#					print 'end',e_wID
#					if leftover>90:
					if leftover>10:
						try:
							windows[e_wID] = hetChar
						except IndexError:
							print "Error:",e_wID,pID,winLen,len(windows)
							windows[e_wID-1] = hetChar
					else:
#						print 'Fail',e_wID
						windows[e_wID] = failChar
				else:
					e_wID=int(pos[-2][end_ID]/winLen)
				for i in range(s_wID,e_wID):
					if windows[i]==homChar:
						windows[i]=failChar
			else:
				wID = int(pos[-2][pID]/winLen)
				try:
					windows[wID] = hetChar
				except IndexError:
					print "Error:",wID,p,winLen,len(windows)
					windows[wID-1] = hetChar
		pID+=1
	if mask_pos[-2]==1:
		if start_ID!=0:
			wID=int(pos[-2][start_ID-1]/winLen)
			leftover = pos[-2][start_ID-1]%winLen
#			if leftover>90:
			if leftover>10:
				try:
					windows[wID] = hetChar
				except IndexError:
					print "Error:",wID,pID,winLen,len(windows)
					windows[wID-1] = hetChar
		if end_ID != len(pos[-2])-1:
			wID=int(pos[-2][end_ID+1]/winLen)
			leftover = winLen-pos[-2][start_ID-1]%winLen
#			if leftover>90:
			if leftover>10:
				try:
					windows[wID] = hetChar
				except IndexError:
					print "Error:",wID,pID,winLen,len(windows)
					windows[wID-1] = hetChar
		for i in range(start_ID,end_ID+1):
			wID=int(pos[-2][i]/winLen)
			if windows[wID]==homChar:
				windows[wID]=failChar
	numHet = 0
	numHom = 0
	numN = 0
	for i in range(len(windows)):
		outFile.write('%s' % (windows[i]))
		if windows[i] == hetChar:
			numHet += 1
		elif windows[i] == failChar:
			numN += 1
		else:
			numHom +=1
		if (i+1) % 50 == 0:
			outFile.write('\n')
	outFile.write('\n')
	print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i NWind %i' % (index,3,len(pos[3]),len(windows),numHom,numHet,numN)

	outFile = f_out[-1]
	outFile.write('>%s\n' % (index))
	windowStart = 0
	windowEnd = windowStart + winLen - 1
	windows = []
	while windowStart < int(args.L):
		windows.append(homChar)
		windowStart = windowEnd + 1
		windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows

	for pID in range(len(pos[-1])):
			wID = int(pos[-1][pID]/winLen)
			if pID == 0:
				old_wID = wID
				uncalled = 0
			if wID !=old_wID:
				if uncalled>=10:
					windows[old_wID] = failChar
				uncalled = 0
			if mask_pos_v2[pID]==1:
				windows[wID] = homChar
				uncalled +=1
			else:
				windows[wID] = hetChar
	numHet = 0
	numHom = 0
	numN = 0
	for i in range(len(windows)):
		outFile.write('%s' % (windows[i]))
		if windows[i] == hetChar:
			numHet += 1
		elif windows[i] == failChar:
			numN += 1
		else:
			numHom +=1
		if (i+1) % 50 == 0:
			outFile.write('\n')
	outFile.write('\n')
	print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i NWind %i' % (index,3,len(pos[3]),len(windows),numHom,numHet,numN)


def write_psmcfa_multiple(pos,f_out,index):
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

def write_msmcfa_multiple(POS,STRING,f_msmc_out,index):
	for i in range(len(POS)):
		if i==0:
			print >>f_msmc_out,'%s\t%i\t%i\t%s' %(index,POS[i],POS[i],STRING[i])
		else:
			print >>f_msmc_out,'%s\t%i\t%i\t%s' %(index,POS[i],POS[i]-POS[i-1],STRING[i])

def write_psmcfa_single(pos,f_out,index):
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

def read_macs_unphase(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_%.3f_leftover10_hetFQ." %(args.freq) 
			f_out = []
			index +=1
			for j in range(5):	# org, random, all_het,skip-mask, skip-unmask
				f_out.append(open(f_out_prefix+str(j),"w"))		
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,mask_pos,mask_pos_v2 = macs_to_record_unphased(f)
			print 'len:',len(pos)
			write_psmcfa_unphased(pos,f_out,i % int(args.nchr),mask_pos,mask_pos_v2)
		else:
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,mask_pos,mask_pos_v2 = macs_to_record_unphased(f)
			print 'len:',len(pos)
			write_psmcfa_unphased(pos,f_out,i % int(args.nchr),mask_pos,mask_pos_v2)

def read_macs_error(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			s=int(args.length)/1000
			print s
			f_out_prefix = macs_file+"_"+str(index)+"_error_%ik_hetFQ." %(int(args.length)/1000)
			f_out = []
			index +=1
			for j in range(4):	# 4*pseudo-diploid with switch errors
				f_out.append(open(f_out_prefix+str(j),"w"))		
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_error(f)
			print 'len:',len(pos)
			write_psmcfa_multiple(pos,f_out,i % int(args.nchr))
		else:
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_error(f)
			print 'len:',len(pos)
			write_psmcfa_multiple(pos,f_out,i % int(args.nchr))

def read_macs_error_with_msmc(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			s=int(args.length)/1000
			print s
			f_out_prefix = macs_file+"_"+str(index)+"_error_%ik_hetFQ." %(int(args.length)/1000)
			f_out = []
			index +=1
			for j in range(4):	# 4*pseudo-diploid with switch errors
				f_out.append(open(f_out_prefix+str(j),"w"))		
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,POS,STRING = macs_to_record_error(f)
			print 'len:',len(pos),len(POS),len(STRING)
			write_psmcfa_multiple(pos,f_out,i % int(args.nchr))
			f_msmc_out = open(macs_file+"_"+str(i)+"_error_%ik.txt" %(int(args.length)/1000),'w')
			write_msmcfa_multiple(POS,STRING,f_msmc_out,i % int(args.nchr))
		else:
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos,POS,STRING = macs_to_record_error(f)
			print 'len:',len(pos),len(POS),len(STRING)
			write_psmcfa_multiple(pos,f_out,i % int(args.nchr))
			f_msmc_out = open(macs_file+"_"+str(i)+"_error_%ik.txt" %(int(args.length)/1000),'w')
			write_msmcfa_multiple(POS,STRING,f_msmc_out,i % int(args.nchr))

def read_macs_single(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ.0" 
			f_out = open(f_out_prefix,"w")
			index +=1
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))

def read_macs_multiple(macs_file):	#diploid+4*pseudo-diploid
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

if __name__=="__main__":

	parser = argparse.ArgumentParser(description='macs to psmcfa file')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--nrep", default=1000,dest='nrep',help="number of repetition")
	parser.add_argument("--nchr", default=100,dest='nchr',help="number of chromosome")
	parser.add_argument("--nhap", default=2,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--version",default=None,dest='version')
	parser.add_argument("--mode",default=None,dest='mode')
	parser.add_argument("--freq",type=float,default=0.05,dest='freq')
	parser.add_argument("--length", type=int,dest='length',help="Expected length (in bp) of correctly phased segments")

	args = parser.parse_args()
	print args.length
	if args.mode=='unphase':
		read_macs_unphase(args.file)
	elif args.mode=='error':   # introduce phase error
		read_macs_error_with_msmc(args.file)
	else:
		read_macs_v2(args.file)
	