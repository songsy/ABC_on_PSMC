import sys
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

def macs_to_msmc(macs_file,L,sample,hap):
	num_compare = sample+sample*(sample-1)/2
	f_out = []
	prev_pos = []
	pos = []
	first = []
	for i in range(num_compare):
		f_out.append(open(macs_file.replace('.txt','_%i.txt' %(i)),'w'))
		pos.append(0)
		prev_pos.append(0)
		first.append(True)
	f = open(macs_file,'r')
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			for i in range(sample):
				substring = col[4][i*hap:(i+1)*hap]
				if have_het(substring) is True:
					if first[index] is True:
						pos[index] = int(L*float(col[2]))
						prev_pos[index] = pos[index]
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index],substring)
						first[index] = False
					else:
						pos[index] = int(L*float(col[2]))
						if pos[index]<=prev_pos[index]:
							pos[index] =prev_pos[index]+1
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index]-prev_pos[index],substring)
						prev_pos[index] = pos[index]
				index +=1
			for i in range(sample):
				for j in range(i+1,sample):
					substring = ''
#					print i,j
					for k in range(hap):
						seq_id=i*hap+k
#						print seq_id
						substring+=col[4][seq_id]
					for k in range(hap):
						seq_id=j*hap+k
#						print seq_id
						substring+=col[4][seq_id]
					if have_het(substring) is True:
						if first[index] is True:
							pos[index] = int(L*float(col[2]))
							prev_pos[index] = pos[index]
							substring=code_to_string(substring)
							print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index],substring)
							first[index] = False
						else:
							pos[index] = int(L*float(col[2]))
							if pos[index]<=prev_pos[index]:
								pos[index] =prev_pos[index]+1
							substring=code_to_string(substring)
							print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index]-prev_pos[index],substring)
							prev_pos[index] = pos[index]
					index +=1

def macs_to_msmc_unphase(macs_file,L,sample,hap):
	f_out = []
	f_out2 = []
	prev_pos = []
	pos = []
	pos2 = []
	prev_pos2 = []
	first = []
	first2 = []
	unphase_num =[]
	unphase_num2 =[]
	freq = [0.01,0.02,0.05,0.1,0.15]
	for i in [0.01,0.02,0.05,0.1,0.15]:
		f_out.append(open(macs_file.replace('.txt','_%.2f.txt' %(i)),'w'))
		pos.append(0)
		prev_pos.append(0)
		first.append(True)
		unphase_num.append(0)
	for i in [0.01,0.02,0.05,0.1,0.15]:
		f_out2.append(open(macs_file.replace('.txt','_%.2f_skip.txt' %(i)),'w'))
		pos2.append(0)
		prev_pos2.append(0)
		first2.append(True)
		unphase_num2.append(0)
	f = open(macs_file,'r')
	tot = 0
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			if col[4][0]==col[4][1] and col[4][2]==col[4][3]:
				s = 1
			else:
				s=random.uniform(0,1)
				tot +=1
			for index in range(len(freq)):
				if s<=freq[index]:
					unphase_num[index]+=1
					substring = col[4]
					unphase_num2[index]+=1
					if col[4][0]!=col[4][1]:
						substring2=col[4][1]+col[4][0]+col[4][2]+col[4][3]
					elif col[4][2]!=col[4][3]:
						substring2=col[4][0]+col[4][1]+col[4][3]+col[4][2]
					else:
						print line
					if first[index] is True:
						pos[index] = int(L*float(col[2]))
						prev_pos[index] = pos[index]
						substring=code_to_string(substring)
						substring2=code_to_string(substring2)
						print >>f_out[index],'1\t%i\t%i\t%s,%s' %(pos[index],pos[index],substring,substring2)
						first[index] = False
					else:
						pos[index] = int(L*float(col[2]))
						if pos[index]<=prev_pos[index]:
							pos[index] =prev_pos[index]+1
						substring=code_to_string(substring)
						substring2=code_to_string(substring2)
						print >>f_out[index],'1\t%i\t%i\t%s,%s' %(pos[index],pos[index]-prev_pos[index],substring,substring2)
						prev_pos[index] = pos[index]
				else:
					substring = col[4]
					if first[index] is True:
						pos[index] = int(L*float(col[2]))
						prev_pos[index] = pos[index]
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index],substring)
						first[index] = False
					else:
						pos[index] = int(L*float(col[2]))
						if pos[index]<=prev_pos[index]:
							pos[index] =prev_pos[index]+1
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index]-prev_pos[index],substring)
						prev_pos[index] = pos[index]
					if first2[index] is True:
						pos2[index] = int(L*float(col[2]))
						prev_pos2[index] = pos2[index]
						print >>f_out2[index],'1\t%i\t%i\t%s' %(pos2[index],pos2[index]-unphase_num2[index],substring)
						first2[index] = False
					else:
						pos2[index] = int(L*float(col[2]))
						if pos2[index]<=prev_pos2[index]:
							pos2[index] =prev_pos2[index]+1
						print >>f_out2[index],'1\t%i\t%i\t%s' %(pos2[index],pos2[index]-prev_pos2[index]-unphase_num2[index],substring)
						prev_pos2[index] = pos2[index]
					unphase_num2[index]=0
	for i in range(len(unphase_num)):
		print '%.2f' %(float(unphase_num[i])/tot)

def macs_to_msmc_multiple(macs_file,L,sample,hap):
	num_compare = 1
	f_out = []
	prev_pos = []
	pos = []
	first = []
	for i in range(num_compare):
		f_out.append(open(macs_file.replace('.txt','_%i.txt' %(i)),'w'))
		pos.append(0)
		prev_pos.append(0)
		first.append(True)
	f = open(macs_file,'r')
	a=random.randint(0,sample/2-1)
	b=random.randint(0,sample/2-1)
	print a,b
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			if True:
				substring = col[4][a*hap]+col[4][a*hap+1]+col[4][sample+b*hap]+col[4][sample+b*hap+1]
				if have_het(substring) is True:
					if first[index] is True:
						pos[index] = int(L*float(col[2]))
						prev_pos[index] = pos[index]
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index],substring)
						first[index] = False
					else:
						pos[index] = int(L*float(col[2]))
						if pos[index]<=prev_pos[index]:
							pos[index] =prev_pos[index]+1
						substring=code_to_string(substring)
						print >>f_out[index],'1\t%i\t%i\t%s' %(pos[index],pos[index]-prev_pos[index],substring)
						prev_pos[index] = pos[index]
				index +=1

if __name__=="__main__":
	macs_file = sys.argv[1]
	length = int(sys.argv[2])
	sample = int(sys.argv[3])    # number of sample
	hap = int(sys.argv[4])		# number of haplotype per individual
	mode = 'phase'
	if mode=='unphase':
		macs_to_msmc_unphase(macs_file,length,sample,hap)
	elif sample==2:
		macs_to_msmc(macs_file,length,sample,hap)
	else:
		macs_to_msmc_multiple(macs_file,length,sample,hap)
	
