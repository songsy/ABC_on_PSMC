# utilities for parsing output of ms
# modified to take more than 2 samples
import sys

##############################################################################
def init_blank_list(listLen,element):
	myList = []
	for i in range(listLen):
		myList.append(element)
	return myList
###############################################################################
def read_msHOT_records(msFileName):
	print msFileName
	msRecords = []
	inFile = open(msFileName)
	# read in to first record
	line = inFile.readline()
	line = inFile.readline()
	line = inFile.readline()	
	while True:
		line = inFile.readline()
		if line == '':
			break
		# start of record
		if line == '//\n':
			while True:
				line = inFile.readline()
				if line[0]!='[':
					break
			rec = {}
#			line = inFile.readline() #num of segsites
			line = inFile.readline() # blank line
			line = inFile.readline()
			line = line.rstrip()
			line = line.split()
			line = line[1:]
			line = [float(i) for i in line]
			rec['pos'] = line
			rec['haps'] = []
			while True:
				line = inFile.readline()
				if line == '\n' or line == '':
					break
				line = line.rstrip()
				rec['haps'].append(line)
			msRecords.append(rec)			 
	inFile.close()
	return msRecords
###############################################################################
# will go from 0 to ... <segmentLen
def convert_record_to_fasta(msRecord,segmentLen):
	varPos = []
	for i in msRecord['pos']:
		p = int(i * segmentLen)
		varPos.append(p)
	mySeqs = []
	for hap in msRecord['haps']:
		seq = init_blank_list(segmentLen+1,0)  # pos can be 1.0 --> gives us the segment len
		for i in range(len(varPos)):
			index = varPos[i]
			if hap[i] == '1':
				seq[index] = 1
		mySeqs.append(seq)
	return mySeqs
###############################################################################
# will go from 0 to ... <segmentLen
def convert_record_to_var_dict(msRecord,segmentLen):
	varPos = []
	for i in msRecord['pos']:
		p = int(i * segmentLen)
		varPos.append(p)
	
	# yikes -- need to check to see if we have SNPs that are at same position
	# i changed the precision of the float ouput in msHOT, but we can still have collisions
	numAdjust = 0
	if varPos[0]==0:
		varPos[0]=1
	for i in range(1,len(varPos)):
		if varPos[i] == varPos[i-1]:
			varPos[i] = varPos[i] + 1
			numAdjust += 1
			if varPos[i] > segmentLen:
				varPos[i] = segmentLen
		elif varPos[i] < varPos[i-1]:
			varPos[i] = varPos[i-1] + 1
			numAdjust += 1
			if varPos[i] > segmentLen:
				varPos[i] = segmentLen
	print 'had to adjust',numAdjust,'SNP pos'
	
	mySeqs = []
	for hap in msRecord['haps']:
		hapVarPos = {}
		for i in range(len(varPos)):
			index = varPos[i]
			if hap[i] == '1':
				hapVarPos[index] = 1
		mySeqs.append(hapVarPos)
	# put the varPos in at the end of mySeqs
	mySeqs.append(varPos)
	return mySeqs
###############################################################################
# use dictionaries to record sites of het
# for single chromosome, pair-wise comparison
def convert_record_to_het_dict(msRecord,segmentLen):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myHets = []
	varPos = mySeqs[-1]
#	 print 'len of varPos is',len(varPos)
#	 t = {}
#	 for i in varPos:
#		 t[i] = 1
#	 print 'len of t is', len(t)
#	 print varPos
#	 if numSeq % 2 != 0:
#		print 'Have %i to combine' % numSeq
#		 print 'NOT EVEN, do not know what to do'
#		 sys.exit()
	for seq_i in range(numSeq):
		for seq_j in range(seq_i+1,numSeq):
			seqA = mySeqs[seq_i]
			seqB = mySeqs[seq_j]		
			myVars = {}
			for i in varPos:
				if i not in seqA and i not in seqB:
					continue
				if i in seqA and i in seqB:
					continue
				myVars[i] = 1
			myHets.append(myVars)
	return myHets
###############################################################################
# use dictionaries to record sites of het
# for diploid individual, het of two chromosome of same individual
def convert_record_to_het_dict_v2(msRecord,segmentLen):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myHets = []
	varPos = mySeqs[-1]
#	 print 'len of varPos is',len(varPos)
#	 t = {}
#	 for i in varPos:
#		 t[i] = 1
#	 print 'len of t is', len(t)
#	 print varPos
	if numSeq % 2 != 0:
		print 'Have %i to combine' % numSeq
		print 'NOT EVEN, do not know what to do'
		sys.exit()
	for seq_i in range(0,numSeq-1,2):
		seqA_i = seq_i
		seqB_i = seq_i + 1
		seqA = mySeqs[seqA_i]
		seqB = mySeqs[seqB_i]
#		 print len(seqA)
#		 print len(seqB)
		
		myVars = {}
		for i in varPos:
			if i not in seqA and i not in seqB:
				continue
			if i in seqA and i in seqB:
					continue
			myVars[i] = 1
		myHets.append(myVars)
	return myHets
	
###############################################################################
# use dictionaries to record sites of het
# for diploid individual, first generate het of two chromosome of same individual, then generate het of two chromosome 
# from different individual
def convert_record_to_het_dict_v3(msRecord,segmentLen):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myHets = []
	varPos = mySeqs[-1]
#	 print 'len of varPos is',len(varPos)
#	 t = {}
#	 for i in varPos:
#		 t[i] = 1
#	 print 'len of t is', len(t)
#	 print varPos
#	 if numSeq % 2 != 0:
#		print 'Have %i to combine' % numSeq
#		 print 'NOT EVEN, do not know what to do'
#		 sys.exit()
	for seq_i in range(0,numSeq-1,2):
		seqA_i = seq_i
		seqB_i = seq_i + 1
		seqA = mySeqs[seqA_i]
		seqB = mySeqs[seqB_i]
#		 print len(seqA)
#		 print len(seqB)
		
		myVars = {}
		for i in varPos:
			if i not in seqA and i not in seqB:
				continue
			if i in seqA and i in seqB:
					continue
			myVars[i] = 1
		myHets.append(myVars)
	for seq_i in range(numSeq/2):
		for seq_j in range(seq_i+1,numSeq/2):
			seqA = mySeqs[2*seq_i]
			seqB = mySeqs[2*seq_j]		
			myVars = {}
			for i in varPos:
				if i not in seqA and i not in seqB:
					continue
				if i in seqA and i in seqB:
					continue
				myVars[i] = 1
			myHets.append(myVars)
	return myHets

###############################################################################
# use dictionaries to record sites of het
# for diploid individual, first generate het of two chromosome of same individual, then generate het of two chromosome 
# from different individual, all combinations
def convert_record_to_het_dict_v4(msRecord,segmentLen):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myHets = []
	varPos = mySeqs[-1]
#	 print 'len of varPos is',len(varPos)
#	 t = {}
#	 for i in varPos:
#		 t[i] = 1
#	 print 'len of t is', len(t)
#	 print varPos
#	 if numSeq % 2 != 0:
#		print 'Have %i to combine' % numSeq
#		 print 'NOT EVEN, do not know what to do'
#		 sys.exit()
	for seq_i in range(0,numSeq-1,2):
		seqA_i = seq_i
		seqB_i = seq_i + 1
		seqA = mySeqs[seqA_i]
		seqB = mySeqs[seqB_i]
#		 print len(seqA)
#		 print len(seqB)
		
		myVars = {}
		for i in varPos:
			if i not in seqA and i not in seqB:
				continue
			if i in seqA and i in seqB:
					continue
			myVars[i] = 1
		myHets.append(myVars)
	for seq_i in range(2):
		for seq_j in range(2):
			seqA = mySeqs[seq_i]
			seqB = mySeqs[2+seq_j]		
			myVars = {}
			for i in varPos:
				if i not in seqA and i not in seqB:
					continue
				if i in seqA and i in seqB:
					continue
				myVars[i] = 1
			myHets.append(myVars)
	return myHets

###############################################################################
# takes seq 0,1 and 2,3, 4,5 etc and returns sequence of het/hom calls
#het =	'e' hom = 'o'
def convert_record_to_fasta_het(msRecord,segmentLen):
	mySeqs = convert_record_to_fasta(msRecord,segmentLen)
	numSeq = len(mySeqs)
	myHets = []
	if numSeq % 2 != 0:
		print 'Have %i to combine' % numSeq
		print 'NOT EVEN, do not know what to do'
		sys.exit()
	for seq_i in range(0,numSeq,2):
		seqA_i = seq_i
		seqB_i = seq_i + 1
		hetGen = init_blank_list(segmentLen+1,'?')
		seqA = mySeqs[seqA_i]
		seqB = mySeqs[seqB_i]
		for i in range(segmentLen):
			if seqA[i] == seqB[i]:
				hetGen[i] = 'o'
			else:
				hetGen[i] = 'e'
		myHets.append(hetGen)
	return myHets
###############################################################################
			
def convert_record_to_var_dict_v2(msRecord,segmentLen,hapnum,samplenum):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myCalls = []
	varPos = mySeqs[-1]
	hom = 0
	for seq_i in range(samplenum):
		seq=[]
		for j in range(hapnum):
			seq_id=seq_i*hapnum+j
			seq.append(mySeqs[seq_id])		
		myVars = []
		oldPos = 0
		for i in varPos:
			m = []
			het = False
			for j in range(hapnum):
				m.append('A')
			for j in range(hapnum):
				if i in seq[j]:
					m[j]='T'
					het = True
			if het is True:
				myVars.append(str(i)+'\t'+str(i-oldPos)+'\t'+''.join(m))
				oldPos = i
		myCalls.append(myVars)
	for seq_i in range(samplenum):
		for seq_j in range(seq_i+1,samplenum):
			seq=[]
			for j in range(hapnum):
				seq_id=seq_i*hapnum+j
				seq.append(mySeqs[seq_id])
			for j in range(hapnum):
				seq_id=seq_j*hapnum+j
				seq.append(mySeqs[seq_id])
			myVars = []
			oldPos = 0
			for i in varPos:
				m = []
				het = False
				for j in range(hapnum*2):
					m.append('A')
				for j in range(hapnum*2):
					if i in seq[j]:
						m[j]='T'
						het = True
				if het is True:
					myVars.append(str(i)+'\t'+str(i-oldPos)+'\t'+''.join(m))
					oldPos = i
			myCalls.append(myVars)
	return myCalls

def convert_record_to_var_dict_v3(msRecord,segmentLen,hapnum,samplenum):
	mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
	numSeq = len(mySeqs) - 1
	myCalls = []
	varPos = mySeqs[-1]
	hom = 0
	for seq_i in range(samplenum):
		seq=[]
		for j in range(hapnum):
			seq_id=seq_i*hapnum+j
			seq.append(mySeqs[seq_id])		
	myVars = []
	for i in varPos:
		code=[]
		if i in mySeqs[samplenum-1]:
			for seq_i in range(samplenum-1):
				seq_id=seq_i*hapnum
				if i in mySeqs[seq_id] and i in mySeqs[seq_id+1]:
					code.append(0)
				elif i in mySeqs[seq_id] and i not in mySeqs[seq_id+1]:
					code.append(1)
				elif i not in mySeqs[seq_id] and i in mySeqs[seq_id+1]:
					code.append(1)
				elif i not in mySeqs[seq_id] and i not in mySeqs[seq_id+1]:
					code.append(2)
		else:
			for seq_i in range(samplenum-1):
				seq_id=seq_i*hapnum
				if i in mySeqs[seq_id] and i in mySeqs[seq_id+1]:
					code.append(2)
				elif i in mySeqs[seq_id] and i not in mySeqs[seq_id+1]:
					code.append(1)
				elif i not in mySeqs[seq_id] and i in mySeqs[seq_id+1]:
					code.append(1)
				elif i not in mySeqs[seq_id] and i not in mySeqs[seq_id+1]:
					code.append(0)
		code.append(0)
		myCalls.append(code)
	return myCalls,varPos




	