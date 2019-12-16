# utilities for parsing output of ms
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
            rec = {}
            line = inFile.readline() #num of segsites
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
    for i in range(1,len(varPos)):
        if varPos[i] == varPos[i-1]:
            varPos[i] = varPos[i] + 1
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
def convert_record_to_het_dict(msRecord,segmentLen):
    mySeqs = convert_record_to_var_dict(msRecord,segmentLen)
    numSeq = len(mySeqs) - 1
    myHets = []
    varPos = mySeqs[-1]
#    print 'len of varPos is',len(varPos)
#    t = {}
#    for i in varPos:
#        t[i] = 1
#    print 'len of t is', len(t)
#    print varPos
    if numSeq % 2 != 0:
        print 'Have %i to combine' % numSeq
        print 'NOT EVEN, do not know what to do'
        sys.exit()
    for seq_i in range(0,numSeq-1,2):
        seqA_i = seq_i
        seqB_i = seq_i + 1
        seqA = mySeqs[seqA_i]
        seqB = mySeqs[seqB_i]
#        print len(seqA)
#        print len(seqB)
        
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
#het =  'e' hom = 'o'
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
            
        

    







    