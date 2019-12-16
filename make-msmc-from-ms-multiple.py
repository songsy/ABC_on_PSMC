# program for making psmc input file from ms output
# can take multiple samples
# python make-msmc-from-ms-multiple.py msfile #sample #hap segment_length

import sys
import msparse_v2

hetChar = 'K'
homChar = 'T'
failChar = 'N'
winLen = 100

segmentLen = int(sys.argv[4])
msFileName = sys.argv[1]
sample = int(sys.argv[2])    # number of sample
hap = int(sys.argv[3])		# number of haplotype per individual
hetOutFileName = msFileName.replace('.txt','.hetFQ')

msRecs = msparse_v2.read_msHOT_records(msFileName)

def ms_to_msmc_mode3(msRecs):     # diploid+pseudo-diploid, 2 per population
	print 'Have %i records' % len(msRecs)
	print 'Have %i samples' % sample
	num_compare = sample+sample*(sample-1)/2
	for comparison in range(num_compare):
		cmds = 'msmc --fixedRecombination --skipAmbiguous -o ' + msFileName.replace('.txt','.'+str(comparison))
		for rec_i in range(len(msRecs)):
			outFile = open(hetOutFileName+'.'+str(comparison)+'.'+str(rec_i+1)+'.txt','w')
			cmds += ' ' + hetOutFileName+'.'+str(comparison)+'.'+str(rec_i+1)+'.txt'
			myHets = msparse_v2.convert_record_to_var_dict_v2(msRecs[rec_i],segmentLen,hap,sample)
			print rec_i,len(myHets),len(myHets[comparison])
			for p in myHets[comparison]:
				outFile.write(str(rec_i+1)+'\t'+p+'\n')
		outFile.close()
		print cmds

if __name__=="__main__":
	ms_to_msmc_mode3(msRecs)
