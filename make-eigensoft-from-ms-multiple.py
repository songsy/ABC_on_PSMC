# program for making psmc input file from ms output
# can take multiple samples
# python make-eigensoft-from-ms-multiple.py msfile #sample #hap segment_length
# Last one is the chimp one
# used for making input for calculating d-Stats
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

def ms_to_eigensoft_mode3(msRecs):     # diploid, 2 per population
	print 'Have %i records' % len(msRecs)
	print 'Have %i samples' % sample
	f_geno=open(msFileName.replace('.txt','.geno'),'w')
	f_snp=open(msFileName.replace('.txt','.snp'),'w')
	for rec_i in range(len(msRecs)):
		myHets,varPos = msparse_v2.convert_record_to_var_dict_v3(msRecs[rec_i],segmentLen,hap,sample)
		for p in range(len(myHets)):
			f_geno.write("".join(map(str,myHets[p]))+'\n')
			f_snp.write('\t%s_%s\t1\t%f\t%s\t%s\t%s\n' %(rec_i,varPos[p],float(varPos[p])/segmentLen,varPos[p],'A','T'))
	

if __name__=="__main__":
	ms_to_eigensoft_mode3(msRecs)
