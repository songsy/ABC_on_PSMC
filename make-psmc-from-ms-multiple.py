# program for making psmc input file from ms output
# can take multiple samples
# python make-psmc-from-ms-multiple.py msfile #sample segment_length mode
import sys
import msparse_v2

hetChar = 'K'
homChar = 'T'
failChar = 'N'
winLen = 100

segmentLen = int(sys.argv[3])
#segmentLen = 10000000
msFileName = sys.argv[1]
sample = int(sys.argv[2])
hetOutFileName = msFileName.replace('.txt','.hetFQ')

msRecs = msparse_v2.read_msHOT_records(msFileName)

def ms_to_psmc_mode3(msRecs):     # diploid+pseudo-diploid, just hap1-hap1, 2 per population
	print 'Have %i records' % len(msRecs)
	print 'Have %i samples' % sample
	num_compare = sample+sample*(sample-1)/2
	for comparison in range(num_compare):
		outFile = open(hetOutFileName+'.'+str(comparison),'w')
		for rec_i in range(len(msRecs)):
			myHets = msparse_v2.convert_record_to_het_dict_v3(msRecs[rec_i],segmentLen)
			print rec_i,len(myHets),len(myHets[comparison])
			outFile.write('>%s\n' % (rec_i))
			windowStart = 0
			windowEnd = windowStart + winLen - 1
			windows = []
			while windowStart < segmentLen:
				windows.append(homChar)
				windowStart = windowEnd + 1
				windowEnd = windowStart + winLen - 1
			print 'Made %i widnows' % len(windows)			  
			#now, drop in the het sites into the windows
			for p in myHets[comparison]:
				wID = int(p/winLen)
				try:
					windows[wID] = hetChar
				except IndexError:
					print 'Error',p,winLen,wID,len(windows)
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
		print 'Record %i homWind %i hetWind %i' % (rec_i,numHom,numHet)
		outFile.close()

def ms_to_psmc_mode4(msRecs):     # diploid+pseudo-diploid*4, 2 per population
	print 'Have %i records' % len(msRecs)
	print 'Have %i samples' % sample
	num_compare = sample+sample*sample
	for comparison in range(num_compare):
		outFile = open(hetOutFileName+'.'+str(comparison),'w')
		for rec_i in range(len(msRecs)):
			myHets = msparse_v2.convert_record_to_het_dict_v4(msRecs[rec_i],segmentLen)
			print rec_i,len(myHets),len(myHets[comparison])
			outFile.write('>%s\n' % (rec_i))
			windowStart = 0
			windowEnd = windowStart + winLen - 1
			windows = []
			while windowStart < segmentLen:
				windows.append(homChar)
				windowStart = windowEnd + 1
				windowEnd = windowStart + winLen - 1
			print 'Made %i widnows' % len(windows)			  
			#now, drop in the het sites into the windows
			for p in myHets[comparison]:
				wID = int(p/winLen)
				try:
					windows[wID] = hetChar
				except IndexError:
					print 'Error',p,winLen,wID,len(windows)
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
		print 'Record %i homWind %i hetWind %i' % (rec_i,numHom,numHet)
		outFile.close()

def ms_to_psmc_mode2(msRecs):     # diploid
	print 'Have %i records' % len(msRecs)
	print 'Have %i samples' % sample
	num_compare = sample
	for comparison in range(num_compare):
		outFile = open(hetOutFileName+'.'+str(comparison),'w')
		for rec_i in range(len(msRecs)):
			myHets = msparse_v2.convert_record_to_het_dict_v2(msRecs[rec_i],segmentLen)
			print rec_i,len(myHets),len(myHets[comparison])
			outFile.write('>%s\n' % (rec_i))
			windowStart = 0
			windowEnd = windowStart + winLen - 1
			windows = []
			while windowStart < segmentLen:
				windows.append(homChar)
				windowStart = windowEnd + 1
				windowEnd = windowStart + winLen - 1
			print 'Made %i widnows' % len(windows)			  
			#now, drop in the het sites into the windows
			for p in myHets[comparison]:
				wID = int(p/winLen)
				windows[wID] = hetChar
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
		print 'Record %i homWind %i hetWind %i' % (rec_i,numHom,numHet)
		outFile.close()

if __name__=="__main__":
	mode = sys.argv[4]
	if mode=="diploid":
		ms_to_psmc_mode2(msRecs)
	elif mode=="all":
		ms_to_psmc_mode4(msRecs)
	else:
		ms_to_psmc_mode3(msRecs)
