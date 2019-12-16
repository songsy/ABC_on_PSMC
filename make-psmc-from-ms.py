# program for making psmc input file from ms output
# assumes that only 2 chromosomes have been sampled from ms
import sys
import msparse
print 'hi'



hetChar = 'K'
homChar = 'T'
failChar = 'N'
winLen = 100

segmentLen = 10000000
#segmentLen = 200

msFileName = sys.argv[1]
#msFileName = 'sim1.mshot.t1'
#msFileName = 'ms.test'
hetOutFileName = msFileName + '.hetFQ'


msRecs = msparse.read_msHOT_records(msFileName)



print 'Have %i records' % len(msRecs)

outFile = open(hetOutFileName,'w')
for rec_i in range(len(msRecs)):
    myHets = msparse.convert_record_to_het_dict(msRecs[rec_i],segmentLen)
    print rec_i,len(myHets),len(myHets[0])
    outFile.write('>%s\n' % (rec_i))
    windowStart = 0
    windowEnd = windowStart + winLen - 1
    windows = []
    while windowStart <= (segmentLen + 1):
        windows.append(homChar)
        windowStart = windowEnd + 1
        windowEnd = windowStart + winLen - 1
    print 'Made %i widnows' % len(windows)            
    #now, drop in the het sites into the windows
    for p in myHets[0]:
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
