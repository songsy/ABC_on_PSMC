import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

pedIndFile = 'H11_all_pop.id'
evalFile = 'H11.eval'
evecFile = 'H11.evec'

print evecFile
print evalFile 
sampleFile = '/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/relationships_w_pops_121708.txt'
# read in color file
'''
nameToPop={}
PopToColor = {'CEU':'orange','TSI':'brown','YRI':'darkgreen','MKK':'greenyellow','LWK':'green','JPT':'purple','CHB':'blue','CHD':'pink','ASW':'gold','MEX':'grey','GIH':'red'}
inFile = open(sampleFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    nameToPop[line[1]] = line[6]

inFile.close()
'''

nameToColor={'NA12878':'orange','NA21302':'greenyellow','NA21732':'greenyellow','NA21733':'greenyellow','NA21737':'greenyellow','NA21767':'greenyellow','NA19240':'darkgreen','HG02799':'green','HG03108':'lightgreen','HG03428':'aquamarine','HGDP01029':'blue','chimp':'red'}
longNames = []
inFile = open(pedIndFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    longName = line[1]
    longNames.append(longName)

inFile.close()

eVals = []
tot = 0.0
inFile = open(evalFile)
for line in inFile:
    line = line.rstrip()
    line = line.split()
    n = float(line[0])
    tot += n
    eVals.append(n)

inFile.close()
for i in range(10):
    eVals[i] = eVals[i] / tot
    print i,eVals[i]

pcCoord = {}    
# read in PC positions   
inFile = open(evecFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[0][0] == '#':
        continue
#    print line
    sn = line[0]
    sn = sn.split(':')[0]
    for i in range(1,6):
        n = float(line[i])
        pcCoord[(sn,i)] = n

inFile.close()
print 'read in PC'
 
fig = plt.figure()
ax = fig.add_subplot(111)


pc1 = 3
pc2 = 4

didLeg = {}

for i in range(len(longNames)):
    ln = longNames[i]
    try:
    	col = nameToColor[ln]
    except KeyError:
    	continue
    x = pcCoord[(ln,pc1)]
    y = pcCoord[(ln,pc2)]
    
    print ln,x,y,col
    
    spec = ln.split('-')[0]
    ax.scatter(x,y,color=col)

xLab = 'PC%i (%.1f' % (pc1,eVals[pc1-1]*100)
xLab += '%)'
plt.xlabel(xLab)

yLab = 'PC%i (%.1f' % (pc2,eVals[pc2-1]*100)
yLab += '%)'
plt.ylabel(yLab)


savePlot = True
if savePlot is False:
    plt.show()
else:
    fileName = 'H11-labeled-pc3_4.png'
    plt.savefig(fileName)    


    