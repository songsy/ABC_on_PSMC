import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

pedIndFile = 'grmjunk.id'
evalFile = 'hapmap3.eval'
evecFile = 'hapmap3.evec'

print evecFile
print evalFile 
sampleFile = '/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/relationships_w_pops_121708.txt'
# read in color file
nameToPop={}
PopToColor = {'CEU':'orange','TSI':'brown','YRI':'darkgreen','MKK':'greenyellow','LWK':'green','JPT':'purple','CHB':'blue','CHD':'pink','ASW':'gold','MEX':'grey','GIH':'red'}
inFile = open(sampleFile,'r')
for line in inFile:
	line = line.rstrip()
	line = line.split('\t')
	nameToPop[line[1]] = line[6]

inFile.close()

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
for i in range(20):
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
#	 print line
	sn = line[0]
	sn = sn.split(':')[0]
	for i in range(1,3):
		n = float(line[i])
		pcCoord[(sn,i)] = n

inFile.close()
print 'read in PC'
 
fig = plt.figure()
ax = fig.add_subplot(111)


pc1 = 1
pc2 = 2

didLeg = {}
didPop=[]
for i in range(len(longNames)):
	ln = longNames[i]
	try:
		col = PopToColor[nameToPop[ln]]
		Label=nameToPop[ln]
	except KeyError:
		continue
	x = pcCoord[(ln,pc1)]
	y = pcCoord[(ln,pc2)]
	
	print ln,x,y,col
	
	spec = ln.split('-')[0]
	if Label not in didPop:
		didPop.append(Label)
		ax.scatter(x,y,color=col,label=Label,marker=".")
	else:
		ax.scatter(x,y,color=col,marker=".")

xLab = 'PC%i (%.1f' % (pc1,eVals[pc1-1]*100)
xLab += '%)'
plt.xlabel(xLab)

yLab = 'PC%i (%.1f' % (pc2,eVals[pc2-1]*100)
yLab += '%)'
plt.ylabel(yLab)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles, labels,scatterpoints=1,prop={'size':8},loc='lower right')

savePlot = True
if savePlot is False:
	plt.show()
else:
	fileName = 'hapmap-labeled-pc1_2.png'
	plt.savefig(fileName)	 


	