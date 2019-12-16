import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mpatches

pedIndFile = '/home/jmkidd/kidd-lab-scratch/shiya-projects/PCA/Illumina_chip_1KG_phase3_HapMap_Af_Eu_no_drop.id'
evalFile = 'Illumina_chip_1KG_phase3_HapMap_Af_Eu_no_drop.eval'
evecFile = 'Illumina_chip_1KG_phase3_HapMap_Af_Eu_no_drop.evec'

print evecFile
print evalFile 
sampleFile = '/home/jmkidd/kidd-lab-scratch/shiya-projects/HGDP/Illumina_chip/Illumina_chip_1KG_phase3_HapMap.ind'
# read in color file
nameToPop={}
#regionToPop={'Bantu':'Afr','Mandenka':'Afr','Yoruba':'YRI','San':'San','Mbuti':'Mbuti','Biaka':'Afr','Mozabite':'Afr','Orcadian':'Eur','Adygei':'Eur','Russian':'Eur','Basque':'Eur','French':'Eur','North':'Eur','Sardinian':'Eur','Tuscan':'Eur','Bedouin':'AsianW','Druze':'AsianW','Palestinian':'AsianW','Balochi':'AsianC','Brahui':'AsianC','Makrani':'AsianC','Sindhi':'AsianC','Pathan':'AsianC','Burusho':'AsianC','Hazara':'AsianC','Uygur':'AsianC','Kalash':'AsianC','Han':'AsianE','Dai':'AsianE','Daur':'AsianE','Hezhen':'AsianE','Lahu':'AsianE','Miaozu':'AsianE','Oroqen':'AsianE','She':'AsianE','Tujia':'AsianE','Tu':'AsianE','Xibo':'AsianE','Yizu':'AsianE','Mongola':'AsianE','Naxi':'AsianE','Cambodian':'AsianE','Japanese':'AsianE','Yakut':'AsianE','NAN':'Ocean','Papuan':'Ocean','Karitiana':'NatAm','Surui':'NatAm','Colombian':'NatAm','Maya':'NatAm','Pima':'NatAm'}
#PopToColor = {'Eur':'orange','YRI':'darkgreen','MKK':'greenyellow','LWK':'green','GIH':'red','Afr':'brown','ESN':'Olive','GWD':'Teal','MSL':'SpringGreen','CEU':'Tomato','GIH':'Yellow','San':'Purple','Mbuti':'Violet','AsianW':'Cyan','AsianC':'DodgerBlue','AsianE':'MediumBlue','Ocean':'Pink','NatAm':'Gray'}
PopToColor = {'Bantu':'Pink','Mandenka':'Magenta','San':'Purple','Mbuti':'Violet','Biaka':'Red','Mozabite':'IndianRed','Orcadian':'Aqua','Adygei':'Navy','Russian':'CadetBlue','Basque':'LightCyan','French':'DodgerBlue','North':'Chocolate','Sardinian':'Maroon','Tuscan':'Brown','YRI':'darkgreen','MKK':'greenyellow','LWK':'green','ESN':'Olive','GWD':'Teal','MSL':'SpringGreen','CEU':'Tomato','GIH':'Yellow'}
inFile = open(sampleFile,'r')
for line in inFile:
	line = line.rstrip()
	line = line.split('\t')
	nameToPop[line[1]] = line[3]

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
didPop=[]
for i in range(len(longNames)):
	ln = longNames[i]
	Label=""
	try:
		col = PopToColor[nameToPop[ln]]
		Label=nameToPop[ln]
	except KeyError:
#		col = PopToColor[regionToPop[nameToPop[ln]]]
#		Label=regionToPop[nameToPop[ln]]
		print ln
	x = pcCoord[(ln,pc1)]
	y = pcCoord[(ln,pc2)]
	
	print ln,x,y,col
	
	spec = ln.split('-')[0]
	if Label not in didPop:
		didPop.append(Label)
		if Label=='North':
			Label='Italian'
		ax.scatter(x,y,color=col,label=Label,marker=".")
	else:
		ax.scatter(x,y,color=col,marker=".")
'''
plt.annotate('NA21302', xy=(-0.0048,-0.0009),  xycoords='data',
                xytext=(-50, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )

plt.annotate('NA19240', xy=(-0.0226,-0.0012),  xycoords='data',
                xytext=(-20, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )

plt.annotate('NA20847', xy=(0.0330,0.0975),  xycoords='data',
                xytext=(20, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )

plt.annotate('NA12878', xy=(0.0405,-0.0245),  xycoords='data',
                xytext=(-50, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )
plt.annotate('HG02799', xy=(-0.0208,-0.0018),  xycoords='data',
                xytext=(-40, -30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )

plt.annotate('HG03108', xy=(-0.0226,-0.0007),  xycoords='data',
                xytext=(10, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )        

plt.annotate('HG03428', xy=(-0.0220,-0.0006),  xycoords='data',
                xytext=(30, -30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )
plt.annotate('HGDP01029', xy=(-0.0202,0.0020),  xycoords='data',
                xytext=(30, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )

plt.annotate('HGDP00456', xy=(-0.0220,0.0042),  xycoords='data',
                xytext=(-25, 25), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=6
                )        
'''

xLab = 'PC%i (%.1f' % (pc1,eVals[pc1-1]*100)
xLab += '%)'
plt.xlabel(xLab)
yLab = 'PC%i (%.1f' % (pc2,eVals[pc2-1]*100)
yLab += '%)'
plt.ylabel(yLab)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles, labels,scatterpoints=1,prop={'size':6},loc='upper left',)

'''
patch=[]
for i in PopToColor.keys():
	red_patch = mpatches.Patch(color=PopToColor[i], label=i)
	patch.append(red_patch)
plt.legend(handles=patch)
'''
savePlot = True
if savePlot is False:
	plt.show()
else:
	fileName = 'Illumina-Af-Eu-no-drop-colorful-pc3_4.pdf'
	plt.savefig(fileName,format='pdf')	 


	