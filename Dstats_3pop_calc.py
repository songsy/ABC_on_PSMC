'''
f=open('../plink3/all101_CTC_HXH_Freedman_SampleInfo_label.txt','r')
breed_list=set()
grey_wolf_list=set()
ancient_list=set()
village_dog_list=set()
Jackel="GoldenJackal"
for line in f:
	line = line.strip()
	col=line.split('\t')
	name="%s_%s" %(col[1],col[0])
#	print name
	if col[2]=='Breed':
		breed_list.add(col[3])
	elif col[2]=='AncientDog':
		ancient_list.add(col[3])
	elif col[2]=='VillageDog':
		village_dog_list.add(col[3])
	elif col[2]=='Wolf':
		grey_wolf_list.add(col[3])
'''
f=open('../plink3/all101_CTC_HXH_Freedman_SampleInfo_label.txt','r')
breed_list=set()
grey_wolf_list=set()
ancient_list=set()
village_dog_list=set()
Jackel="GoldenJackal"
RedWolf="wolfRed"
for line in f:
	line = line.strip()
	col=line.split('\t')
	name="%s_%s" %(col[1],col[0])
#	print name
#	if col[2]=='Breed':
#		breed_list.add(col[3])
	if col[2]=='AncientDog':
		ancient_list.add(col[3])
#	elif col[2]=='VillageDog':
#		village_dog_list.add(col[3])
	elif col[2]=='Wolf':
		grey_wolf_list.add(col[3])

f=open('../plink4/ShannonBoyko_All_label.txt','r')
for line in f:
	line = line.strip()
	col=line.split('\t')
	if col[2]=='Breed':
		breed_list.add(col[3])
	elif col[2]=='VillageDog':
		village_dog_list.add(col[3])
	elif col[2]=='Wolf':
		grey_wolf_list.add(col[3])

breed_list=list(breed_list)
village_dog_list=list(village_dog_list)
ancient_list=list(ancient_list)
grey_wolf_list=list(grey_wolf_list)

dog_list=ancient_list+breed_list+village_dog_list

f=open('Shannon_GoldenJackal_qp3pop_dstats.txt','r')
all_list=dog_list+grey_wolf_list
matrix=[[0 for i in range(len(all_list))] for j in range(len(all_list))]
#print len(matrix[0]),len(matrix),len(all_list)
for line in f:
	line = line.rstrip().split('\t')
	a=all_list.index(line[0])
	b=all_list.index(line[1])
	matrix[a][b]=float(line[3])
	matrix[b][a]=float(line[3])

print "\t".join(all_list)
for i in range(len(matrix)):
	print "\t".join(map(str,matrix[i]))
	