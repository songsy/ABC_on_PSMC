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


breed_list=list(breed_list)
village_dog_list=list(village_dog_list)
ancient_list=list(ancient_list)
grey_wolf_list=list(grey_wolf_list)
print grey_wolf_list
dog_list=ancient_list+breed_list+village_dog_list
'''
for dog in dog_list:
	f=open('Jackel_%s_wolfs_dstats.txt' %(dog),'r')
	matrix=[[0 for i in range(len(grey_wolf_list))] for j in range(len(grey_wolf_list))]
#print len(matrix[0]),len(matrix),len(all_list)
	for line in f:
		line = line.rstrip().split('\t')
		a=grey_wolf_list.index(line[0])
		b=grey_wolf_list.index(line[1])
		matrix[a][b]=float(line[5])
		matrix[b][a]=-float(line[5])
	fout=open('Jackel_%s_wolfs_dstats_Zscore_matrix.txt' %(dog),'w')
	print >>fout,"\t".join(grey_wolf_list)
	for i in range(len(matrix)):
		print >>fout,"\t".join(map(str,matrix[i]))
'''
for wolf in grey_wolf_list:
	f=open('Jackel_%s_dogs_dstats.txt' %(wolf),'r')
	matrix=[[0 for i in range(len(dog_list))] for j in range(len(dog_list))]
	#print len(matrix[0]),len(matrix),len(all_list)
	for line in f:
		line = line.rstrip().split('\t')
		a=dog_list.index(line[0])
		b=dog_list.index(line[1])
		matrix[a][b]=float(line[5])
		matrix[b][a]=-float(line[5])
	fout=open('Jackel_%s_dogs_dstats_Zscore_matrix.txt' %(wolf),'w')
	print >>fout,"\t".join(dog_list)
	for i in range(len(matrix)):
		print >>fout,"\t".join(map(str,matrix[i]))
