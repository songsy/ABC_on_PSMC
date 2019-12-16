f=open('../plink2/all101_CTC_HXH_Freedman_SampleInfo_label.txt','r')
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

dog_list=ancient_list+breed_list+village_dog_list


def get_order(wolf_compare):
	value=list(set([wolf_compare[i] for i in wolf_compare.keys()]))
	value.sort(reverse=True)
	order=[]
	for i in value:
		m=[]
		for j in wolf_compare.keys():
			if wolf_compare[j]==i:
				m.append(j)
		order.append(m)
	return order

def get_order2(order,compare_result):
	i=0
	order2=[]
	while i<len(order):
		if len(order[i])==1:
			order2.append([order[i]])
			m=order[i].pop()
		i+=1
def get_closest_wolf():
	wolf_list=grey_wolf_list
	wolf_list.sort()
	compare_result={}
	Zscore_list=[]
	for dog in dog_list:
		f=open('Jackel_%s_wolfs_dstats.txt' %(dog),'r')
		wolf_order=[]
		wolf_compare={i:0 for i in wolf_list}
		for line in f:
			line = line.strip().split('\t')
			Zscore=float(line[5])
			wolf1=line[0]
			wolf2=line[1]
			pair1=[wolf1,wolf2]
			pair2=[wolf2,wolf1]
			if Zscore<-3:  # favor wolf2
				compare_result["_".join(pair1)]=-1
				compare_result["_".join(pair2)]=1
				Zscore_list.append([abs(Zscore),wolf2,wolf1])
				wolf_compare[wolf2]+=1
			elif Zscore>3: # favor wolf1
				compare_result["_".join(pair1)]=1
				compare_result["_".join(pair2)]=-1
				Zscore_list.append([abs(Zscore),wolf1,wolf2])
				wolf_compare[wolf1]+=1
			else:
				compare_result["_".join(pair1)]=0
				compare_result["_".join(pair2)]=0
				Zscore_list.append([abs(Zscore),wolf1,wolf2])
		Zscore_list.sort( key=lambda w: int(w[0]),reverse=True)  # sort Zscore from large to small
#		for i in sorted(compare_result.keys()):
#			print i,compare_result[i]
#		print Zscore_list
		wolf_ordered=[]
		for i in Zscore_list:
			if wolf_ordered==len(wolf_list):
				break
			wolf=i[1]
			if wolf in wolf_ordered:
				continue
			if len(wolf_order)==0:
				wolf_order.append([wolf])
				wolf_order.append([i[2]])
				wolf_ordered.append(wolf)
				wolf_ordered.append(i[2])
			else:
				j = 0
				while j<len(wolf_order):
					wolf_to_be_compare=wolf_order[j]
#					print wolf,'compare:',j,wolf_order,wolf_to_be_compare
					same = False
					for m in wolf_to_be_compare:
						pair=[wolf,m]
						if compare_result["_".join(pair)]>0:
#							print 'insert',j,wolf_order,wolf_to_be_compare
							wolf_order.insert(j,[wolf])
							wolf_ordered.append(wolf)
							same = False
							break
						elif compare_result["_".join(pair)]==0:
							same=True
					if same:
						wolf_order[j].append(wolf)
						wolf_ordered.append(wolf)
						break
					if wolf in wolf_ordered:
						break
					j+=1
				if wolf not in wolf_ordered:
					wolf_order.append([wolf])
					wolf_ordered.append(wolf)
#			print i,wolf_order
		order=get_order(wolf_compare)
#		order2=get_order2(order,compare_result)
#		print dog,wolf_order,order
#		";".join([",".join(i) for i in wolf_order])
#		print "%s\t%s\t%s\t%s\t%s" %(dog,",".join(wolf_order[0]),",".join(order[0]),",".join(wolf_order[1]),",".join(order[1]))
		print "%s\t%s\t%s" %(dog,";".join([",".join(i) for i in wolf_order]),";".join([",".join(i) for i in order]))
if __name__=="__main__":		
	get_closest_wolf()