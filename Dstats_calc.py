
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

def get_closest_wolf():
	wolf_list=['chw', 'glw', 'IBW', 'inw', 'irw', 'ita', 'mxa', 'mxb', 'ptw', 'spw', 'ysa', 'ysb', 'ysc','ChineseWolf', 'CroatianWolf', 'IsraeliWolf']
	wolf_list.sort()
	dog_list=['1233', '1735', '1756', '2972', '4669', '6610', '8542', '10442', '13131', '14529', '14566', '15630', '16145', '20576', '21270', '23356', 'ali', 'BA19', 'box', 'cec', 'CTCdog', 'dlr', 'Dog01', 'Dog02', 'Dog03', 'Dog04', 'Dog05', 'Dog06', 'Dog07', 'Dog08', 'Dog09', 'Dog10', 'Dog11', 'Dog12', 'Dog13', 'Dog14', 'Dog15', 'EG44', 'EG49', 'HR85', 'HR93', 'ID60', 'ID91', 'ID125', 'ID137', 'ID165', 'ID168', 'IN18', 'IN23', 'IN29', 'jcc', 'LB74', 'LB79', 'LB85', 'mag', 'mba', 'NA8', 'NA63', 'NA89', 'NGSD1', 'NGSD2', 'NGSD3', 'osp', 'PG84', 'PG115', 'PG122', 'PT49', 'PT61', 'PT71', 'QA5', 'QA27', 'sai', 'tch', 'tpw', 'TW04', 'VN4', 'VN21', 'VN37', 'VN42', 'VN59', 'VN76', 'wpw', 'zoey', 'Basenji', 'Dingo']
	compare_result={}
	Zscore_list=[]
	for dog in dog_list:
		if dog=="CTCdog":
			dog="CTC"
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
#		print dog,wolf_order,order
		print "%s\t%s\t%s\t%s\t%s" %(dog,",".join(wolf_order[0]),",".join(order[0]),",".join(wolf_order[1]),",".join(order[1]))

if __name__=="__main__":		
	get_closest_wolf()