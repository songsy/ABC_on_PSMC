
import argparse
class BedIteratorv2:
	def __init__(self, filename):
		if filename[-2:]=='gz':
			self.file = gzip.open(filename, "r")
		else:
			self.file = open(filename, "r")
	def __iter__(self):
		return self
	def __next__(self):
		return self.next()
	def next(self):
		line = next(self.file)
		fields = line.rstrip().split('\t')
		return tuple(fields)

def subset(file_geno,file_ind,regionToPop):
	index_list =[]
	index=0
	f_ind=open(file_ind.replace('.ind','_Af_NA21302.ind'),'w')
	for i in BedIteratorv2(file_ind):
		if i[3] in regionToPop.keys():
			index_list.append(index)
			print >>f_ind,'\t'.join(i)
		elif i[3] in ['YRI','LWK','ESN','GWD','MSL']:
			index_list.append(index)
			print >>f_ind,'\t'.join(i)
		elif i[1]=='NA21302':
			index_list.append(index)
			print >>f_ind,'\t'.join(i)
		'''
		if i[3] in ['YRI','MKK','LWK','GIH','ESN','GWD','MSL','CEU','CHB','CHD','JPT','TSI']:
			index_list.append(index)
			print >>f_ind,'\t'.join(i)

		elif regionToPop[i[3]] in ['Afr','Eur','YRI','MKK','LWK','GIH','ESN','GWD','MSL','CEU','San','Mbuti']:
			index_list.append(index)
			print >>f_ind,'\t'.join(i)
		'''
		index+=1
	f_geno=open(file_geno.replace('.geno','_Af_NA21302.geno'),'w')
	for i in BedIteratorv2(file_geno):
		m=list(i[0])
		m_new=[m[j] for j in index_list]
		print >>f_geno,''.join(m_new)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='make input for eigensoft packages')
	parser.add_argument("--ind", dest='file_ind',help="ind")
	parser.add_argument("--geno", dest='file_geno',help="ind")
	args = parser.parse_args()
#	regionToPop={'Bantu':'Afr','Mandenka':'Afr','Yoruba':'YRI','San':'San','Mbuti':'Mbuti','Biaka':'Afr','Mozabite':'Afr','Orcadian':'Eur','Adygei':'Eur','Russian':'Eur','Basque':'Eur','French':'Eur','North':'Eur','Sardinian':'Eur','Tuscan':'Eur','Bedouin':'AsianW','Druze':'AsianW','Palestinian':'AsianW','Balochi':'AsianC','Brahui':'AsianC','Makrani':'AsianC','Sindhi':'AsianC','Pathan':'AsianC','Burusho':'AsianC','Hazara':'AsianC','Uygur':'AsianC','Kalash':'AsianC','Han':'AsianE','Dai':'AsianE','Daur':'AsianE','Hezhen':'AsianE','Lahu':'AsianE','Miaozu':'AsianE','Oroqen':'AsianE','She':'AsianE','Tujia':'AsianE','Tu':'AsianE','Xibo':'AsianE','Yizu':'AsianE','Mongola':'AsianE','Naxi':'AsianE','Cambodian':'AsianE','Japanese':'AsianE','Yakut':'AsianE','NAN':'Ocean','Papuan':'Ocean','Karitiana':'NatAm','Surui':'NatAm','Colombian':'NatAm','Maya':'NatAm','Pima':'NatAm'}
#	PopToColor = {'Eur':'orange','YRI':'darkgreen','MKK':'greenyellow','LWK':'green','Afr':'brown','ESN':'Olive','GWD':'Teal','MSL':'SpringGreen','CEU':'Tomato','GIH':'Yellow','San':'Purple','Mbuti':'Violet','AsianW':'Cyan','AsianC':'DodgerBlue','AsianE':'MediumBlue','Ocean':'Pink','NatAm':'Gray'}
	regionToPop={'Bantu':'Afr','Mandenka':'Afr','Yoruba':'YRI','San':'San','Mbuti':'Mbuti','Biaka':'Afr','Mozabite':'Afr'}
	subset(args.file_geno,args.file_ind,regionToPop)
	