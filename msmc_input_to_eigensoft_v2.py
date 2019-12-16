import argparse
import genutils
from NGS_utils import *
import numpy as np
import pickle
chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = genutils.read_chrom_len(chromLenFile)

def msmc_input_to_eigensoft(file):
	f=open(file,'r')
	f_geno=open('humanAll.geno','w')
	f_snp=open('humanAll.snp','w')
	for line in f:
		line=line.strip().split('\t')
		if line[4]=='-':
			continue
		code=[]
		ref=line[3]
		alt=[]
		chimp=line[4]
		seq=line[9].split(',')[0]
		for i in range(len(seq)):
			if seq[i]!=ref and seq[i] not in alt:
				alt.append(seq[i])
		if len(alt)>1:
			print line
			continue
		for i in range(len(seq)/2):
			if seq[2*i:2*i+2]==ref+ref:
				code.append(0)
			elif seq[2*i:2*i+2]==ref+alt[0] or seq[2*i:2*i+2]==alt[0]+ref:
				code.append(1)
			elif seq[2*i:2*i+2]==alt[0]+alt[0]:
				code.append(2)
			else:
				print 'three allele',line
		if chimp==ref:
			code.append(0)
		elif chimp==alt[0]:
			code.append(1)
		else:
			print 'chimp third allele',line
			continue
		f_geno.write("".join(map(str,code))+'\n')
		f_snp.write('\t%s_%s\t1\t%f\t%s\t%s\t%s\n' %(line[0],line[2],float(line[2])/chromLens[line[0]],line[2],ref,alt[0]))

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
		fields = line.strip().split('\t')
		return tuple(fields)

def read_ped(locus_dict,snp_ref,snp_alt,locus_present,sample_id,pop,pop_ind_num,hg19_to_hg18):
	nuc_map={'1':'A','2':'C','3':'G','4':'T','0':'N'}
	for chr in range(1,23):
		print chr
		f1='1KG_phase3_chr%i.info' %(chr)
		locus_list=[]
		for i in BedIterator(f1):
			locus_list.append(i[0])
		locus_index=[]
		for locus_id in range(len(locus_list)):
			try:
				ref = snp_ref[locus_list[locus_id]]
				alt = snp_alt[locus_list[locus_id]]
			except KeyError:
				try:
					locus=hg19_to_hg18[locus_list[locus_id]]
					ref = snp_ref[locus]
					alt = snp_alt[locus]
				except KeyError:
					print 'could not find locus',locus_list[locus_id],locus
					continue
			locus_index.append(locus_id+6)
		f2='1KG_phase3_chr%i.ped' %(chr)
		snp_matrix=[]
		for i in BedIteratorv2(f2):
			if chr==1:
				sample_id.append(i[1])
				pop[i[1]]=i[0].split('_')[0]
				pop_ind_num[1]+=1
			snp_matrix.append([nuc_map[i[m][0]]+nuc_map[i[m][2]] for m in locus_index])
		snp_matrix=np.matrix(snp_matrix)
		print snp_matrix.shape, len(locus_index)
		for i in range(len(locus_index)):
			locus=locus_list[locus_index[i]-6]
			try:
				ref = snp_ref[locus]
				alt = snp_alt[locus]
			except KeyError:
				try:
					locus=hg19_to_hg18[locus]
					ref = snp_ref[locus]
					alt = snp_alt[locus]
				except KeyError:
					print 'could not find LOCUS',locus
					continue
			alt_allele=[]
			for j in range(pop_ind_num[1]):
				if snp_matrix[j,i][0]!=ref and snp_matrix[j,i][0] not in alt_allele and snp_matrix[j,i][0]!='N':
					alt_allele.append(snp_matrix[j,i][0])
				if snp_matrix[j,i][1]!=ref and snp_matrix[j,i][1] not in alt_allele and snp_matrix[j,i][1]!='N':
					alt_allele.append(snp_matrix[j,i][1])
			if len(alt_allele)>1:
				print 'three allele1',locus,alt_allele,ref
				continue
			if len(alt)>1:
				snp_alt[locus]=alt_allele
				alt=alt_allele
			if len(alt)==1:
				if len(alt_allele)==1:
					if alt[0]!=alt_allele[0]:
						print 'three allele2',locus,alt_allele,alt,ref
						continue
			if len(alt_allele)==0:
				alt=["N"]
			for j in range(pop_ind_num[1]):
				seq=snp_matrix[j,i]
				if seq==ref+ref:
					locus_dict[locus]+="0"
				elif seq==ref+alt[0] or seq==alt[0]+ref:
					locus_dict[locus]+="1"
				elif seq==alt[0]+alt[0]:
					locus_dict[locus]+="2"
				elif seq[0]=='N' or seq[1]=='N':
					locus_dict[locus]+="9"
				else:
					print 'three allele2',locus,locus_index[i],j,ref,alt,seq
					continue
				locus_present[locus][1]+=1
	return locus_dict,locus_present,sample_id,pop,pop_ind_num

def read_ped_v2(locus_dict,snp_ref,snp_alt,locus_present,sample_id,pop,pop_ind_num):
	nuc_map={'1':'A','2':'C','3':'G','4':'T'}
	index=1
	for POP in ['CEU','GIH','LWK','MKK','YRI','CHB','CHD','JPT','TSI']:
		print POP
		index +=1
		f1='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.%s.qc.poly.map' %(POP)
		locus_list=[]
		for i in BedIterator(f1):
			locus_list.append(i[1])
		locus_index=[]
		for locus_id in range(len(locus_list)):
			try:
				ref = snp_ref[locus_list[locus_id]]
				alt = snp_alt[locus_list[locus_id]]
			except KeyError:
				try:
					locus=hg19_to_hg18[locus_list[locus_id]]
					ref = snp_ref[locus]
					alt = snp_alt[locus]
				except KeyError:
#					print 'could not find locus',locus_list[locus_id]
					continue
			locus_index.append(locus_id+6)
		f2='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.%s.qc.poly.ped' %(POP)
		snp_matrix=[]
		for i in BedIteratorv2(f2):
			sample_id.append(i[1])
			pop[i[1]]=POP
			pop_ind_num[index]+=1
			snp_matrix.append([i[m][0]+i[m][2] for m in locus_index])
		snp_matrix=np.matrix(snp_matrix)
		print snp_matrix.shape, len(locus_index)
		for i in range(len(locus_index)):
			locus=locus_list[locus_index[i]-6]
			try:
				ref = snp_ref[locus]
				alt = snp_alt[locus]
			except KeyError:
				try:
					locus=hg19_to_hg18[locus]
					ref = snp_ref[locus]
					alt = snp_alt[locus]
				except KeyError:
					print 'could not find LOCUS',locus
					continue
			alt_allele=[]
			for j in range(pop_ind_num[index]):
				if snp_matrix[j,i][0]!=ref and snp_matrix[j,i][0] not in alt_allele and snp_matrix[j,i][0]!='0':
					alt_allele.append(snp_matrix[j,i][0])
				if snp_matrix[j,i][1]!=ref and snp_matrix[j,i][1] not in alt_allele and snp_matrix[j,i][1]!='0':
					alt_allele.append(snp_matrix[j,i][1])
			if len(alt_allele)>1:
				print 'three allele1',locus,alt_allele,ref
				continue
			if len(alt)>1:
				snp_alt[locus]=alt_allele
				alt=alt_allele
			if len(alt)==1:
				if len(alt_allele)==1:
					if alt[0]!=alt_allele[0]:
						print 'three allele2',locus,alt_allele,alt,ref
						continue
			if len(alt_allele)==0:
				alt=["N"]
			for j in range(pop_ind_num[index]):
				seq=snp_matrix[j,i]
				if seq==ref+ref:
					locus_dict[locus]+="0"
				elif seq==ref+alt[0] or seq==alt[0]+ref:
					locus_dict[locus]+="1"
				elif seq==alt[0]+alt[0]:
					locus_dict[locus]+="2"
				elif seq[0]=='0' or seq[1]=='0':
					locus_dict[locus]+="9"
				else:
					print 'three allele2',locus,locus_index[i],j,ref,alt,seq
					continue
				locus_present[locus][index]+=1
	return locus_dict,locus_present,sample_id,pop,pop_ind_num

def read_store_locus(file):
	locus=[]
	locus_dict={}
	locus_present={}
	snp_ref={}
	snp_alt={}
	hg19_to_hg18={}
	strand={}
	for i in BedIterator(file):
		if i[0]=='chrX':
			break
		if len(i[4])!=1:
			continue
		try:
			m=snp_ref[i[3]]
			continue
		except KeyError:
			a=0
		snp_ref[i[3]]=i[4]
		snp_alt[i[3]]=i[5].split(',')
		locus.append([i[0],i[2],i[3]])
		strand[i[3]]=i[7]
		locus_dict[i[3]]=""
		locus_present[i[3]]=[0,0,0,0,0,0,0,0,0,0,0]
		if i[6]!=i[3]:
			hg19_to_hg18[i[6]]=i[3]
	pop_ind_num=[0,0,0,0,0,0,0,0,0,0,0]
	return locus,locus_dict,snp_ref,snp_alt,strand,locus_present,pop_ind_num,hg19_to_hg18

def read_illumina_chip(file,locus_dict,snp_ref,snp_alt,strand,locus_present,pop_ind_num):
	sample_id=[]
	three_allele=0
	good_allele=0
	map={'A':'T','T':'A','C':'G','G':'C','-':'-'}
	reverse_strand = False
	first = True
	for i in BedIterator(file):
		reverse_strand = False
		if first is True:
			sample_id=list(i)
			pop_ind_num[0]=len(sample_id)
			first = False
			print sample_id
			print pop_ind_num[0]
			continue
		try:
			m=locus_dict[i[0]]
		except KeyError:
			continue
		ref = snp_ref[i[0]]
		alt=[]
		for j in range(1,len(i)):
			a=i[j][0]
			b=i[j][1]
			if a!=ref and a not in alt and a!='-':
				alt.append(a)
			if b!=ref and b not in alt and b!='-':
				alt.append(b)
		if len(alt)>1:
			allele_remove=[]
			if strand[i[0]]=='-':
				for j in range(len(alt)):
					if map[alt[j]]==ref:
						allele_remove.append(alt[j])
				for j in allele_remove:
					alt.remove(j)
				if len(alt)>1:
					print 'three allele',i[0],alt,ref
					three_allele+=1
					continue
				else:
					reverse_strand=True
			else:
				print 'three allele',i[0],alt,ref
				three_allele+=1
				continue
		if len(alt)==1:
			if reverse_strand:
				snp_alt[i[0]]=[map[alt[0]]]
			else:
				snp_alt[i[0]]=[alt[0]]
		elif len(alt)==0:
			alt.append('N')
		drop_allele=False
		for j in range(1,len(i)):
			a=i[j][0]
			b=i[j][1]
			if reverse_strand:
				a=map[a]
				b=map[b]
			if a+b==ref+ref:
				locus_dict[i[0]]+="0"
			elif a+b==ref+snp_alt[i[0]][0] or a+b==snp_alt[i[0]][0]+ref:
				locus_dict[i[0]]+="1"
			elif a+b==snp_alt[i[0]][0]+snp_alt[i[0]][0]:
				locus_dict[i[0]]+="2"
			elif a+b=='--':
				locus_dict[i[0]]+="9"
			else:
				print 'three allele',i[0],reverse_strand,i[j][0],i[j][1],a,b,ref,snp_alt[i[0]]
				drop_allele=True
				break
		if not drop_allele:
			locus_present[i[0]][0]=1
			good_allele+=1
	print 'summary',three_allele,good_allele
	return locus_dict,sample_id,locus_present,snp_alt,pop_ind_num

def read_pop(file):
	pop={}
	for i in BedIteratorv2(file):
		if len(i)==5:
			pop[i[1]]=i[3].split(' ')[0]
	return pop

def write_output(locus,sample_id,pop,locus_dict,locus_present,snp_ref,snp_alt,pop_ind_num):
	f_geno=open('Illumina_chip_1KG_phase3_HapMap_all.geno','w')
	f_snp=open('Illumina_chip_1KG_phase3_HapMap_all.snp','w')
	f_ind=open('Illumina_chip_1KG_phase3_HapMap_all.ind','w')
	for i in range(len(locus)):
		locus_id=locus[i][2]
		chr=locus[i][0][3:]
		present = True
		if locus_present[locus_id][0]==1:
			for j in range(1,11):
				if locus_present[locus_id][j]!=pop_ind_num[j]:
					present = False
		else:
			present = False
		if present:
			print >>f_geno, locus_dict[locus_id]
			print >>f_snp,'\t%s\t%s\t%f\t%s\t%s\t%s' %(locus_id,chr,float(locus[i][1])/chromLens[locus[i][0]],locus[i][1],snp_ref[locus_id],snp_alt[locus_id][0])
	for sample in sample_id:
		print >>f_ind,'\t%s\tF\t%s' %(sample,pop[sample])

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='make input for eigensoft packages')
	parser.add_argument("--file", dest='file',help="file")
	args = parser.parse_args()
	file='/home/jmkidd/kidd-lab-scratch/shiya-projects/HGDP/Illumina_chip/hgdp/HGDP_FinalReport_Forward.txt'
	locus,locus_dict,snp_ref,snp_alt,strand,locus_present,pop_ind_num,hg19_to_hg18=read_store_locus('HGDP_Map_hg19_sort_ref_strand.bed')
	pop=read_pop('HGDP_sample_pop.txt')
	locus_dict,sample_id,locus_present,snp_alt,pop_ind_num=read_illumina_chip(file,locus_dict,snp_ref,snp_alt,strand,locus_present,pop_ind_num)
	locus_dict,locus_present,sample_id,pop,pop_ind_num=read_ped(locus_dict,snp_ref,snp_alt,locus_present,sample_id,pop,pop_ind_num,hg19_to_hg18)
	print 'part2'
	locus_dict,locus_present,sample_id,pop,pop_ind_num=read_ped_v2(locus_dict,snp_ref,snp_alt,locus_present,sample_id,pop,pop_ind_num)
	dbfile=open('locus_dict_pickle_v2','wb')
	pickle.dump(locus,dbfile) 
	pickle.dump(sample_id,dbfile)
	pickle.dump(pop,dbfile)
	pickle.dump(locus_dict,dbfile)
	pickle.dump(locus_present,dbfile)
	pickle.dump(snp_ref,dbfile)
	pickle.dump(snp_alt,dbfile)
	pickle.dump(pop_ind_num,dbfile)
	write_output(locus,sample_id,pop,locus_dict,locus_present,snp_ref,snp_alt,pop_ind_num)