#!/usr/bin/env python
# python phasing_compare.py
# Shiya Song
# 5th Dec 2014
# Comparing phasing statistics in terms of switch error and mean distance between switch errors

import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
import math
import gzip
import pickle
from scipy.interpolate import interp1d

def find_interval(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = (right + left)/2
		if pos < list[midpoint]:
			right = midpoint
		elif pos > list[midpoint]:
			left = midpoint+1
		elif pos == list[midpoint]:
			left = midpoint
			find = midpoint
			break
		if pos<list[midpoint]:
			midpoint=midpoint-1
	return midpoint

def genetic_map(chr):
	POS = []
	MAP = []
	file = "/home/jmkidd/kidd-lab/genomes/hg19/genetic-map/genetic_map_GRCh37_CHROM.txt"
	file = file.replace("CHROM",chr)
	f = open(file,"r")
	for line in f:
		if line[:5]=="Chrom":
			continue
		line=line.strip().split("\t")
		pos = int(line[1])
		map = float(line[3])
		POS.append(pos)
		MAP.append(map)
	POS=np.array(POS)
	MAP=np.array(MAP)
	interpolate=interp1d(POS,MAP)
	return interpolate,POS,MAP	

def read_snp(SNP,vcf_file,bed_file,ind):		# read in hapmap phasing after liftOver
	for i in VcfIterator(vcf_file):
		try:
			info = SNP["phase"][i[1]]
		except KeyError:
			continue
		ref=SNP["ref"][i[1]]
		alt=SNP["alt"][i[1]]
		if ref!=i[2][0]:
			print 'ref does not match',i,ref
			continue
		if len(alt)==1:
			if i[2][1]!=alt[0]:
				print 'third allele',i,ref,alt
				continue
		else:
			SNP["alt"][i[1]]=[i[2][1]]
			alt=i[2][1]
		ref_alt=i[2][0]+i[2][1]
		if i[4] is True or i[3][0]==i[3][1]:
			SNP["phase"][i[1]].append(i[3][0])
			SNP["phase"][i[1]].append(i[3][1])
		else:
			SNP["phase"][i[1]].append(9)
			SNP["phase"][i[1]].append(9)
	list = SNP.index
	left_pos=[]
	right_pos=[]
	for i in BedIterator(bed_file):
		left_pos.append(int(i[1]))
		right_pos.append(int(i[2])-1)
	for i in range(len(list)):
		pos = list[i]
		if len(SNP["phase"][list[i]])<2*(len(ind)-1):
			print len(SNP["phase"][list[i]]),2*(len(ind)-1)
			continue
		elif len(SNP["phase"][list[i]])==2*len(ind):
			continue
		else:
			if 9 in SNP["phase"][list[i]]:
				continue
			find=find_interval(pos,left_pos)
			try:
				assert pos>=left_pos[find]
			except AssertionError:
#				print pos,find,left_pos[find]
				continue
			if find<len(left_pos)-1:
				assert pos<=left_pos[find+1]
			if pos<=right_pos[find]:
				print 'Find',pos,find,left_pos[find],right_pos[find]
				SNP["phase"][pos].append(0)
				SNP["phase"][pos].append(0)
	return SNP

def read_sample(file):
	sample_list1=[]
	sample_list2=[]
	f = open(file,'r')
	index = 0
	for line in f:
		line = line.rstrip().split(' ')
		if line[0]=='sample':
			continue
		if line[1]=='YRI':
			sample_list1.append(index)
		if line[1]=='CEU':
			sample_list2.append(index)
		index+=1
	return sample_list1,sample_list2
	
def read_samplename(file):
	sample_list1=[]
	sample_list2=[]
	f = open(file,'r')
	index = 0
	for line in f:
		line = line.rstrip().split(' ')
		if line[0]=='sample':
			continue
		if line[1]=='YRI':
			sample_list1.append(line[0])
		if line[1]=='CEU':
			sample_list2.append(line[0])
		index+=1
	return sample_list1,sample_list2

def read_reference_SNP(SNP,hap_file,YRI_list,CEU_list,pos_index,index_to_pos):
	index = 0
	list = SNP.index
	for i in BedIterator(hap_file):
		if pos_index[index]==0:
			index+=1
			continue
		pos = index_to_pos[index]
		for j in YRI_list:
			SNP["phase"][pos].append(i[j*2])
			SNP["phase"][pos].append(i[j*2+1])
		for j in CEU_list:
			SNP["phase"][pos].append(i[j*2])
			SNP["phase"][pos].append(i[j*2+1])
		index+=1
		if index%10000==0:
			print index,'done'
	return SNP

def write_output(SNP,ind,chr,interpolate,POS,MAP):
	f_geno=open('alleles.'+chr+'.txt','w')
	f_snp=open('positions.'+chr+'.detail.txt','w')
	f_snp2=open('positions.'+chr+'.txt','w')
	list = SNP.index
	for i in range(len(list)):
		pos = list[i]
		ref_alt=SNP["ref"][pos]+SNP["alt"][pos][0]
		try:
			cm = interpolate(pos)
		except ValueError:
			if pos<=POS[0]:
				cm = min(MAP)
			else:
				cm = max(MAP)
		if len(SNP["phase"][pos])==2*len(ind):
			if 9 not in SNP["phase"][pos]:
				print >>f_geno,''.join(map(str,SNP["phase"][pos]))
				print >>f_snp,pos,cm,ref_alt
				print >>f_snp2,cm

def write_output_v2(SNP,ind,chr):
	f_snp=open('mapfile.'+chr+'.txt','w')
	f1 = open('admixed.'+chr+'.txt','w')
	f2 = open('yri.'+chr+'.txt','w')
	f3 = open('ceu.'+chr+'.txt','w')
	list = SNP.index
	het =0
	for i in range(len(list)):
		pos = list[i]
		ref_alt=SNP["ref"][pos]+SNP["alt"][pos][0]
		if len(SNP["phase"][pos])==2*len(ind):
			if 9 not in SNP["phase"][pos]:
				a=[]
				b=[]
				c=[]
				for j in range(len(ind)):
					if ind[j] ==0:
						if int(SNP["phase"][pos][2*j])!=int(SNP["phase"][pos][2*j+1]):
#							print 'het',pos,int(SNP["phase"][pos][2*j]),int(SNP["phase"][pos][2*j+1])
							het+=1
						a.append(ref_alt[int(SNP["phase"][pos][2*j])])
						a.append(ref_alt[int(SNP["phase"][pos][2*j+1])])
					if ind[j]==1:
						b.append(ref_alt[int(SNP["phase"][pos][2*j])])
						b.append(ref_alt[int(SNP["phase"][pos][2*j+1])])
					if ind[j]==2:
						c.append(ref_alt[int(SNP["phase"][pos][2*j])])
						c.append(ref_alt[int(SNP["phase"][pos][2*j+1])])
				print >>f1,'M %s %s' %(SNP["rs"][pos],' '.join(a))
				print >>f2,'M %s %s' %(SNP["rs"][pos],' '.join(b))
				print >>f3,'M %s %s' %(SNP["rs"][pos],' '.join(c))
				print >>f_snp,'%s %s %s %s' %(chr[3:],SNP["rs"][pos],0,pos)
	print 'tot het',het

def read_hapmap_ped(chr):
	f_map=open('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/hapmap3_r3_b36_fwd.consensus.qc.poly.map','r')
	f_hg19=open('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/hapmap3_r3_b36_fwd.consensus.qc.poly.hg19.dbSNP.bed','r')
	rs_id = []
	hg18_to_hg19={}
	hg19_ref={}
	hg19_alt={}
	for line in f_hg19:
		line = line.rstrip().split('\t')
		if line[0]==chr:
			hg18_to_hg19[int(line[6])]=int(line[2])
			hg19_ref[int(line[2])]=line[4]
			hg19_alt[int(line[2])]=line[5].split(',')
		elif int(line[0][3:])> int(chr[3:]):
			break
	index = 0
	pos_hg19=[]
	pos_chr=[]
	index_list=[]
	for line in f_map:
		line = line.rstrip().split('\t')
		hg19_pos = None
		if int(line[0])> int(chr[3:]):
			break
		elif line[0]==chr[3:]:
			try:
				hg19_pos=hg18_to_hg19[int(line[3])]
			except KeyError:
				index +=1
				continue
			rs_id.append(line[1])
			pos_hg19.append(hg19_pos)
			index_list.append(index+6)
			index +=1
		else:
			index +=1
	print len(pos_hg19),len(index_list)
	CEU_list=[]
	YRI_list=[]
	MKK_list=[]
	for i in BedIterator('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/CEU_list.txt'):
		if i[1]=='0' and i[2]=='0':
			CEU_list.append(i[0])
		else:
			print i
	for i in BedIterator('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/YRI_list.txt'):
		if i[1]=='0' and i[2]=='0':
			YRI_list.append(i[0])
		else:
			print i
	for i in BedIterator('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/MKK_list.txt'):
		if i[1]=='0' and i[2]=='0':
			if i[0]!='NA21302' or 'NA21301' or 'NA21301':
				MKK_list.append(i[0])
		else:
			print i
	MKK_list.append('NA21302')
	print len(YRI_list),len(CEU_list),len(MKK_list)
	sample_list=[]
	pop_list=[]
	f_ped=open('/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/PLINK/hapmap3_r3_b36_fwd.consensus.qc.poly.ped','r')
	snp_matrix=[]
	for line in f_ped:
		line = line.rstrip().split('\t')
		if line[1] in YRI_list:
			sample_list.append(line[1])
			pop_list.append(1)
		elif line[1] in CEU_list:
			sample_list.append(line[1])
			pop_list.append(2)
		elif line[1] in MKK_list:
			sample_list.append(line[1])
			pop_list.append(0)
		else:
			continue
		snp_matrix.append([line[m][0]+line[m][2] for m in index_list])
	snp_matrix=np.matrix(snp_matrix)
	print snp_matrix.shape,len(sample_list),len(pop_list),len(index_list)
	haplotype = []  # This stores the haplotype
	ref_list=[]
	alt_list=[]
	for i in range(len(index_list)):
		pos = pos_hg19[i]
		ref = hg19_ref[pos]
		alt = hg19_alt[pos]
		alt_allele=[]
		for j in range(snp_matrix.shape[0]):
			if snp_matrix[j,i][0]!=ref and snp_matrix[j,i][0] not in alt_allele and snp_matrix[j,i][0] in ['A','T','C','G']:
				alt_allele.append(snp_matrix[j,i][0])
			if snp_matrix[j,i][1]!=ref and snp_matrix[j,i][1] not in alt_allele and snp_matrix[j,i][1] in ['A','T','C','G']:
				alt_allele.append(snp_matrix[j,i][1])
		if len(alt_allele)>1:
			print 'three allele1',pos,ref,alt,alt_allele
			haplotype.append([])
			ref_list.append('N')
			alt_list.append(['N'])
			continue
		if len(alt_allele)==1:
			hg19_alt[pos]=alt_allele
		alt = hg19_alt[pos]
		m = []
		for j in range(snp_matrix.shape[0]):
			seq=snp_matrix[j,i]
			if seq[0]==ref:
				m.append(0)
			elif seq[0]==alt[0]:
				m.append(1)
			else:
				m.append(9)
			if seq[1]==ref:
				m.append(0)
			elif seq[1]==alt[0]:
				m.append(1)
			else:
				m.append(9)
		haplotype.append(m)
		ref_list.append(ref)
		alt_list.append(alt)
	print len(haplotype),len(index_list)
	snp = {}
	snp["rs"] = pd.Series(rs_id,index=pos_hg19)
	snp["ref"] = pd.Series(ref_list,index=pos_hg19)
	snp["alt"] = pd.Series(alt_list,index=pos_hg19)
	snp["phase"] = pd.Series(haplotype,index=pos_hg19)
	snp = pd.DataFrame(snp)
	return snp,sample_list,pop_list

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--name",default='NA21302',help="Input VCF files")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--ref_dir", dest='ref_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/',help="directory for 1KG reference file")
	parser.add_argument("--chr", dest='chr',help="directory for 1KG reference file")

	args = parser.parse_args()
	SNP,sample_list,pop_list=read_hapmap_ped(args.chr)

#	sample_list.append('NA21302')
#	pop_list.append(0)

#	vcf_file = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,args.chr)
#	bed_file = '%s%s/gVCF_calls/%s_%s.mask.bed.gz' %(args.wgs_dir,args.name,args.name,args.chr)
#	SNP = read_snp(SNP,vcf_file,bed_file,sample_list)

	dbfile=open('hapmap_%s_pickle' %(args.chr),'wb')
	pickle.dump(SNP,dbfile) 
	interpolate,POS,MAP=genetic_map(args.chr)
	write_output(SNP,sample_list,args.chr,interpolate,POS,MAP)  	# for RFMix
	write_output_v2(SNP,pop_list,args.chr)			# for PCAdmix
#	f=open('classes.txt','w')
#	print >>f,' '.join(map(str,pop_list))
#	f=open('samples.txt','w')
#	print >>f,' '.join(map(str,sample_list))




		
