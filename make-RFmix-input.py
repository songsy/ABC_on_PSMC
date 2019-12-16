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

def create_snp_list(file1):
	snp_pos = []
	snp_geno = []
	snp_phase = []
	f=gzip.open(legend_file)
	pos_index =[]
	index_to_pos=[]
	rs_id = []
	for line in f:
		line = line.rstrip().split(' ')
		if line[0]=='id':
			continue
		pos=int(line[1])
		if line[4]!='SNP':
			pos_index.append(0)
			index_to_pos.append(pos)   # map row number to pos
			continue
		rs_id.append(line[0])
		index_to_pos.append(pos)
		pos_index.append(1)
		snp_pos.append(pos)				# store the position of the snp, used for fast search 
		snp_geno.append(line[2]+line[3])		# store the ref and alt allele of the snp
		snp_phase.append([])			# store the phasing information, like [0,1]
	snp = {}
	snp["rs"] = pd.Series(rs_id,index=snp_pos)
	snp["geno"] = pd.Series(snp_geno,index=snp_pos)
	snp["phase"] = pd.Series(snp_phase,index=snp_pos)
	snp = pd.DataFrame(snp)
	return snp,pos_index,index_to_pos

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
	pos_list = SNP.index
	pos_list = map(int,pos_list)
	for i in VcfIterator(vcf_file):
		find_real = find_pos(int(i[1]),pos_list)
		if find_real >=0:
			info = SNP["geno"][i[1]]
		else:
			continue
		'''
		try:
			info = SNP["geno"][i[1]]
		except KeyError:
			continue
		'''
		ref_alt=i[2][0]+i[2][1]
		if SNP["geno"][i[1]]!=ref_alt:
#			print i,SNP["geno"][i[1]]
			SNP["phase"][i[1]].append(2)
			SNP["phase"][i[1]].append(2)
		else:
			if i[4] is True or i[3][0]==i[3][1]:
#				print 'het',i[1],i[3]
				SNP["phase"][i[1]].append(i[3][0])
				SNP["phase"][i[1]].append(i[3][1])
			else:
				SNP["phase"][i[1]].append(2)
				SNP["phase"][i[1]].append(2)
	list = SNP.index
	left_pos=[]
	right_pos=[]
	for i in BedIterator(bed_file):
		left_pos.append(int(i[1]))
		right_pos.append(int(i[2])-1)
	for i in range(len(list)):
		pos = list[i]
		if len(SNP["phase"][list[i]])<2*(len(ind)-1):
			continue
		elif len(SNP["phase"][list[i]])==2*len(ind):
			continue
		else:
			if 2 in SNP["phase"][list[i]]:
				continue
			find=find_interval(pos,left_pos)
			try:
				assert pos>=left_pos[find]
			except AssertionError:
				print pos,find,left_pos[find]
				continue
			if find<len(left_pos)-1:
				assert pos<=left_pos[find+1]
			if pos<=right_pos[find]:
#				print 'Find',pos,find,left_pos[find],right_pos[find]
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
		try:
			cm = interpolate(pos)
		except ValueError:
			if pos<=POS[0]:
				cm = min(MAP)
			else:
				cm = max(MAP)
		if len(SNP["phase"][pos])==2*len(ind):
			if 2 not in SNP["phase"][pos]:
				print >>f_geno,''.join(map(str,SNP["phase"][pos]))
				print >>f_snp,pos,cm,SNP["geno"][pos]
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
		if len(SNP["phase"][pos])==2*len(ind):
			if 2 not in SNP["phase"][pos]:
				a=[]
				b=[]
				c=[]
				for j in range(len(ind)):
					if ind[j] ==0:
						if int(SNP["phase"][pos][2*j])!=int(SNP["phase"][pos][2*j+1]):
#							print 'het',pos,int(SNP["phase"][pos][2*j]),int(SNP["phase"][pos][2*j+1])
							het+=1
						a.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j])])
						a.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j+1])])
					if ind[j]==1:
						b.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j])])
						b.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j+1])])
					if ind[j]==2:
						c.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j])])
						c.append(SNP["geno"][pos][int(SNP["phase"][pos][2*j+1])])
				print >>f1,'M %s %s' %(SNP["rs"][pos],' '.join(a))
				print >>f2,'M %s %s' %(SNP["rs"][pos],' '.join(b))
				print >>f3,'M %s %s' %(SNP["rs"][pos],' '.join(c))
				print >>f_snp,'%s %s %s %s' %(chr[3:],SNP["rs"][pos],0,pos)
	print 'tot het',het

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--name",default='NA21302',help="Input VCF files")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--ref_dir", dest='ref_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/',help="directory for 1KG reference file")
	parser.add_argument("--chr", dest='chr',help="directory for 1KG reference file")

	args = parser.parse_args()

	YRI_list,CEU_list=read_sample(args.ref_dir+'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample')
	ind=[0]
	chr = 'chr'+args.chr
	legend_file=args.ref_dir+'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(chr)
	SNP,pos_index,index_to_pos = create_snp_list(legend_file)

	vcf_file = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
	bed_file = '%s%s/gVCF_calls/%s_%s.mask.bed.gz' %(args.wgs_dir,args.name,args.name,chr)
	SNP = read_snp(SNP,vcf_file,bed_file,ind)

	for sample in ['NA21732','NA21733','NA21737','NA21767']:
		print sample
		ind.append(0)
		vcf_file = '/share/jmkidd/songsy/complete-genomics/MKK/%s.%s.phased.vcf.gz' %(sample,chr)
		bed_file = '/share/jmkidd/songsy/complete-genomics/MKK/%s_%s.mask.bed.gz' %(sample,chr)
		SNP = read_snp(SNP,vcf_file,bed_file,ind)

	hap_file = args.ref_dir+'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz' %(chr)
	SNP=read_reference_SNP(SNP,hap_file,YRI_list,CEU_list,pos_index,index_to_pos)
	for i in range(len(YRI_list)):
		ind.append(1)
	for i in range(len(CEU_list)):
		ind.append(2)
	dbfile=open(chr+'_pickle','wb')
	pickle.dump(SNP,dbfile) 
	interpolate,POS,MAP=genetic_map(chr)
	write_output(SNP,ind,chr,interpolate,POS,MAP)
	write_output_v2(SNP,ind,chr)
#	f=open('classes.txt','w')
#	print >>f,' '.join(map(str,ind))

		
