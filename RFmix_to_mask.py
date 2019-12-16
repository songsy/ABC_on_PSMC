
import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def RFMix_to_mask(chr,RFMix_file,snp_file,f):
	ancestry=[]
	for i in BedIterator(RFMix_file):
		ancestry.append(int(i[0]))
	index = 0
	start_pos = None
	end_pos = None
	for i in BedIterator(snp_file):
		pos = int(i[0])
		if ancestry[index]==2:
			if start_pos is None:
				start_pos = pos
			end_pos = pos
		if ancestry[index]==1:
			if start_pos is not None and end_pos is not None:
				print >>f,'%s\t%s\t%s' %(chr,start_pos-1,end_pos)   # in bed format
				start_pos = None
				end_pos = None
		index +=1
	if ancestry[-1]==2:
		if start_pos is not None and end_pos is not None:
			print >>f,'%s\t%s\t%s' %(chr,start_pos-1,end_pos)	 # in bed format
			start_pos = None
			end_pos = None

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

def mask_to_karyogram(file):
	chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
	chromLens = genutils.read_chrom_len(chromLenFile)
	done = []
	for i in BedIterator(file):
		chr = 'chr'+i[0]
		if chr not in done:
			if len(done)!=0:
				last_chr=done[-1]
				try:
					cm=interpolate(float(last))
					print '%s\t%s\t%i\t%s\t%f\t%f' %(last_chr[3:],last,POS[-1],'AFR',cm,max(MAP))
				except ValueError:
					if float(last)<=POS[0]:
						cm = min(MAP)
					else:
						cm = max(MAP)
			interpolate,POS,MAP=genetic_map(chr)
			try:
				cm=interpolate(float(i[1]))
				print '%s\t%i\t%s\t%s\t%f\t%f' %(i[0],POS[0],i[1],'AFR',min(MAP),cm)
			except ValueError:
				if float(i[1])<=POS[0]:
					cm = min(MAP)
				else:
					cm = max(MAP)
			try:
				cm2=interpolate(float(i[2]))
			except ValueError:
				if float(i[2])<=POS[0]:
					cm2 = min(MAP)
				else:
					cm2 = max(MAP)
			print '%s\t%s\t%s\t%s\t%f\t%f' %(i[0],i[1],i[2],'EUR',cm,cm2)
			done.append(chr)
			last = i[2]
		else:
			try:
				cm1=interpolate(float(i[1]))
			except ValueError:
				if float(i[1])<=POS[0]:
					cm1 = min(MAP)
				else:
					cm1 = max(MAP)
			print '%s\t%s\t%s\t%s\t%f\t%f' %(i[0],last,i[1],'AFR',interpolate(float(last)),cm1)
			try:
				cm2=interpolate(float(i[2]))
			except ValueError:
				if float(i[2])<=POS[0]:
					cm2 = min(MAP)
				else:
					cm2 = max(MAP)
			print '%s\t%s\t%s\t%s\t%f\t%f' %(i[0],i[1],i[2],'EUR',cm1,cm2)
			last = i[2]
	last_chr=done[-1]
	print '%s\t%s\t%f\t%s\t%f\t%f' %(last_chr[3:],last,POS[-1],'AFR',interpolate(float(last)),max(MAP))

		
def PCA_input_mask(snp_file,geno_file,mask_file):
	mask=[]
	new_snp_file=open('Illumina_chip_1KG_phase3_HapMap_Af_Eu_NA21302_0.2cm_G200_masked.snp','w')
	new_geno_file=open('Illumina_chip_1KG_phase3_HapMap_Af_Eu_NA21302_0.2cm_G200_masked.geno','w')
	for i in BedIterator(mask_file):
		mask.append(int(i[9]))
	index = 0
	f = open(snp_file,'r')
	for line in f:
		if mask[index]==0:
			new_snp_file.write(line)
		index+=1
	index = 0
	f = open(geno_file,'r')
	for line in f:
		if mask[index]==0:
			new_geno_file.write(line)
		index+=1
		
if __name__=="__main__":
	'''
	out_file = '/mnt/EXT/Kidd-scratch/shiya-projects/admixture/MKK/1000G/RFMix/NA21302_0.2cm_G88_0_mask.bed'
	f=open(out_file,'w')
	for i in range(1,23):
		RFMix_file ='/mnt/EXT/Kidd-scratch/shiya-projects/admixture/MKK/1000G/RFMix/MKK.chr%s.0.2cm.G88.0.Viterbi.txt' %(i)
		snp_file = '/mnt/EXT/Kidd-scratch/shiya-projects/admixture/MKK/1000G/positions.chr%s.detail.txt' %(i)
		RFMix_to_mask(i,RFMix_file,snp_file,f)
		print i

	snp_file = 'Illumina_chip_1KG_phase3_HapMap.snp'
	geno_file = 'Illumina_chip_1KG_phase3_HapMap_Af_Eu.geno'
	mask_file = 'Illumina_chip_1KG_phase3_HapMap.snp.NA21302.0.2cm.G200.masked.bed'
	PCA_input_mask(snp_file,geno_file,mask_file)
	'''
#	f.close()
	mask_to_karyogram('NA21302_0.2cm_G88_0_mask.bed')