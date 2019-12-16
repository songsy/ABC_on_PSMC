#!/usr/bin/env python
# Shiya Song
# 3rd July 2014
# generate_psmcfa.py

import sys
import signal

def compare_seq(seq1,seq2,h,tot):
	try:
		assert len(seq1)==len(seq2)
	except AssertionError:
		print seq1
		print seq2
	non_N = 0
	het = 0
	symbol = ""
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			non_N +=1
			tot += 1
			if seq1[i]!=seq2[i]:
				het +=1
				h +=1
	if non_N >=10:
		if het >=1:
			symbol = "K"
		else:
			symbol = "T"
	else:
		symbol="N"
	return symbol,h,tot
	
def fa_to_psmcfa(file1,file2,outfile):
	f_out=open(outfile,"w")
	f1=open(file1,"r")
	f2=open(file2,"r")
	symbol = ""
	head = False
	h = 0
	tot = 0
	while True:
		line1=f1.readline().strip()
		line2=f2.readline().strip()
		if line1=="":
			break
		if line1[0]==">":
			if symbol !="":
				print >>f_out,symbol
			assert line1==line2
			print >>f_out,line1
			symbol = ""
		else:
			seq1=line1
			seq2=line2
			line1=f1.readline().strip()
			line2=f2.readline().strip()
			if line1=="":
				s,h,tot= compare_seq(seq1,seq2,h,tot)
				symbol += s
				print >>f_out,symbol
				break
			if line1[0]==">":
				s,h,tot = compare_seq(seq1,seq2,h,tot)
				symbol += s
				print >>f_out,symbol
				assert line1==line2
				print >>f_out,line1
				symbol = ""
			else:
				seq1+=line1
				seq2+=line2
				s,h,tot = compare_seq(seq1,seq2,h,tot)
				symbol += s
				if len(symbol)==60:
					print >>f_out,symbol
					symbol = ""

def fq_to_psmcfa(file1,file2,outfile):
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	f_out=open(outfile,"w")
	f1=os.popen('gunzip -c ' + file1, 'r')
	f2=os.popen('gunzip -c ' + file2, 'r')
	symbol = ""
	head = False
	h = 0
	tot = 0
	while True:
		line1=f1.readline().strip()
		line2=f2.readline().strip()
		if line1=="":
			break
		if line1[0]==">":
			if symbol !="":
				print >>f_out,symbol
			assert line1==line2
			print >>f_out,line1
			symbol = ""
		else:
			seq1=line1
			seq2=line2
			line1=f1.readline().strip()
			line2=f2.readline().strip()
			if line1=="":
				s,h,tot= compare_seq(seq1,seq2,h,tot)
				symbol += s
				print >>f_out,symbol
				break
			if line1[0]==">":
				s,h,tot = compare_seq(seq1,seq2,h,tot)
				symbol += s
				print >>f_out,symbol
				assert line1==line2
				print >>f_out,line1
				symbol = ""
			else:
				seq1+=line1
				seq2+=line2
				s,h,tot = compare_seq(seq1,seq2,h,tot)
				symbol += s
				if len(symbol)==60:
					print >>f_out,symbol
					symbol = ""

if __name__=="__main__":
	file1=sys.argv[1]
	file2=sys.argv[2]
	outfile=sys.argv[3]
	fa_to_psmcfa(file1,file2,outfile)
			
	