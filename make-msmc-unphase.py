import sys
import random
name=sys.argv[1]
freq=float(sys.argv[2])
print freq
for i in range(100):
	file = name+'_%s_sites_2.txt' %(i)
	f=open(file,'r')
	f_out=open(name+'_%s_sites_%.2f.txt' %(i,freq),'w')
	for line in f:
		col=line.rstrip().split('\t')
		phase=col[3]
		s=random.uniform(0,1)
		if phase[0]==phase[1] and phase[2]==phase[3]:
			print >>f_out,'\t'.join(col)
		else:
			if s<=freq:
				alternative=[col[3]]
				if phase[0]!=phase[1]:
					alternative.append(phase[1]+phase[0]+phase[2]+phase[3])
				if phase[2]!=phase[3]:
					alternative.append(phase[0]+phase[1]+phase[3]+phase[2])
				col[3]=','.join(alternative)
				print >>f_out,'\t'.join(col)
			else:
				print >>f_out,'\t'.join(col)