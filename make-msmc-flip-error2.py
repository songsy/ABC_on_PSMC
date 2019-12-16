import sys
import random
name=sys.argv[1]
freq=float(sys.argv[2])
mode=sys.argv[3]
print freq,mode
for i in range(100):
	file = name+'_%s_sites_2.txt' %(i)
	f=open(file,'r')
	f_out=open(name+'_%s_sites_%s_%.2f.txt' %(i,mode,freq),'w')
	for line in f:
		col=line.rstrip().split('\t')
		phase=col[3]
		s=random.uniform(0,1)
		if phase[0]==phase[1] and phase[2]==phase[3]:
			print >>f_out,'\t'.join(col)
		else:
			if phase[0]!=phase[1] and phase[2]==phase[3]:
				if s<=freq and mode=='singleton':
					phase_flip=phase[1]+phase[0]+phase[2]+phase[3]
					col[3]=phase_flip
#					print 'singleton flip', col,phase
			elif phase[2]!=phase[3] and phase[0]==phase[1]:
				if s<=freq and mode=='singleton':
					phase_flip=phase[0]+phase[1]+phase[3]+phase[2]
					col[3]=phase_flip
#					print 'singleton flip', col,phase
			elif phase[0]!=phase[1] and phase[2]!=phase[3]:
				if s<=freq and mode=='doubleton':
					if random.uniform(0,1)<=0.5:
						phase_flip=phase[0]+phase[1]+phase[3]+phase[2]
						col[3]=phase_flip
					else:
						phase_flip=phase[1]+phase[0]+phase[2]+phase[3]
						col[3]=phase_flip
#					print 'doubleton flip', col,phase
			print >>f_out,'\t'.join(col)