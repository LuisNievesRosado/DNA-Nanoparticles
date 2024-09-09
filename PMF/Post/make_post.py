import io

READFILE = 'toread.dat'
OUTFILE = 'post.sh'

MIN = '12.0'
MAX = '50.0'
N = '39'
TOL = '0.0000001'

# Read data file
NUMBERS = []

F1 = open(READFILE,'r')
for line in F1:
	spl = line.split()
	NUMBERS.append(spl[0])
F1.close()	

# Write post file 
F2 = io.open(OUTFILE,'w',newline='\n')
F2.write('#!/bin/bash \n')
F2.write('\n')


F2.write('mkdir sims \n')
for i in range(0,len(NUMBERS)):

	F2.write('cp makehistH.py ' + NUMBERS[i] + '\n')
	F2.write('cp makeprob.py ' + NUMBERS[i] + '\n')
	F2.write('cp PMF2pot.py ' + NUMBERS[i] + '\n')
	F2.write('cp wham ' + NUMBERS[i] + '\n')
	F2.write('cp simrun.sh ' + NUMBERS[i] + '\n')
	F2.write('cp PNP.in ' + NUMBERS[i] + '\n')
	F2.write('cp writerun.py ' + NUMBERS[i] + '\n')
	F2.write('cd ' + NUMBERS[i] + '\n')
	F2.write('python3 makehistH.py \n')
	F2.write('python3 makeprob.py \n')
	F2.write('./wham ' + MIN + ' ' + MAX + ' ' + N + ' ' + TOL + ' ' + '1.0 0.0 MetaAA.wham AA.out \n')
	F2.write('./wham ' + MIN + ' ' + MAX + ' ' + N + ' ' + TOL + ' ' + '1.0 0.0 MetaAB.wham AB.out \n')
	F2.write('./wham ' + MIN + ' ' + MAX + ' ' + N + ' ' + TOL + ' ' + '1.0 0.0 MetaBB.wham BB.out \n')
	F2.write('python3 PMF2pot.py \n')
	F2.write('mkdir ' + NUMBERS[i] + '\n')
	F2.write('mv simrun.sh ' + NUMBERS[i] + '\n')
	F2.write('mv PNP.in ' + NUMBERS[i] + '\n')
	F2.write('mv writerun.py '+ NUMBERS[i] + '\n')
	F2.write('cp AAP.pot '+ NUMBERS[i] + '\n')
	F2.write('cp ABP.pot ' + NUMBERS[i] + '\n')
	F2.write('cp BBP.pot ' + NUMBERS[i] + '\n')
	F2.write('mv ' + NUMBERS[i] + ' ../sims \n')
	F2.write('cd .. \n')
	F2.write('\n')