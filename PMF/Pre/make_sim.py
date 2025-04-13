import io

READFILE = 'toread.dat'
DESIGNFILE = 'Designs.out'
OUTFILE = 'params.out'
SIMFILE = 'localsimset.sh'
RUNFILE = 'simrun.sh'


# Read data file
NUMBERS = []

F1 = open(READFILE,'r')
for line in F1:
	spl = line.split()
	NUMBERS.append(spl[0])
F1.close()		
	

# Read design file
PARAMS = []
F2 = open(DESIGNFILE,'r')
for line in F2:
	spl = line.split()
	if spl[0] in NUMBERS:
		PARAMS.append([spl[1],spl[2],spl[3],spl[4]])
F2.close()

# Output parameters
F3 = open(OUTFILE,'w')
for i in range(0,len(NUMBERS)):
	F3.write(NUMBERS[i] + ' ' + PARAMS[i][0] + ' ' + PARAMS[i][1] + ' ' + PARAMS[i][2] + ' ' + PARAMS[i][3] + '\n')
F3.close()	

# Write sim file 
F4 = io.open(SIMFILE,'w',newline='\n')
F4.write('#!/bin/bash \n')
F4.write('\n')

for i in range(0,len(NUMBERS)):
	F4.write('mkdir ' + NUMBERS[i] + '\n')
	F4.write('\n')
	F4.write('python3 write_run.py ' + NUMBERS[i] + '\n')
	F4.write('mv hamilAA.in ' + NUMBERS[i] + '\n')
	F4.write('mv hamilAB.in ' + NUMBERS[i] + '\n')	
	F4.write('mv hamilBB.in ' + NUMBERS[i] + '\n')	
	F4.write('mv runAA.sh ' + NUMBERS[i] + '\n')
	F4.write('mv runAB.sh ' + NUMBERS[i] + '\n')	
	F4.write('mv runBB.sh ' + NUMBERS[i] + '\n')	
	F4.write('cp simset.sh ' + NUMBERS[i] + '\n')	
	F4.write('\n')
	F4.write('mkdir ' + NUMBERS[i] + '/AA \n')
	F4.write('mkdir ' + NUMBERS[i] + '/AB \n')
	F4.write('mkdir ' + NUMBERS[i] + '/BB \n')
	F4.write('\n')
	F4.write('python3 makedatafile.py AA ' + PARAMS[i][0] + ' ' + PARAMS[i][1] + ' ' + PARAMS[i][2] + ' ' + PARAMS[i][3] + '\n')
	F4.write('mv PNP.data ' + NUMBERS[i] + '/AA \n')
	F4.write('cp base.in ' + NUMBERS[i] + '/AA \n')
	F4.write('\n')
	F4.write('python3 makedatafile.py AB ' + PARAMS[i][0] + ' ' + PARAMS[i][1] + ' ' + PARAMS[i][2] + ' ' + PARAMS[i][3] + '\n')
	F4.write('mv PNP.data ' + NUMBERS[i] + '/AB \n')
	F4.write('cp base.in ' + NUMBERS[i] + '/AB \n')	
	F4.write('\n')
	F4.write('python3 makedatafile.py BB ' + PARAMS[i][0] + ' ' + PARAMS[i][1] + ' ' + PARAMS[i][2] + ' ' + PARAMS[i][3] + '\n')
	F4.write('mv PNP.data ' + NUMBERS[i] + '/BB \n')
	F4.write('cp base.in ' + NUMBERS[i] + '/BB \n')
	F4.write('\n')
F4.close()

# Write run file 
F5 = io.open(RUNFILE,'w',newline='\n')
F5.write('#!/bin/bash \n')
F5.write('\n')
for i in range(0,len(NUMBERS)):
	F5.write('cd ' + NUMBERS[i] + '/AA \n')
	F5.write('lmp_stable -in base.in \n')
	F5.write('cd ../.. \n')
	F5.write('cd ' + NUMBERS[i] + '/AB \n')
	F5.write('lmp_stable -in base.in \n')
	F5.write('cd ../.. \n')
	F5.write('cd ' + NUMBERS[i] + '/BB \n')
	F5.write('lmp_stable -in base.in \n')
	F5.write('cd ../.. \n')	
	
	
	

F5.close()


