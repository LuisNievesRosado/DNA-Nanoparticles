import sys as s
import numpy as np
import os
import math as m

LENGTH_SHORT_POLY_MIN = 1
LENGTH_LONG_POLY_MIN = 2
LENGTH_LONG_POLY_MAX = 15
NUMBER_SHORT_POLY_MIN = 5
NUMBER_LONG_POLY_MIN = 5
NUMBER_SHORT_POLY_MAX = 70
NUMBER_LONG_POLY_MAX = 70

NUMBER_POLY_STEP = 1
NUMBER_POLY_MAX = 100

OUTNAME = 'Designs.out'

F = open(OUTNAME,'w')

count = 0
for k in range(NUMBER_SHORT_POLY_MIN,NUMBER_SHORT_POLY_MAX + NUMBER_POLY_STEP,NUMBER_POLY_STEP):
	for m in range(NUMBER_LONG_POLY_MIN,NUMBER_LONG_POLY_MAX + NUMBER_POLY_STEP,NUMBER_POLY_STEP):
		for i in range(LENGTH_LONG_POLY_MIN,LENGTH_LONG_POLY_MAX + 1):
			for j in range(LENGTH_SHORT_POLY_MIN,i):	
					if m + k > NUMBER_POLY_MAX:
						break
					else:
						count = count + 1
						F.write(str(count) + ' ' + str(k) + ' ' + str(m) + ' ' + str(i) + ' ' + str(j) +'\n')	
F.close()