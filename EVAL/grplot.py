import pickle
import matplotlib.pyplot as plt
import sys as s


INFILE1 = 'gr_' + s.argv[1] + '.pickle'

Rmax = 5
Rmin = 0.1


with open(INFILE1, 'rb') as f1:
    DGR = pickle.load(f1)
	
R = DGR[0]
g = DGR[1]

plt.figure()
plt.plot(R,g,'k')
plt.xlim([Rmin,Rmax])
plt.show()	
	
	