import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
from pyscal.trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
import sys as s
import statistics
import math as m
import scipy.spatial.distance as dist
import pickle

INFILE = s.argv[1]
TT = -1 # Timestep to analyze

ALL = 0 # Compute g(r) for all particles if = 1
if ALL != 1: # Compute g(r) for a specific pair if ALL != 1
    TYPE1 = 2
    TYPE2 = 2
NHISTO = 500
Rmax = 5
Rmin = 0.1

Traj = Trajectory(INFILE + '.lammpstrj')

SL = Traj[TT].to_dict()

BOXXa = SL[0]["box"][0][0]
BOXXb = SL[0]["box"][0][1]

BOXYa = SL[0]["box"][1][0]
BOXYb = SL[0]["box"][1][1]

BOXZa = SL[0]["box"][2][0]
BOXZb = SL[0]["box"][2][1]

BOXXL = BOXXb - BOXXa
BOXYL = BOXYb - BOXYa
BOXZL = BOXZb - BOXZa


ID = SL[0]["atoms"]["id"]
XS = SL[0]["atoms"]["xs"]
YS = SL[0]["atoms"]["ys"]
ZS = SL[0]["atoms"]["zs"]
#XR = SL[0]["atoms"]["x"]
#YR = SL[0]["atoms"]["y"]
#ZR = SL[0]["atoms"]["z"]

TY = SL[0]["atoms"]["type"]

print('Data read for timestep ' + str(TT))

N = len(TY)
rho = N/(BOXXL*BOXYL*BOXZL)

X = []
Y = []
Z = []

for i in range(0,N):
	X.append([int(ID[i]),int(TY[i]),float(XS[i])*BOXXL + BOXXa])
	Y.append([int(ID[i]),int(TY[i]),float(YS[i])*BOXYL + BOXYa])
	Z.append([int(ID[i]),int(TY[i]),float(ZS[i])*BOXZL + BOXZa])
    #X.append([int(ID[i]),int(TY[i]),float(XR[i])])
    #Y.append([int(ID[i]),int(TY[i]),float(YR[i])])
    #Z.append([int(ID[i]),int(TY[i]),float(ZR[i])])
    
X.sort()
Y.sort()
Z.sort()

PX = []
PY = []
PZ = []	

for i in range(0,N):
	PX.append([X[i][2]]) 
	PY.append([Y[i][2]]) 
	PZ.append([Z[i][2]]) 
	
	


# Find all distances:

PXA = np.array(PX)
PYA = np.array(PY)
PZA = np.array(PZ)
PXD = dist.pdist(PXA)
PYD = dist.pdist(PYA)			
PZD = dist.pdist(PZA)

PXD[PXD > BOXXL/2] -= BOXXL
PYD[PYD > BOXYL/2] -= BOXYL			
PZD[PZD > BOXZL/2] -= BOXZL



DVEC = np.sqrt(np.power(PXD,2) + np.power(PYD,2) + np.power(PZD,2))

if ALL == 1:
    g = []
    R = []
    [counts,edges] = np.histogram(DVEC,NHISTO) 
    for i in range(0,len(counts)):
        R.append((edges[i]+edges[i+1])/2)
        DR = edges[i+1]-edges[i]
        g.append(2*counts[i]/(4*m.pi*R[i]**2*DR)/rho/N) # 2 beacuse we need to count the i-j and j-i pair
else:
    g = []
    R = []
    DFIL = []
    T1L = []
    T2L = []
    for i in range(0,N):
        if X[i][1] == TYPE1:        
            T1L.append(i)
        if X[i][1] == TYPE2:
            T2L.append(i)

    for i in T1L:
        for j in T2L:
            if i < j:
                IND = N*i + j - ((i + 2)*(i + 1))// 2
                DFIL.append(DVEC[IND])
        
    [counts,edges] = np.histogram(DFIL,NHISTO) 
    for i in range(0,len(counts)):
        R.append((edges[i]+edges[i+1])/2)
        DR = edges[i+1]-edges[i]
        g.append(2*counts[i]/(4*m.pi*R[i]**2*DR)/rho/N/(len(DFIL)/len(DVEC))) # 2 because we need to count the i-j and j-i pair
    


plt.figure()
plt.plot(R,g,'k')
plt.xlim([Rmin,Rmax])
plt.show()	

DATA = [R,g]
with open('gr_' + INFILE + '.pickle','wb') as handle:
    pickle.dump(DATA, handle, protocol=pickle.HIGHEST_PROTOCOL)
				










