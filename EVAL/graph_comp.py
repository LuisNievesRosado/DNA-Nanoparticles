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
import scipy.special as sp
import pickle
import networkx as nx

INFILE = s.argv[1] 
TT = -1 # Timestep to analyze
TYPE = 2 # 0 if all the particles, otherwise should be the type of particles to cluster 

CUT = float(s.argv[2])


Traj = Trajectory(INFILE+'.lammpstrj')

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


DXY = np.power(PXD,2) + np.power(PYD,2)
DVEC = np.sqrt(DXY + np.power(PZD,2))

# Make neighbor lists:
NEIGH = []
    
for i in range(0,N):
    NEIGH.append([])
    
    
if TYPE == 0:
    TYL = range(0,N)
else:
    TYL = []
    for i in range(0,N):
        if X[i][1] == TYPE:
            TYL.append(i)
    
    
G = nx.Graph()  
G.add_nodes_from(TYL)


for i in TYL:
    for j in TYL:
        if i < j:
            IND = N*i + j - ((i + 2)*(i + 1))// 2
            if DVEC[IND] < CUT:
                NEIGH[i].append(j)
                NEIGH[j].append(i)
                G.add_edge(i,j)
print('Neighbor lists complete')    

#plt.figure()
#nx.draw(G)
#plt.show()

print('Separating largest cluster')

# Take largest cluster only

CLUS = []
for i in TYL:
    if i == 0:
        CLUS.append([])
        CLUS[-1].append(i)
        Ccl = 1
        while Ccl > 0:
            Ccl = 0
            for j in CLUS[-1]:
                for k in NEIGH[j]:
                    if k in CLUS[-1]:
                        continue
                    else:
                        Ccl = Ccl + 1
                        CLUS[-1].append(k)
    else:
        BCcl = 0
        for j in CLUS:
            if i in j:
                BCcl = BCcl + 1
        if BCcl == 0:
            CLUS.append([])
            CLUS[-1].append(i)
            Ccl = 1
            while Ccl > 0:
                Ccl = 0
                for k in CLUS[-1]:
                    for n in NEIGH[k]:
                        if n in CLUS[-1]:
                            continue
                        else:
                            Ccl = Ccl + 1
                            CLUS[-1].append(n)

L = []
for i in CLUS:
    L.append(len(i))

MAXI = L.index(max(L))


GL = nx.Graph()

for i in CLUS[MAXI]:
    GL.add_node(i)
    GL.add_edges_from(G.edges(i))
    
#plt.figure()
#nx.draw(GL)
#plt.show() 



DATA = [G,GL]
with open(INFILE + '.pickle','wb') as handle:
    pickle.dump(DATA, handle, protocol=pickle.HIGHEST_PROTOCOL)



#print(len(CLUS[MAXI]))
#
#F = open('test.xyz','w')
#F.write(str(len(CLUS[MAXI])) + '\n')
#F.write('\n')
#for i in CLUS[MAXI]:
#    F.write('1 ' + str(X[i][2]) + ' ' + str(Y[i][2]) + ' ' + str(Z[i][2]) + '\n')
#F.close()