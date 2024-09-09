import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
from pyscal.trajectory import Trajectory
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import numpy as np
import sys as s
import statistics
import math as m
import scipy.spatial.distance as dist
import pickle

INFILE = s.argv[1] + '.lammpstrj'
TT = -1 # Timestep to analyze


Traj = Trajectory(INFILE)

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



    
X.sort()
Y.sort()
Z.sort()


X1 = []
Y1 = []
Z1 = []

X2 = []
Y2 = []
Z2 = []

for i in range(0,N):
    if X[i][1] == 1:
        X1.append(X[i][2])
        Y1.append(Y[i][2])
        Z1.append(Z[i][2])
    else:    
        X2.append(X[i][2])
        Y2.append(Y[i][2])
        Z2.append(Z[i][2])







# Define cube to section off
L = 10
Lx = L
Ly = L
Lz = L
C = np.array([BOXXa+BOXXL/2,BOXYa+BOXYL/2,BOXZa+BOXZL/2])


a = 0*m.pi/180
b = 0*m.pi/180
g = 0*m.pi/180




R1 = [m.cos(a)*m.cos(b),m.cos(a)*m.sin(b)*m.sin(g)-m.sin(a)*m.cos(g),m.cos(a)*m.sin(b)*m.cos(g)+m.sin(a)*m.sin(g)]
R2 = [m.sin(a)*m.cos(b),m.sin(a)*m.sin(b)*m.sin(g)+m.cos(a)*m.cos(g),m.sin(a)*m.sin(b)*m.cos(g)-m.cos(a)*m.sin(g)]
R3 = [-m.sin(b),m.cos(b)*m.sin(g),m.cos(b)*m.cos(g)]
R = np.array([R1,R2,R3])



P1o = np.array([Lx/2,Ly/2,Lz/2]) 
P2o = np.array([-Lx/2,Ly/2,Lz/2]) 
P3o = np.array([Lx/2,-Ly/2,Lz/2]) 
P4o = np.array([-Lx/2,-Ly/2,Lz/2]) 
P5o = np.array([Lx/2,Ly/2, -Lz/2]) 
P6o = np.array([-Lx/2,Ly/2, -Lz/2]) 
P7o = np.array([Lx/2,-Ly/2, -Lz/2]) 
P8o = np.array([-Lx/2,-Ly/2, -Lz/2]) 

P1 = np.matmul(R,P1o) + C
P2 = np.matmul(R,P2o) + C
P3 = np.matmul(R,P3o) + C
P4 = np.matmul(R,P4o) + C
P5 = np.matmul(R,P5o) + C
P6 = np.matmul(R,P6o) + C
P7 = np.matmul(R,P7o) + C
P8 = np.matmul(R,P8o) + C

P1P2x = [P1[0],P2[0]]
P1P2y = [P1[1],P2[1]]
P1P2z = [P1[2],P2[2]]

P1P3x = [P1[0],P3[0]]
P1P3y = [P1[1],P3[1]]
P1P3z = [P1[2],P3[2]]

P1P5x = [P1[0],P5[0]]
P1P5y = [P1[1],P5[1]]
P1P5z = [P1[2],P5[2]]

P2P4x = [P2[0],P4[0]]
P2P4y = [P2[1],P4[1]]
P2P4z = [P2[2],P4[2]]

P2P6x = [P2[0],P6[0]]
P2P6y = [P2[1],P6[1]]
P2P6z = [P2[2],P6[2]]

P3P4x = [P3[0],P4[0]]
P3P4y = [P3[1],P4[1]]
P3P4z = [P3[2],P4[2]]

P3P7x = [P3[0],P7[0]]
P3P7y = [P3[1],P7[1]]
P3P7z = [P3[2],P7[2]]

P4P8x = [P4[0],P8[0]]
P4P8y = [P4[1],P8[1]]
P4P8z = [P4[2],P8[2]]

P7P8x = [P7[0],P8[0]]
P7P8y = [P7[1],P8[1]]
P7P8z = [P7[2],P8[2]]

P6P8x = [P6[0],P8[0]]
P6P8y = [P6[1],P8[1]]
P6P8z = [P6[2],P8[2]]

P5P6x = [P5[0],P6[0]]
P5P6y = [P5[1],P6[1]]
P5P6z = [P5[2],P6[2]]

P5P7x = [P5[0],P7[0]]
P5P7y = [P5[1],P7[1]]
P5P7z = [P5[2],P7[2]]


fig = plt.figure(figsize=(10,10))
ax= fig.add_subplot(2,1,1,projection='3d', proj_type = 'ortho')

l1, = ax.plot(P1P2x,P1P2y,P1P2z,'k-')
l2, = ax.plot(P1P3x,P1P3y,P1P3z,'k-')
l3, = ax.plot(P1P5x,P1P5y,P1P5z,'k-')
l4, = ax.plot(P2P4x,P2P4y,P2P4z,'k-')
l5, = ax.plot(P2P6x,P2P6y,P2P6z,'k-')
l6, = ax.plot(P3P4x,P3P4y,P3P4z,'k-')
l7, = ax.plot(P3P7x,P3P7y,P3P7z,'k-')
l8, = ax.plot(P4P8x,P4P8y,P4P8z,'k-')
l9, = ax.plot(P7P8x,P7P8y,P7P8z,'k-')
l10, = ax.plot(P6P8x,P6P8y,P6P8z,'k-')
l11, = ax.plot(P5P6x,P5P6y,P5P6z,'k-')
l12, = ax.plot(P5P7x,P5P7y,P5P7z,'k-')

ax.plot(X2,Y2,Z2,'bo',markersize=0.5)
#ax.scatter(X1,Y1,Z1,'r')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


P1o = np.array([Lx/2,Ly/2,Lz/2]) 
P2o = np.array([-Lx/2,Ly/2,Lz/2]) 
P3o = np.array([Lx/2,-Ly/2,Lz/2]) 
P4o = np.array([-Lx/2,-Ly/2,Lz/2]) 
P5o = np.array([Lx/2,Ly/2, -Lz/2]) 
P6o = np.array([-Lx/2,Ly/2, -Lz/2]) 
P7o = np.array([Lx/2,-Ly/2, -Lz/2]) 
P8o = np.array([-Lx/2,-Ly/2, -Lz/2]) 

P1 = np.matmul(R,P1o) + C
P2 = np.matmul(R,P2o) + C
P3 = np.matmul(R,P3o) + C
P4 = np.matmul(R,P4o) + C
P5 = np.matmul(R,P5o) + C
P6 = np.matmul(R,P6o) + C
P7 = np.matmul(R,P7o) + C
P8 = np.matmul(R,P8o) + C


oP1P2x = [P1o[0],P2o[0]]
oP1P2y = [P1o[1],P2o[1]]
oP1P2z = [P1o[2],P2o[2]]

oP1P3x = [P1o[0],P3o[0]]
oP1P3y = [P1o[1],P3o[1]]
oP1P3z = [P1o[2],P3o[2]]

oP1P5x = [P1o[0],P5o[0]]
oP1P5y = [P1o[1],P5o[1]]
oP1P5z = [P1o[2],P5o[2]]

oP2P4x = [P2o[0],P4o[0]]
oP2P4y = [P2o[1],P4o[1]]
oP2P4z = [P2o[2],P4o[2]]

oP2P6x = [P2o[0],P6o[0]]
oP2P6y = [P2o[1],P6o[1]]
oP2P6z = [P2o[2],P6o[2]]

oP3P4x = [P3o[0],P4o[0]]
oP3P4y = [P3o[1],P4o[1]]
oP3P4z = [P3o[2],P4o[2]]

oP3P7x = [P3o[0],P7o[0]]
oP3P7y = [P3o[1],P7o[1]]
oP3P7z = [P3o[2],P7o[2]]

oP4P8x = [P4o[0],P8o[0]]
oP4P8y = [P4o[1],P8o[1]]
oP4P8z = [P4o[2],P8o[2]]

oP7P8x = [P7o[0],P8o[0]]
oP7P8y = [P7o[1],P8o[1]]
oP7P8z = [P7o[2],P8o[2]]

oP6P8x = [P6o[0],P8o[0]]
oP6P8y = [P6o[1],P8o[1]]
oP6P8z = [P6o[2],P8o[2]]

oP5P6x = [P5o[0],P6o[0]]
oP5P6y = [P5o[1],P6o[1]]
oP5P6z = [P5o[2],P6o[2]]

oP5P7x = [P5o[0],P7o[0]]
oP5P7y = [P5o[1],P7o[1]]
oP5P7z = [P5o[2],P7o[2]]


Rinv = np.linalg.inv(R)

# Go point by point and check if it is inside cube


Xl = (P5 - P6)/Lx
Yl = (P8 - P6)/Ly
Zl = (P2 - P6)/Lz

X2in = []
Y2in = []
Z2in = []
for i in range(0,len(X2)):
    V = np.array([X2[i]-C[0],Y2[i]-C[1],Z2[i]-C[2]])
    px = abs(np.dot(V,Xl))
    py = abs(np.dot(V,Yl))
    pz = abs(np.dot(V,Zl))
    if 2*px <= Lx and 2*py <= Ly and 2*pz <= Lz:
        X2in.append(X2[i])
        Y2in.append(Y2[i])
        Z2in.append(Z2[i])
                     
X1in = []
Y1in = []
Z1in = []                     
for i in range(0,len(X1)):
    V = np.array([X1[i]-C[0],Y1[i]-C[1],Z1[i]-C[2]])
    px = abs(np.dot(V,Xl))
    py = abs(np.dot(V,Yl))
    pz = abs(np.dot(V,Zl))
    if 2*px <= Lx and 2*py <= Ly and 2*pz <= Lz:
        X1in.append(X1[i])
        Y1in.append(Y1[i])
        Z1in.append(Z1[i])        
        
# Rotate points inside cube back to base

oX2in = []
oY2in = []
oZ2in = []

for i in range(0,len(X2in)):
    Pt = [X2in[i]-C[0],Y2in[i]-C[1],Z2in[i]-C[2]]
    Pto = np.matmul(Rinv,Pt)
    oX2in.append(Pto[0])
    oY2in.append(Pto[1])
    oZ2in.append(Pto[2])

oX1in = []
oY1in = []
oZ1in = []

for i in range(0,len(X1in)):
    Pt = [X1in[i]-C[0],Y1in[i]-C[1],Z1in[i]-C[2]]
    Pto = np.matmul(Rinv,Pt)
    oX1in.append(Pto[0])
    oY1in.append(Pto[1])
    oZ1in.append(Pto[2])



# Plot points with periodic BC
bX1  = []
bX2  = []
bX3  = []
bX4  = []
bX5  = []
bX6  = []
bX7  = []
bX8  = []
bX9  = []
bX10 = []
bX11 = []
bX12 = []
bX13 = []
bX14 = []
bX15 = []
bX16 = []
bX17 = []
bX18 = []
bX19 = []
bX20 = []
bX21 = []
bX22 = []
bX23 = []
bX24 = []
bX25 = []
bX26 = []
bX27 = []

bY1  = []
bY2  = []
bY3  = []
bY4  = []
bY5  = []
bY6  = []
bY7  = []
bY8  = []
bY9  = []
bY10 = []
bY11 = []
bY12 = []
bY13 = []
bY14 = []
bY15 = []
bY16 = []
bY17 = []
bY18 = []
bY19 = []
bY20 = []
bY21 = []
bY22 = []
bY23 = []
bY24 = []
bY25 = []
bY26 = []
bY27 = []

bZ1  = []
bZ2  = []
bZ3  = []
bZ4  = []
bZ5  = []
bZ6  = []
bZ7  = []
bZ8  = []
bZ9  = []
bZ10 = []
bZ11 = []
bZ12 = []
bZ13 = []
bZ14 = []
bZ15 = []
bZ16 = []
bZ17 = []
bZ18 = []
bZ19 = []
bZ20 = []
bZ21 = []
bZ22 = []
bZ23 = []
bZ24 = []
bZ25 = []
bZ26 = []
bZ27 = []

for i in range(0,len(X2in)):
    bX1.append(oX2in[i]+Lx)   
    bX2.append(oX2in[i]+Lx)
    bX3.append(oX2in[i]+Lx)
    bX4.append(oX2in[i])
    bX5.append(oX2in[i])
    bX6.append(oX2in[i])
    bX7.append(oX2in[i]-Lx)
    bX8.append(oX2in[i]-Lx)
    bX9.append(oX2in[i]-Lx)
    bX10.append(oX2in[i]+Lx)
    bX11.append(oX2in[i]+Lx)
    bX12.append(oX2in[i]+Lx)
    bX13.append(oX2in[i])
    bX14.append(oX2in[i])
    bX15.append(oX2in[i])
    bX16.append(oX2in[i]-Lx)
    bX17.append(oX2in[i]-Lx)
    bX18.append(oX2in[i]-Lx)
    bX19.append(oX2in[i]+Lx)
    bX20.append(oX2in[i]+Lx)
    bX21.append(oX2in[i]+Lx)
    bX22.append(oX2in[i])
    bX23.append(oX2in[i])
    bX24.append(oX2in[i])
    bX25.append(oX2in[i]-Lx)
    bX26.append(oX2in[i]-Lx)
    bX27.append(oX2in[i]-Lx)
    bY1.append(oY2in[i]+Ly)
    bY2.append(oY2in[i])
    bY3.append(oY2in[i]-Ly)
    bY4.append(oY2in[i]+Ly)
    bY5.append(oY2in[i])
    bY6.append(oY2in[i]-Ly)
    bY7.append(oY2in[i]+Ly)
    bY8.append(oY2in[i])
    bY9.append(oY2in[i]-Ly)
    bY10.append(oY2in[i]+Ly)
    bY11.append(oY2in[i])
    bY12.append(oY2in[i]-Ly)
    bY13.append(oY2in[i]+Ly)
    bY14.append(oY2in[i])
    bY15.append(oY2in[i]-Ly)
    bY16.append(oY2in[i]+Ly)
    bY17.append(oY2in[i])
    bY18.append(oY2in[i]-Ly)
    bY19.append(oY2in[i]+Ly)
    bY20.append(oY2in[i])
    bY21.append(oY2in[i]-Ly)
    bY22.append(oY2in[i]+Ly)
    bY23.append(oY2in[i])
    bY24.append(oY2in[i]-Ly)
    bY25.append(oY2in[i]+Ly)
    bY26.append(oY2in[i])
    bY27.append(oY2in[i]-Ly)
    bZ1.append(oZ2in[i]-Lz)
    bZ2.append(oZ2in[i]-Lz)
    bZ3.append(oZ2in[i]-Lz)
    bZ4.append(oZ2in[i]-Lz)
    bZ5.append(oZ2in[i]-Lz)
    bZ6.append(oZ2in[i]-Lz)
    bZ7.append(oZ2in[i]-Lz)
    bZ8.append(oZ2in[i]-Lz)
    bZ9.append(oZ2in[i]-Lz)
    bZ10.append(oZ2in[i])
    bZ11.append(oZ2in[i])
    bZ12.append(oZ2in[i])
    bZ13.append(oZ2in[i])
    bZ14.append(oZ2in[i])
    bZ15.append(oZ2in[i])
    bZ16.append(oZ2in[i])
    bZ17.append(oZ2in[i])
    bZ18.append(oZ2in[i])
    bZ19.append(oZ2in[i]+Lz)
    bZ20.append(oZ2in[i]+Lz)
    bZ21.append(oZ2in[i]+Lz)
    bZ22.append(oZ2in[i]+Lz)
    bZ23.append(oZ2in[i]+Lz)
    bZ24.append(oZ2in[i]+Lz)
    bZ25.append(oZ2in[i]+Lz)
    bZ26.append(oZ2in[i]+Lz)
    bZ27.append(oZ2in[i]+Lz)
    

ax2= fig.add_subplot(2,1,2,projection='3d', proj_type = 'ortho')
b1, = ax2.plot(oP1P2x,oP1P2y,oP1P2z,'k-')
b2, = ax2.plot(oP1P3x,oP1P3y,oP1P3z,'k-')
b3, = ax2.plot(oP1P5x,oP1P5y,oP1P5z,'k-')
b4, = ax2.plot(oP2P4x,oP2P4y,oP2P4z,'k-')
b5, = ax2.plot(oP2P6x,oP2P6y,oP2P6z,'k-')
b6, = ax2.plot(oP3P4x,oP3P4y,oP3P4z,'k-')
b7, = ax2.plot(oP3P7x,oP3P7y,oP3P7z,'k-')
b8, = ax2.plot(oP4P8x,oP4P8y,oP4P8z,'k-')
b9, = ax2.plot(oP7P8x,oP7P8y,oP7P8z,'k-')
b10, = ax2.plot(oP6P8x,oP6P8y,oP6P8z,'k-')
b11, = ax2.plot(oP5P6x,oP5P6y,oP5P6z,'k-')
b12, = ax2.plot(oP5P7x,oP5P7y,oP5P7z,'k-')
p1, =  ax2.plot(bX1,bY1,bZ1,'bo',markersize=0.5)
p2, =  ax2.plot(bX2,bY2,bZ2,'bo',markersize=0.5)
p3, =  ax2.plot(bX3,bY3,bZ3,'bo',markersize=0.5)
p4, =  ax2.plot(bX4,bY4,bZ4,'bo',markersize=0.5)
p5, =  ax2.plot(bX5,bY5,bZ5,'bo',markersize=0.5)
p6, =  ax2.plot(bX6,bY6,bZ6,'bo',markersize=0.5)
p7, =  ax2.plot(bX7,bY7,bZ7,'bo',markersize=0.5)
p8, =  ax2.plot(bX8,bY8,bZ8,'bo',markersize=0.5)
p9, =  ax2.plot(bX9,bY9,bZ9,'bo',markersize=0.5)
p10, = ax2.plot(bX10,bY10,bZ10,'bo',markersize=0.5)
p11, = ax2.plot(bX11,bY11,bZ11,'bo',markersize=0.5)
p12, = ax2.plot(bX12,bY12,bZ12,'bo',markersize=0.5)
p13, = ax2.plot(bX13,bY13,bZ13,'bo',markersize=0.5)
p14, = ax2.plot(bX14,bY14,bZ14,'bo',markersize=0.5)
p15, = ax2.plot(bX15,bY15,bZ15,'bo',markersize=0.5)
p16, = ax2.plot(bX16,bY16,bZ16,'bo',markersize=0.5)
p17, = ax2.plot(bX17,bY17,bZ17,'bo',markersize=0.5)
p18, = ax2.plot(bX18,bY18,bZ18,'bo',markersize=0.5)
p19, = ax2.plot(bX19,bY19,bZ19,'bo',markersize=0.5)
p20, = ax2.plot(bX20,bY20,bZ20,'bo',markersize=0.5)
p21, = ax2.plot(bX21,bY21,bZ21,'bo',markersize=0.5)
p22, = ax2.plot(bX22,bY22,bZ22,'bo',markersize=0.5)
p23, = ax2.plot(bX23,bY23,bZ23,'bo',markersize=0.5)
p24, = ax2.plot(bX24,bY24,bZ24,'bo',markersize=0.5)
p25, = ax2.plot(bX25,bY25,bZ25,'bo',markersize=0.5)
p26, = ax2.plot(bX26,bY26,bZ26,'bo',markersize=0.5)
p27, = ax2.plot(bX27,bY27,bZ27,'bo',markersize=0.5) 
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')






















fig.subplots_adjust(left=-0.25)
ax_L = fig.add_axes([0.65, 0.1, 0.03, 0.65])
Lsli = Slider(ax_L,label='L',valmin=0,valmax=2*L,valinit = L,valstep=0.1,orientation='vertical')

ax_Cx = fig.add_axes([0.7, 0.1, 0.03,0.65])
Cxsli = Slider(ax_Cx,label='Cx',valmin=BOXXa,valmax=BOXXb,valinit = C[0],valstep=0.1,orientation='vertical')

ax_Cy = fig.add_axes([0.75, 0.1, 0.03, 0.65])
Cysli = Slider(ax_Cy,label='Cy',valmin=BOXYa,valmax=BOXYb,valinit = C[1],valstep=0.1,orientation='vertical')

ax_Cz = fig.add_axes([0.8, 0.1, 0.03, 0.65])
Czsli = Slider(ax_Cz,label='Cz',valmin=BOXZa,valmax=BOXZb,valinit = C[2],valstep=0.1,orientation='vertical')

ax_a = fig.add_axes([0.85, 0.1, 0.03, 0.65])
asli = Slider(ax_a,label='alpha',valmin=0,valmax=90,valinit = 0,valstep=0.1,orientation='vertical')

ax_b = fig.add_axes([0.9, 0.1, 0.03, 0.65])
bsli = Slider(ax_b,label='beta',valmin=0,valmax=90,valinit = 0,valstep=0.1,orientation='vertical')

ax_g = fig.add_axes([0.95, 0.1, 0.03, 0.65])
gsli = Slider(ax_g,label='gamma',valmin=0,valmax=90,valinit = 0,valstep=0.1,orientation='vertical')





def update(val):
    global a
    global b
    global g 
    global C
    global Lx
    global Ly
    global Lz
    
    global P1P2x
    global P1P2y
    global P1P2z
    global P1P3x
    global P1P3y
    global P1P3z
    global P1P5x
    global P1P5y
    global P1P5z    
    global P2P4x
    global P2P4y
    global P2P4z    
    global P2P6x
    global P2P6y
    global P2P6z    
    global P3P4x
    global P3P4y
    global P3P4z    
    global P3P7x
    global P3P7y
    global P3P7z    
    global P4P8x
    global P4P8y
    global P4P8z    
    global P7P8x
    global P7P8y
    global P7P8z    
    global P6P8x
    global P6P8y
    global P6P8z    
    global P5P6x
    global P5P6y
    global P5P6z    
    global P5P7x
    global P5P7y
    global P5P7z    
    
    global R
    
    Lx = Lsli.val
    Ly = Lsli.val
    Lz = Lsli.val
    
    C = np.array([Cxsli.val,Cysli.val,Czsli.val])
    
    a = asli.val*m.pi/180
    b = bsli.val*m.pi/180
    g = gsli.val*m.pi/180

    R1 = [m.cos(a)*m.cos(b),m.cos(a)*m.sin(b)*m.sin(g)-m.sin(a)*m.cos(g),m.cos(a)*m.sin(b)*m.cos(g)+m.sin(a)*m.sin(g)]
    R2 = [m.sin(a)*m.cos(b),m.sin(a)*m.sin(b)*m.sin(g)+m.cos(a)*m.cos(g),m.sin(a)*m.sin(b)*m.cos(g)-m.cos(a)*m.sin(g)]
    R3 = [-m.sin(b),m.cos(b)*m.sin(g),m.cos(b)*m.cos(g)]
    R = np.array([R1,R2,R3])
    
    Rinv = np.linalg.inv(R)

    P1o = np.array([Lx/2,Ly/2,Lz/2]) 
    P2o = np.array([-Lx/2,Ly/2,Lz/2]) 
    P3o = np.array([Lx/2,-Ly/2,Lz/2]) 
    P4o = np.array([-Lx/2,-Ly/2,Lz/2]) 
    P5o = np.array([Lx/2,Ly/2, -Lz/2]) 
    P6o = np.array([-Lx/2,Ly/2, -Lz/2]) 
    P7o = np.array([Lx/2,-Ly/2, -Lz/2]) 
    P8o = np.array([-Lx/2,-Ly/2, -Lz/2]) 
    
    P1 = np.matmul(R,P1o) + C
    P2 = np.matmul(R,P2o) + C
    P3 = np.matmul(R,P3o) + C
    P4 = np.matmul(R,P4o) + C
    P5 = np.matmul(R,P5o) + C
    P6 = np.matmul(R,P6o) + C
    P7 = np.matmul(R,P7o) + C
    P8 = np.matmul(R,P8o) + C
    
    P1P2x = [P1[0],P2[0]]
    P1P2y = [P1[1],P2[1]]
    P1P2z = [P1[2],P2[2]]
    
    P1P3x = [P1[0],P3[0]]
    P1P3y = [P1[1],P3[1]]
    P1P3z = [P1[2],P3[2]]
    
    P1P5x = [P1[0],P5[0]]
    P1P5y = [P1[1],P5[1]]
    P1P5z = [P1[2],P5[2]]
    
    P2P4x = [P2[0],P4[0]]
    P2P4y = [P2[1],P4[1]]
    P2P4z = [P2[2],P4[2]]
    
    P2P6x = [P2[0],P6[0]]
    P2P6y = [P2[1],P6[1]]
    P2P6z = [P2[2],P6[2]]
    
    P3P4x = [P3[0],P4[0]]
    P3P4y = [P3[1],P4[1]]
    P3P4z = [P3[2],P4[2]]
    
    P3P7x = [P3[0],P7[0]]
    P3P7y = [P3[1],P7[1]]
    P3P7z = [P3[2],P7[2]]
    
    P4P8x = [P4[0],P8[0]]
    P4P8y = [P4[1],P8[1]]
    P4P8z = [P4[2],P8[2]]
    
    P7P8x = [P7[0],P8[0]]
    P7P8y = [P7[1],P8[1]]
    P7P8z = [P7[2],P8[2]]
    
    P6P8x = [P6[0],P8[0]]
    P6P8y = [P6[1],P8[1]]
    P6P8z = [P6[2],P8[2]]
    
    P5P6x = [P5[0],P6[0]]
    P5P6y = [P5[1],P6[1]]
    P5P6z = [P5[2],P6[2]]
    
    P5P7x = [P5[0],P7[0]]
    P5P7y = [P5[1],P7[1]]
    P5P7z = [P5[2],P7[2]]
    
    
    oP1P2x = [P1o[0],P2o[0]]
    oP1P2y = [P1o[1],P2o[1]]
    oP1P2z = [P1o[2],P2o[2]]
    
    oP1P3x = [P1o[0],P3o[0]]
    oP1P3y = [P1o[1],P3o[1]]
    oP1P3z = [P1o[2],P3o[2]]
    
    oP1P5x = [P1o[0],P5o[0]]
    oP1P5y = [P1o[1],P5o[1]]
    oP1P5z = [P1o[2],P5o[2]]
    
    oP2P4x = [P2o[0],P4o[0]]
    oP2P4y = [P2o[1],P4o[1]]
    oP2P4z = [P2o[2],P4o[2]]
    
    oP2P6x = [P2o[0],P6o[0]]
    oP2P6y = [P2o[1],P6o[1]]
    oP2P6z = [P2o[2],P6o[2]]
    
    oP3P4x = [P3o[0],P4o[0]]
    oP3P4y = [P3o[1],P4o[1]]
    oP3P4z = [P3o[2],P4o[2]]
    
    oP3P7x = [P3o[0],P7o[0]]
    oP3P7y = [P3o[1],P7o[1]]
    oP3P7z = [P3o[2],P7o[2]]
    
    oP4P8x = [P4o[0],P8o[0]]
    oP4P8y = [P4o[1],P8o[1]]
    oP4P8z = [P4o[2],P8o[2]]
    
    oP7P8x = [P7o[0],P8o[0]]
    oP7P8y = [P7o[1],P8o[1]]
    oP7P8z = [P7o[2],P8o[2]]
    
    oP6P8x = [P6o[0],P8o[0]]
    oP6P8y = [P6o[1],P8o[1]]
    oP6P8z = [P6o[2],P8o[2]]
    
    oP5P6x = [P5o[0],P6o[0]]
    oP5P6y = [P5o[1],P6o[1]]
    oP5P6z = [P5o[2],P6o[2]]
    
    oP5P7x = [P5o[0],P7o[0]]
    oP5P7y = [P5o[1],P7o[1]]
    oP5P7z = [P5o[2],P7o[2]]

    l1.set_data_3d(P1P2x,P1P2y,P1P2z)
    l2.set_data_3d(P1P3x,P1P3y,P1P3z)    
    l3.set_data_3d(P1P5x,P1P5y,P1P5z)    
    l4.set_data_3d(P2P4x,P2P4y,P2P4z)    
    l5.set_data_3d(P2P6x,P2P6y,P2P6z)    
    l6.set_data_3d(P3P4x,P3P4y,P3P4z)    
    l7.set_data_3d(P3P7x,P3P7y,P3P7z)    
    l8.set_data_3d(P4P8x,P4P8y,P4P8z)    
    l9.set_data_3d(P7P8x,P7P8y,P7P8z)    
    l10.set_data_3d(P6P8x,P6P8y,P6P8z) 
    l11.set_data_3d(P5P6x,P5P6y,P5P6z)
    l12.set_data_3d(P5P7x,P5P7y,P5P7z)


    Xl = (P5 - P6)/Lx
    Yl = (P8 - P6)/Ly
    Zl = (P2 - P6)/Lz
    
    X2in = []
    Y2in = []
    Z2in = []
    for i in range(0,len(X2)):
        V = np.array([X2[i]-C[0],Y2[i]-C[1],Z2[i]-C[2]])
        px = abs(np.dot(V,Xl))
        py = abs(np.dot(V,Yl))
        pz = abs(np.dot(V,Zl))
        if 2*px <= Lx and 2*py <= Ly and 2*pz <= Lz:
            X2in.append(X2[i])
            Y2in.append(Y2[i])
            Z2in.append(Z2[i])
            
    # Rotate points inside cube back to base
    
    oX2in = []
    oY2in = []
    oZ2in = []
    
    for i in range(0,len(X2in)):
        Pt = [X2in[i]-C[0],Y2in[i]-C[1],Z2in[i]-C[2]]
        Pto = np.matmul(Rinv,Pt)
        oX2in.append(Pto[0])
        oY2in.append(Pto[1])
        oZ2in.append(Pto[2])
    
    
    # Plot points with periodic BC
    bX1  = []
    bX2  = []
    bX3  = []
    bX4  = []
    bX5  = []
    bX6  = []
    bX7  = []
    bX8  = []
    bX9  = []
    bX10 = []
    bX11 = []
    bX12 = []
    bX13 = []
    bX14 = []
    bX15 = []
    bX16 = []
    bX17 = []
    bX18 = []
    bX19 = []
    bX20 = []
    bX21 = []
    bX22 = []
    bX23 = []
    bX24 = []
    bX25 = []
    bX26 = []
    bX27 = []
    
    bY1  = []
    bY2  = []
    bY3  = []
    bY4  = []
    bY5  = []
    bY6  = []
    bY7  = []
    bY8  = []
    bY9  = []
    bY10 = []
    bY11 = []
    bY12 = []
    bY13 = []
    bY14 = []
    bY15 = []
    bY16 = []
    bY17 = []
    bY18 = []
    bY19 = []
    bY20 = []
    bY21 = []
    bY22 = []
    bY23 = []
    bY24 = []
    bY25 = []
    bY26 = []
    bY27 = []
    
    bZ1  = []
    bZ2  = []
    bZ3  = []
    bZ4  = []
    bZ5  = []
    bZ6  = []
    bZ7  = []
    bZ8  = []
    bZ9  = []
    bZ10 = []
    bZ11 = []
    bZ12 = []
    bZ13 = []
    bZ14 = []
    bZ15 = []
    bZ16 = []
    bZ17 = []
    bZ18 = []
    bZ19 = []
    bZ20 = []
    bZ21 = []
    bZ22 = []
    bZ23 = []
    bZ24 = []
    bZ25 = []
    bZ26 = []
    bZ27 = []
    
    for i in range(0,len(X2in)):
        bX1.append(oX2in[i]+Lx)   
        bX2.append(oX2in[i]+Lx)
        bX3.append(oX2in[i]+Lx)
        bX4.append(oX2in[i])
        bX5.append(oX2in[i])
        bX6.append(oX2in[i])
        bX7.append(oX2in[i]-Lx)
        bX8.append(oX2in[i]-Lx)
        bX9.append(oX2in[i]-Lx)
        bX10.append(oX2in[i]+Lx)
        bX11.append(oX2in[i]+Lx)
        bX12.append(oX2in[i]+Lx)
        bX13.append(oX2in[i])
        bX14.append(oX2in[i])
        bX15.append(oX2in[i])
        bX16.append(oX2in[i]-Lx)
        bX17.append(oX2in[i]-Lx)
        bX18.append(oX2in[i]-Lx)
        bX19.append(oX2in[i]+Lx)
        bX20.append(oX2in[i]+Lx)
        bX21.append(oX2in[i]+Lx)
        bX22.append(oX2in[i])
        bX23.append(oX2in[i])
        bX24.append(oX2in[i])
        bX25.append(oX2in[i]-Lx)
        bX26.append(oX2in[i]-Lx)
        bX27.append(oX2in[i]-Lx)
        bY1.append(oY2in[i]+Ly)
        bY2.append(oY2in[i])
        bY3.append(oY2in[i]-Ly)
        bY4.append(oY2in[i]+Ly)
        bY5.append(oY2in[i])
        bY6.append(oY2in[i]-Ly)
        bY7.append(oY2in[i]+Ly)
        bY8.append(oY2in[i])
        bY9.append(oY2in[i]-Ly)
        bY10.append(oY2in[i]+Ly)
        bY11.append(oY2in[i])
        bY12.append(oY2in[i]-Ly)
        bY13.append(oY2in[i]+Ly)
        bY14.append(oY2in[i])
        bY15.append(oY2in[i]-Ly)
        bY16.append(oY2in[i]+Ly)
        bY17.append(oY2in[i])
        bY18.append(oY2in[i]-Ly)
        bY19.append(oY2in[i]+Ly)
        bY20.append(oY2in[i])
        bY21.append(oY2in[i]-Ly)
        bY22.append(oY2in[i]+Ly)
        bY23.append(oY2in[i])
        bY24.append(oY2in[i]-Ly)
        bY25.append(oY2in[i]+Ly)
        bY26.append(oY2in[i])
        bY27.append(oY2in[i]-Ly)
        bZ1.append(oZ2in[i]-Lz)
        bZ2.append(oZ2in[i]-Lz)
        bZ3.append(oZ2in[i]-Lz)
        bZ4.append(oZ2in[i]-Lz)
        bZ5.append(oZ2in[i]-Lz)
        bZ6.append(oZ2in[i]-Lz)
        bZ7.append(oZ2in[i]-Lz)
        bZ8.append(oZ2in[i]-Lz)
        bZ9.append(oZ2in[i]-Lz)
        bZ10.append(oZ2in[i])
        bZ11.append(oZ2in[i])
        bZ12.append(oZ2in[i])
        bZ13.append(oZ2in[i])
        bZ14.append(oZ2in[i])
        bZ15.append(oZ2in[i])
        bZ16.append(oZ2in[i])
        bZ17.append(oZ2in[i])
        bZ18.append(oZ2in[i])
        bZ19.append(oZ2in[i]+Lz)
        bZ20.append(oZ2in[i]+Lz)
        bZ21.append(oZ2in[i]+Lz)
        bZ22.append(oZ2in[i]+Lz)
        bZ23.append(oZ2in[i]+Lz)
        bZ24.append(oZ2in[i]+Lz)
        bZ25.append(oZ2in[i]+Lz)
        bZ26.append(oZ2in[i]+Lz)
        bZ27.append(oZ2in[i]+Lz)

    b1.set_data_3d(oP1P2x,oP1P2y,oP1P2z)
    b2.set_data_3d(oP1P3x,oP1P3y,oP1P3z)
    b3.set_data_3d(oP1P5x,oP1P5y,oP1P5z)
    b4.set_data_3d(oP2P4x,oP2P4y,oP2P4z)
    b5.set_data_3d(oP2P6x,oP2P6y,oP2P6z)
    b6.set_data_3d(oP3P4x,oP3P4y,oP3P4z)
    b7.set_data_3d(oP3P7x,oP3P7y,oP3P7z)
    b8.set_data_3d(oP4P8x,oP4P8y,oP4P8z)
    b9.set_data_3d(oP7P8x,oP7P8y,oP7P8z)
    b10.set_data_3d(oP6P8x,oP6P8y,oP6P8z)
    b11.set_data_3d(oP5P6x,oP5P6y,oP5P6z)
    b12.set_data_3d(oP5P7x,oP5P7y,oP5P7z)

    p1.set_data_3d(bX1,bY1,bZ1)
    p2.set_data_3d(bX2,bY2,bZ2)
    p3.set_data_3d(bX3,bY3,bZ3)
    p4.set_data_3d(bX4,bY4,bZ4)
    p5.set_data_3d(bX5,bY5,bZ5)
    p6.set_data_3d(bX6,bY6,bZ6)
    p7.set_data_3d(bX7,bY7,bZ7)
    p8.set_data_3d(bX8,bY8,bZ8)
    p9.set_data_3d(bX9,bY9,bZ9)
    p10.set_data_3d(bX10,bY10,bZ10)
    p11.set_data_3d(bX11,bY11,bZ11)
    p12.set_data_3d(bX12,bY12,bZ12)
    p13.set_data_3d(bX13,bY13,bZ13)
    p14.set_data_3d(bX14,bY14,bZ14)
    p15.set_data_3d(bX15,bY15,bZ15)
    p16.set_data_3d(bX16,bY16,bZ16)
    p17.set_data_3d(bX17,bY17,bZ17)
    p18.set_data_3d(bX18,bY18,bZ18)
    p19.set_data_3d(bX19,bY19,bZ19)
    p20.set_data_3d(bX20,bY20,bZ20)
    p21.set_data_3d(bX21,bY21,bZ21)
    p22.set_data_3d(bX22,bY22,bZ22)
    p23.set_data_3d(bX23,bY23,bZ23)
    p24.set_data_3d(bX24,bY24,bZ24)
    p25.set_data_3d(bX25,bY25,bZ25)
    p26.set_data_3d(bX26,bY26,bZ26)
    p27.set_data_3d(bX27,bY27,bZ27) 































Lsli.on_changed(update)
Cxsli.on_changed(update)
Cysli.on_changed(update)
Czsli.on_changed(update)
asli.on_changed(update)
bsli.on_changed(update)
gsli.on_changed(update)


plt.show()


P1o = np.array([Lx/2,Ly/2,Lz/2]) 
P2o = np.array([-Lx/2,Ly/2,Lz/2]) 
P3o = np.array([Lx/2,-Ly/2,Lz/2]) 
P4o = np.array([-Lx/2,-Ly/2,Lz/2]) 
P5o = np.array([Lx/2,Ly/2, -Lz/2]) 
P6o = np.array([-Lx/2,Ly/2, -Lz/2]) 
P7o = np.array([Lx/2,-Ly/2, -Lz/2]) 
P8o = np.array([-Lx/2,-Ly/2, -Lz/2]) 

P1 = np.matmul(R,P1o) + C
P2 = np.matmul(R,P2o) + C
P3 = np.matmul(R,P3o) + C
P4 = np.matmul(R,P4o) + C
P5 = np.matmul(R,P5o) + C
P6 = np.matmul(R,P6o) + C
P7 = np.matmul(R,P7o) + C
P8 = np.matmul(R,P8o) + C


oP1P2x = [P1o[0],P2o[0]]
oP1P2y = [P1o[1],P2o[1]]
oP1P2z = [P1o[2],P2o[2]]

oP1P3x = [P1o[0],P3o[0]]
oP1P3y = [P1o[1],P3o[1]]
oP1P3z = [P1o[2],P3o[2]]

oP1P5x = [P1o[0],P5o[0]]
oP1P5y = [P1o[1],P5o[1]]
oP1P5z = [P1o[2],P5o[2]]

oP2P4x = [P2o[0],P4o[0]]
oP2P4y = [P2o[1],P4o[1]]
oP2P4z = [P2o[2],P4o[2]]

oP2P6x = [P2o[0],P6o[0]]
oP2P6y = [P2o[1],P6o[1]]
oP2P6z = [P2o[2],P6o[2]]

oP3P4x = [P3o[0],P4o[0]]
oP3P4y = [P3o[1],P4o[1]]
oP3P4z = [P3o[2],P4o[2]]

oP3P7x = [P3o[0],P7o[0]]
oP3P7y = [P3o[1],P7o[1]]
oP3P7z = [P3o[2],P7o[2]]

oP4P8x = [P4o[0],P8o[0]]
oP4P8y = [P4o[1],P8o[1]]
oP4P8z = [P4o[2],P8o[2]]

oP7P8x = [P7o[0],P8o[0]]
oP7P8y = [P7o[1],P8o[1]]
oP7P8z = [P7o[2],P8o[2]]

oP6P8x = [P6o[0],P8o[0]]
oP6P8y = [P6o[1],P8o[1]]
oP6P8z = [P6o[2],P8o[2]]

oP5P6x = [P5o[0],P6o[0]]
oP5P6y = [P5o[1],P6o[1]]
oP5P6z = [P5o[2],P6o[2]]

oP5P7x = [P5o[0],P7o[0]]
oP5P7y = [P5o[1],P7o[1]]
oP5P7z = [P5o[2],P7o[2]]


Rinv = np.linalg.inv(R)

# Go point by point and check if it is inside cube


Xl = (P5 - P6)/Lx
Yl = (P8 - P6)/Ly
Zl = (P2 - P6)/Lz

X2in = []
Y2in = []
Z2in = []
for i in range(0,len(X2)):
    V = np.array([X2[i]-C[0],Y2[i]-C[1],Z2[i]-C[2]])
    px = abs(np.dot(V,Xl))
    py = abs(np.dot(V,Yl))
    pz = abs(np.dot(V,Zl))
    if 2*px <= Lx and 2*py <= Ly and 2*pz <= Lz:
        X2in.append(X2[i])
        Y2in.append(Y2[i])
        Z2in.append(Z2[i])
                     
X1in = []
Y1in = []
Z1in = []                     
for i in range(0,len(X1)):
    V = np.array([X1[i]-C[0],Y1[i]-C[1],Z1[i]-C[2]])
    px = abs(np.dot(V,Xl))
    py = abs(np.dot(V,Yl))
    pz = abs(np.dot(V,Zl))
    if 2*px <= Lx and 2*py <= Ly and 2*pz <= Lz:
        X1in.append(X1[i])
        Y1in.append(Y1[i])
        Z1in.append(Z1[i])        
        
# Rotate points inside cube back to base

oX2in = []
oY2in = []
oZ2in = []

for i in range(0,len(X2in)):
    Pt = [X2in[i]-C[0],Y2in[i]-C[1],Z2in[i]-C[2]]
    Pto = np.matmul(Rinv,Pt)
    oX2in.append(Pto[0])
    oY2in.append(Pto[1])
    oZ2in.append(Pto[2])

oX1in = []
oY1in = []
oZ1in = []

for i in range(0,len(X1in)):
    Pt = [X1in[i]-C[0],Y1in[i]-C[1],Z1in[i]-C[2]]
    Pto = np.matmul(Rinv,Pt)
    oX1in.append(Pto[0])
    oY1in.append(Pto[1])
    oZ1in.append(Pto[2])

print('Total majority component:')
print(len(oX1in))
print('Total minority component:')
print(len(oX2in))
print('Total')
print(len(oX1in)+len(oX2in))


fig = plt.figure()
ax= fig.add_subplot(1,2,1,projection='3d', proj_type = 'ortho')
ax.plot(P1P2x,P1P2y,P1P2z,'k-')
ax.plot(P1P3x,P1P3y,P1P3z,'k-')
ax.plot(P1P5x,P1P5y,P1P5z,'k-')
ax.plot(P2P4x,P2P4y,P2P4z,'k-')
ax.plot(P2P6x,P2P6y,P2P6z,'k-')
ax.plot(P3P4x,P3P4y,P3P4z,'k-')
ax.plot(P3P7x,P3P7y,P3P7z,'k-')
ax.plot(P4P8x,P4P8y,P4P8z,'k-')
ax.plot(P7P8x,P7P8y,P7P8z,'k-')
ax.plot(P6P8x,P6P8y,P6P8z,'k-')
ax.plot(P5P6x,P5P6y,P5P6z,'k-')
ax.plot(P5P7x,P5P7y,P5P7z,'k-')
ax.scatter(X2in,Y2in,Z2in,'b')
#ax.scatter(X1in,Y1in,Z1in,'r')   
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')               
ax = fig.add_subplot(1,2,2,projection='3d', proj_type = 'ortho')
ax.plot(oP1P2x,oP1P2y,oP1P2z,'k-')
ax.plot(oP1P3x,oP1P3y,oP1P3z,'k-')
ax.plot(oP1P5x,oP1P5y,oP1P5z,'k-')
ax.plot(oP2P4x,oP2P4y,oP2P4z,'k-')
ax.plot(oP2P6x,oP2P6y,oP2P6z,'k-')
ax.plot(oP3P4x,oP3P4y,oP3P4z,'k-')
ax.plot(oP3P7x,oP3P7y,oP3P7z,'k-')
ax.plot(oP4P8x,oP4P8y,oP4P8z,'k-')
ax.plot(oP7P8x,oP7P8y,oP7P8z,'k-')
ax.plot(oP6P8x,oP6P8y,oP6P8z,'k-')
ax.plot(oP5P6x,oP5P6y,oP5P6z,'k-')
ax.plot(oP5P7x,oP5P7y,oP5P7z,'k-')
ax.scatter(oX2in,oY2in,oZ2in,'b')
#ax.scatter(oX1in,oY1in,oZ1in,'r')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()   

# Plot points with periodic BC
aX1  = []
aX2  = []
aX3  = []
aX4  = []
aX5  = []
aX6  = []
aX7  = []
aX8  = []
aX9  = []
aX10 = []
aX11 = []
aX12 = []
aX13 = []
aX14 = []
aX15 = []
aX16 = []
aX17 = []
aX18 = []
aX19 = []
aX20 = []
aX21 = []
aX22 = []
aX23 = []
aX24 = []
aX25 = []
aX26 = []
aX27 = []

aY1  = []
aY2  = []
aY3  = []
aY4  = []
aY5  = []
aY6  = []
aY7  = []
aY8  = []
aY9  = []
aY10 = []
aY11 = []
aY12 = []
aY13 = []
aY14 = []
aY15 = []
aY16 = []
aY17 = []
aY18 = []
aY19 = []
aY20 = []
aY21 = []
aY22 = []
aY23 = []
aY24 = []
aY25 = []
aY26 = []
aY27 = []

aZ1  = []
aZ2  = []
aZ3  = []
aZ4  = []
aZ5  = []
aZ6  = []
aZ7  = []
aZ8  = []
aZ9  = []
aZ10 = []
aZ11 = []
aZ12 = []
aZ13 = []
aZ14 = []
aZ15 = []
aZ16 = []
aZ17 = []
aZ18 = []
aZ19 = []
aZ20 = []
aZ21 = []
aZ22 = []
aZ23 = []
aZ24 = []
aZ25 = []
aZ26 = []
aZ27 = []


bX1  = []
bX2  = []
bX3  = []
bX4  = []
bX5  = []
bX6  = []
bX7  = []
bX8  = []
bX9  = []
bX10 = []
bX11 = []
bX12 = []
bX13 = []
bX14 = []
bX15 = []
bX16 = []
bX17 = []
bX18 = []
bX19 = []
bX20 = []
bX21 = []
bX22 = []
bX23 = []
bX24 = []
bX25 = []
bX26 = []
bX27 = []

bY1  = []
bY2  = []
bY3  = []
bY4  = []
bY5  = []
bY6  = []
bY7  = []
bY8  = []
bY9  = []
bY10 = []
bY11 = []
bY12 = []
bY13 = []
bY14 = []
bY15 = []
bY16 = []
bY17 = []
bY18 = []
bY19 = []
bY20 = []
bY21 = []
bY22 = []
bY23 = []
bY24 = []
bY25 = []
bY26 = []
bY27 = []

bZ1  = []
bZ2  = []
bZ3  = []
bZ4  = []
bZ5  = []
bZ6  = []
bZ7  = []
bZ8  = []
bZ9  = []
bZ10 = []
bZ11 = []
bZ12 = []
bZ13 = []
bZ14 = []
bZ15 = []
bZ16 = []
bZ17 = []
bZ18 = []
bZ19 = []
bZ20 = []
bZ21 = []
bZ22 = []
bZ23 = []
bZ24 = []
bZ25 = []
bZ26 = []
bZ27 = []

for i in range(0,len(X2in)):
    bX1.append(oX2in[i]+Lx)   
    bX2.append(oX2in[i]+Lx)
    bX3.append(oX2in[i]+Lx)
    bX4.append(oX2in[i])
    bX5.append(oX2in[i])
    bX6.append(oX2in[i])
    bX7.append(oX2in[i]-Lx)
    bX8.append(oX2in[i]-Lx)
    bX9.append(oX2in[i]-Lx)
    bX10.append(oX2in[i]+Lx)
    bX11.append(oX2in[i]+Lx)
    bX12.append(oX2in[i]+Lx)
    bX13.append(oX2in[i])
    bX14.append(oX2in[i])
    bX15.append(oX2in[i])
    bX16.append(oX2in[i]-Lx)
    bX17.append(oX2in[i]-Lx)
    bX18.append(oX2in[i]-Lx)
    bX19.append(oX2in[i]+Lx)
    bX20.append(oX2in[i]+Lx)
    bX21.append(oX2in[i]+Lx)
    bX22.append(oX2in[i])
    bX23.append(oX2in[i])
    bX24.append(oX2in[i])
    bX25.append(oX2in[i]-Lx)
    bX26.append(oX2in[i]-Lx)
    bX27.append(oX2in[i]-Lx)
    bY1.append(oY2in[i]+Ly)
    bY2.append(oY2in[i])
    bY3.append(oY2in[i]-Ly)
    bY4.append(oY2in[i]+Ly)
    bY5.append(oY2in[i])
    bY6.append(oY2in[i]-Ly)
    bY7.append(oY2in[i]+Ly)
    bY8.append(oY2in[i])
    bY9.append(oY2in[i]-Ly)
    bY10.append(oY2in[i]+Ly)
    bY11.append(oY2in[i])
    bY12.append(oY2in[i]-Ly)
    bY13.append(oY2in[i]+Ly)
    bY14.append(oY2in[i])
    bY15.append(oY2in[i]-Ly)
    bY16.append(oY2in[i]+Ly)
    bY17.append(oY2in[i])
    bY18.append(oY2in[i]-Ly)
    bY19.append(oY2in[i]+Ly)
    bY20.append(oY2in[i])
    bY21.append(oY2in[i]-Ly)
    bY22.append(oY2in[i]+Ly)
    bY23.append(oY2in[i])
    bY24.append(oY2in[i]-Ly)
    bY25.append(oY2in[i]+Ly)
    bY26.append(oY2in[i])
    bY27.append(oY2in[i]-Ly)
    bZ1.append(oZ2in[i]-Lz)
    bZ2.append(oZ2in[i]-Lz)
    bZ3.append(oZ2in[i]-Lz)
    bZ4.append(oZ2in[i]-Lz)
    bZ5.append(oZ2in[i]-Lz)
    bZ6.append(oZ2in[i]-Lz)
    bZ7.append(oZ2in[i]-Lz)
    bZ8.append(oZ2in[i]-Lz)
    bZ9.append(oZ2in[i]-Lz)
    bZ10.append(oZ2in[i])
    bZ11.append(oZ2in[i])
    bZ12.append(oZ2in[i])
    bZ13.append(oZ2in[i])
    bZ14.append(oZ2in[i])
    bZ15.append(oZ2in[i])
    bZ16.append(oZ2in[i])
    bZ17.append(oZ2in[i])
    bZ18.append(oZ2in[i])
    bZ19.append(oZ2in[i]+Lz)
    bZ20.append(oZ2in[i]+Lz)
    bZ21.append(oZ2in[i]+Lz)
    bZ22.append(oZ2in[i]+Lz)
    bZ23.append(oZ2in[i]+Lz)
    bZ24.append(oZ2in[i]+Lz)
    bZ25.append(oZ2in[i]+Lz)
    bZ26.append(oZ2in[i]+Lz)
    bZ27.append(oZ2in[i]+Lz)
    
    
for i in range(0,len(X1in)):
    aX1.append(oX1in[i]+Lx)   
    aX2.append(oX1in[i]+Lx)
    aX3.append(oX1in[i]+Lx)
    aX4.append(oX1in[i])
    aX5.append(oX1in[i])
    aX6.append(oX1in[i])
    aX7.append(oX1in[i]-Lx)
    aX8.append(oX1in[i]-Lx)
    aX9.append(oX1in[i]-Lx)
    aX10.append(oX1in[i]+Lx)
    aX11.append(oX1in[i]+Lx)
    aX12.append(oX1in[i]+Lx)
    aX13.append(oX1in[i])
    aX14.append(oX1in[i])
    aX15.append(oX1in[i])
    aX16.append(oX1in[i]-Lx)
    aX17.append(oX1in[i]-Lx)
    aX18.append(oX1in[i]-Lx)
    aX19.append(oX1in[i]+Lx)
    aX20.append(oX1in[i]+Lx)
    aX21.append(oX1in[i]+Lx)
    aX22.append(oX1in[i])
    aX23.append(oX1in[i])
    aX24.append(oX1in[i])
    aX25.append(oX1in[i]-Lx)
    aX26.append(oX1in[i]-Lx)
    aX27.append(oX1in[i]-Lx)
    aY1.append(oY1in[i]+Ly)
    aY2.append(oY1in[i])
    aY3.append(oY1in[i]-Ly)
    aY4.append(oY1in[i]+Ly)
    aY5.append(oY1in[i])
    aY6.append(oY1in[i]-Ly)
    aY7.append(oY1in[i]+Ly)
    aY8.append(oY1in[i])
    aY9.append(oY1in[i]-Ly)
    aY10.append(oY1in[i]+Ly)
    aY11.append(oY1in[i])
    aY12.append(oY1in[i]-Ly)
    aY13.append(oY1in[i]+Ly)
    aY14.append(oY1in[i])
    aY15.append(oY1in[i]-Ly)
    aY16.append(oY1in[i]+Ly)
    aY17.append(oY1in[i])
    aY18.append(oY1in[i]-Ly)
    aY19.append(oY1in[i]+Ly)
    aY20.append(oY1in[i])
    aY21.append(oY1in[i]-Ly)
    aY22.append(oY1in[i]+Ly)
    aY23.append(oY1in[i])
    aY24.append(oY1in[i]-Ly)
    aY25.append(oY1in[i]+Ly)
    aY26.append(oY1in[i])
    aY27.append(oY1in[i]-Ly)
    aZ1.append(oZ1in[i]-Lz)
    aZ2.append(oZ1in[i]-Lz)
    aZ3.append(oZ1in[i]-Lz)
    aZ4.append(oZ1in[i]-Lz)
    aZ5.append(oZ1in[i]-Lz)
    aZ6.append(oZ1in[i]-Lz)
    aZ7.append(oZ1in[i]-Lz)
    aZ8.append(oZ1in[i]-Lz)
    aZ9.append(oZ1in[i]-Lz)
    aZ10.append(oZ1in[i])
    aZ11.append(oZ1in[i])
    aZ12.append(oZ1in[i])
    aZ13.append(oZ1in[i])
    aZ14.append(oZ1in[i])
    aZ15.append(oZ1in[i])
    aZ16.append(oZ1in[i])
    aZ17.append(oZ1in[i])
    aZ18.append(oZ1in[i])
    aZ19.append(oZ1in[i]+Lz)
    aZ20.append(oZ1in[i]+Lz)
    aZ21.append(oZ1in[i]+Lz)
    aZ22.append(oZ1in[i]+Lz)
    aZ23.append(oZ1in[i]+Lz)
    aZ24.append(oZ1in[i]+Lz)
    aZ25.append(oZ1in[i]+Lz)
    aZ26.append(oZ1in[i]+Lz)
    aZ27.append(oZ1in[i]+Lz)    
    
fig = plt.figure()
ax= fig.add_subplot(projection='3d', proj_type = 'ortho')
ax.plot(oP1P2x,oP1P2y,oP1P2z,'k-')
ax.plot(oP1P3x,oP1P3y,oP1P3z,'k-')
ax.plot(oP1P5x,oP1P5y,oP1P5z,'k-')
ax.plot(oP2P4x,oP2P4y,oP2P4z,'k-')
ax.plot(oP2P6x,oP2P6y,oP2P6z,'k-')
ax.plot(oP3P4x,oP3P4y,oP3P4z,'k-')
ax.plot(oP3P7x,oP3P7y,oP3P7z,'k-')
ax.plot(oP4P8x,oP4P8y,oP4P8z,'k-')
ax.plot(oP7P8x,oP7P8y,oP7P8z,'k-')
ax.plot(oP6P8x,oP6P8y,oP6P8z,'k-')
ax.plot(oP5P6x,oP5P6y,oP5P6z,'k-')
ax.plot(oP5P7x,oP5P7y,oP5P7z,'k-')
#ax.scatter(aX1,aY1,aZ1,color='r')
#ax.scatter(aX2,aY2,aZ2,color='r')
#ax.scatter(aX3,aY3,aZ3,color='r')
#ax.scatter(aX4,aY4,aZ4,color='r')
#ax.scatter(aX5,aY5,aZ5,color='r')
#ax.scatter(aX6,aY6,aZ6,color='r')
#ax.scatter(aX7,aY7,aZ7,color='r')
#ax.scatter(aX8,aY8,aZ8,color='r')
#ax.scatter(aX9,aY9,aZ9,color='r')
#ax.scatter(aX10,aY10,aZ10,color='r')
#ax.scatter(aX11,aY11,aZ11,color='r')
#ax.scatter(aX12,aY12,aZ12,color='r')
#ax.scatter(aX13,aY13,aZ13,color='r')
#ax.scatter(aX14,aY14,aZ14,color='r')
#ax.scatter(aX15,aY15,aZ15,color='r')
#ax.scatter(aX16,aY16,aZ16,color='r')
#ax.scatter(aX17,aY17,aZ17,color='r')
#ax.scatter(aX18,aY18,aZ18,color='r')
#ax.scatter(aX19,aY19,aZ19,color='r')
#ax.scatter(aX20,aY20,aZ20,color='r')
#ax.scatter(aX21,aY21,aZ21,color='r')
#ax.scatter(aX22,aY22,aZ22,color='r')
#ax.scatter(aX23,aY23,aZ23,color='r')
#ax.scatter(aX24,aY24,aZ24,color='r')
#ax.scatter(aX25,aY25,aZ25,color='r')
#ax.scatter(aX26,aY26,aZ26,color='r')
#ax.scatter(aX27,aY27,aZ27,color='r') 
ax.scatter(bX1,bY1,bZ1,color='b')
ax.scatter(bX2,bY2,bZ2,color='b')
ax.scatter(bX3,bY3,bZ3,color='b')
ax.scatter(bX4,bY4,bZ4,color='b')
ax.scatter(bX5,bY5,bZ5,color='b')
ax.scatter(bX6,bY6,bZ6,color='b')
ax.scatter(bX7,bY7,bZ7,color='b')
ax.scatter(bX8,bY8,bZ8,color='b')
ax.scatter(bX9,bY9,bZ9,color='b')
ax.scatter(bX10,bY10,bZ10,color='b')
ax.scatter(bX11,bY11,bZ11,color='b')
ax.scatter(bX12,bY12,bZ12,color='b')
ax.scatter(bX13,bY13,bZ13,color='b')
ax.scatter(bX14,bY14,bZ14,color='b')
ax.scatter(bX15,bY15,bZ15,color='b')
ax.scatter(bX16,bY16,bZ16,color='b')
ax.scatter(bX17,bY17,bZ17,color='b')
ax.scatter(bX18,bY18,bZ18,color='b')
ax.scatter(bX19,bY19,bZ19,color='b')
ax.scatter(bX20,bY20,bZ20,color='b')
ax.scatter(bX21,bY21,bZ21,color='b')
ax.scatter(bX22,bY22,bZ22,color='b')
ax.scatter(bX23,bY23,bZ23,color='b')
ax.scatter(bX24,bY24,bZ24,color='b')
ax.scatter(bX25,bY25,bZ25,color='b')
ax.scatter(bX26,bY26,bZ26,color='b')
ax.scatter(bX27,bY27,bZ27,color='b') 
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()       
    
    
# Output X Y and Z of new particles

F = open('gyr.data','w')
F.write('LAMMPS Data File \n')
F.write('\n')
F.write(str(len(oX1in)+len(oX2in))+' atoms\n')
F.write('2 atom types \n')
F.write('\n')
F.write(str(-Lx/2) + ' ' + str(Lx/2) + ' xlo xhi \n')
F.write(str(-Ly/2) + ' ' + str(Ly/2) + ' ylo yhi \n')
F.write(str(-Lz/2) + ' ' + str(Lz/2) + ' zlo zhi \n')
F.write('\n')
F.write('Masses\n')
F.write('\n')
F.write('1 1\n')
F.write('2 1\n')
F.write('\n')
F.write('Atoms # atomic\n')
F.write('\n')
count = 1
for i in range(0,len(oX1in)):
    F.write(str(count) + ' 1 ' + str(oX1in[i]) + ' ' + str(oY1in[i]) + ' ' + str(oZ1in[i]) +'\n')
    count = count + 1
for i in range(0,len(oX2in)):
    F.write(str(count) + ' 2 ' + str(oX2in[i]) + ' ' + str(oY2in[i]) + ' ' + str(oZ2in[i]) +'\n')
    count = count + 1    
F.close()