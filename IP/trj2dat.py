import math as m
import random as r
import sys as s
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
from pyscal.trajectory import Trajectory

NTRJ = int(s.argv[1])   # Timestep chosen for data

NpolyS = int(s.argv[2]) # Number of long polymer chains in each particle
NpolyL = int(s.argv[3]) # Number of short polymer chains in each particle
NspaceL = int(s.argv[4]) # Number of spacer beads in each long polymer
NspaceS = int(s.argv[5]) # Number of spacer beads in each  short polymer

#Variables

phi = m.pi*(3-m.sqrt(5))

Nbig = 100  # Number of spheres in each core

NinterS = 2 # Number of interactive beads in each short polymer 
NinterL = 4 # Number of interactive beads in each long polymer

D = 0.8 # Distance between linker and interactive/flanking bead
CRad = 5 # Radius of Core Particle
# Set types

AtTN = 10 # Number of Atom Types

Cor1 = 1 # Core Particle for Particle 1
SpaT = 2 # Spacer Particle
P1STA = [5, 3] # Particle 1 Short Polymer Interactive, Length equal NinterS
P1LTA = [3, 5, 4, 6] # Particle 1 Long Polymer Interactive, Length equal NinterL
P1STB = [3, 5] # Particle 2 Short Polymer Interactive, Length equal NinterS
P1LTB = [5, 3, 6, 4] # Particle 2 Long Polymer Interactive, Length equal NinterL
FlaT = 7 # Flanking Particle
Cor2 = 8 # Core Particle for Particle 2
TraA = 9  # Core particle for easier COM tracking
TraB = 10 # Core particle for easier COM tracking



BTN = 4 # Number of Bond Types
BT1 = 1 # Spacer Bonds
BT2 = 2 # Base Bonds
BT3 = 3 # Flanking Bonds
BT4 = 4 # Flanking-Linker Bonds

ATN = 3 # Number of Angle Types
AT1 = 1 # Spacer Angles
AT2 = 2 # Base Angles
AT3 = 3 # Flanking Angles

omega = 1   # Size of beads

Cscale = 1 # Scaling on core
Pscale = 0.125/10 # Scaling on polymer

GScale = 1 # Global scaling 

Mass = 1.0   # Mass of beads


OUTNAME = 'LargeOUT.data'

# Read data file





traj = Trajectory('PNP.lammpstrj')


SL = traj[NTRJ].to_dict()


BOXXL = SL[0]["box"][0][0]
BOXXH = SL[0]["box"][0][1]

BOXYL = SL[0]["box"][1][0]
BOXYH = SL[0]["box"][1][1]

BOXZL = SL[0]["box"][2][0]
BOXZH = SL[0]["box"][2][1]


ID = SL[0]["atoms"]["id"]
XS = SL[0]["atoms"]["xs"]
YS = SL[0]["atoms"]["ys"]
ZS = SL[0]["atoms"]["zs"]
TYPT = SL[0]["atoms"]["type"]


XI = []
YI = []
ZI = []

TYP = []

for i in range(0,len(XS)):
	if TYPT[i] == 9:
		TYP.append(1)
		XI.append(XS[i]*(BOXXH-BOXXL) + BOXXL)
		YI.append(YS[i]*(BOXYH-BOXYL) + BOXYL)
		ZI.append(ZS[i]*(BOXZH-BOXZL) + BOXZL)	
	elif TYPT[i] == 10:
		TYP.append(2)
		XI.append(XS[i]*(BOXXH-BOXXL) + BOXXL)
		YI.append(YS[i]*(BOXYH-BOXYL) + BOXYL)
		ZI.append(ZS[i]*(BOXZH-BOXZL) + BOXZL)	



Atom = []
Bond = []
Angl = []



### Make Cores ###
Core = []

## Make points in unit sphere ##
for i in range(0,Nbig):        
    # Choose point in sphere

    Ys = 1 - (i/(float(Nbig)-1))*2
   
    radius = m.sqrt(1-Ys**2)
    
    theta = phi*i
    Xs = m.cos(theta)*radius
    Zs = m.sin(theta)*radius
    Core.append([CRad*Xs,CRad*Ys,CRad*Zs])

# Find maximum distance between points, = diameter
Dmax = -1
for i in range(0,Nbig):
    for j in range(i+1,Nbig):
        Dnew = m.sqrt((Core[i][0]-Core[j][0])**2+(Core[i][1]-Core[j][1])**2+(Core[i][2]-Core[j][2])**2)
        if Dnew > Dmax:
            Imax = i
            Jmax = j
            Dmax = Dnew
print('Core radius is ' + str(Dmax))



CoreL = []
CoreS = []
CoreST = []

for i in range(0,NpolyL):        
    Ys = 1 - (i/(float(NpolyL)-1))*2
   
    radius = m.sqrt(1-Ys**2)
    
    theta = phi*i
    Xs = m.cos(theta)*radius
    Zs = m.sin(theta)*radius
    CoreL.append([CRad*(1+Pscale)*Xs,CRad*(1+Pscale)*Ys,CRad*(1+Pscale)*Zs])


Srot = np.linspace(-m.pi,m.pi,50)
Drot = []

for k in Srot:
	for i in range(0,NpolyS):        	
		Xs = 1 - (i/(float(NpolyS)-1))*2
	
		radius = m.sqrt(1-Xs**2)
		
		theta = phi*i + k
		Ys = m.cos(theta)*radius
		Zs = m.sin(theta)*radius
		CoreST.append([CRad*(1+Pscale)*Xs,CRad*(1+Pscale)*Ys,CRad*(1+Pscale)*Zs])
	Dmin = 10000
	for i in CoreL:
		for j in CoreST:
			Dnew = m.sqrt((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)
			if Dnew < Dmin:
				Dmin = Dnew
	
	Drot.append(Dmin)


SrotF = Srot[Drot.index(max(Drot))]
print('Mimimum distance between short and long polymer core beads is ' + str(max(Drot)))

for i in range(0,NpolyS):        	
		Xs = 1 - (i/(float(NpolyS)-1))*2
	
		radius = m.sqrt(1-Xs**2)
		
		theta = phi*i + SrotF
		Ys = m.cos(theta)*radius
		Zs = m.sin(theta)*radius
		CoreS.append([CRad*(1+Pscale)*Xs,CRad*(1+Pscale)*Ys,CRad*(1+Pscale)*Zs])


## Make Atoms of Core 1 and 2 ##
# Atom = [AtID, MolID, AType, X, Y ,Z]  For Lammps it is [AtID, MolID, AType, Charge, X, Y ,Z] but all charges here are 0 so done below
# Bond = [BondType, AtID1, AtID2] same for LAMMPS
# Angle = [AngleType, AtID1, AtID2, AtID3] same for LAMMPS



GAtomC = 0


for k in range(0,len(TYP)):
	C1S = []
	C1L = []
	C2S = []
	C2L = []

	# Make particle core MolID k + 2
	if TYP[k] == 1:
		GAtomC = GAtomC + 1
		Atom.append([GAtomC,k+2,TraA,GScale*XI[k],GScale*YI[k],GScale*ZI[k]]) 
		for i in Core:
			GAtomC = GAtomC + 1
			Atom.append([GAtomC,k+2,Cor1,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]]) 
		
		for i in CoreL:
			GAtomC = GAtomC + 1
			C1L.append(GAtomC-1)
			Atom.append([GAtomC,k+2,Cor1,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]]) 
		
		for i in CoreS:
			GAtomC = GAtomC + 1
			C1S.append(GAtomC-1)
			Atom.append([GAtomC,k+2,Cor1,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]]) 
	
	elif TYP[k] == 2:	
		GAtomC = GAtomC + 1
		Atom.append([GAtomC,k+2,TraB,GScale*XI[k],GScale*YI[k],GScale*ZI[k]]) 	
		for i in Core:
			GAtomC = GAtomC + 1
			Atom.append([GAtomC,k+2,Cor2,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]]) 
		
		for i in CoreL:
			GAtomC = GAtomC + 1
			C1L.append(GAtomC-1)
			Atom.append([GAtomC,k+2,Cor2,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]]) 
		
		for i in CoreS:
			GAtomC = GAtomC + 1
			C1S.append(GAtomC-1)
			Atom.append([GAtomC,k+2,Cor2,i[0] + GScale*XI[k],i[1]+ GScale*YI[k],i[2]+ GScale*ZI[k]])

	# Make polymer chains MolID 1
	MolID = 1
	# Choose Atom 
	
	
	if TYP[k] == 1:
	
	
		for i in C1S:
		
			BaseA = Atom[i][0]
			BaseX = Atom[i][3] - GScale*XI[k]
			BaseY = Atom[i][4] - GScale*YI[k]
			BaseZ = Atom[i][5] - GScale*ZI[k]
			
			BL = m.sqrt(BaseX**2+BaseY**2+BaseZ**2)
			B = [BaseX/BL,BaseY/BL,BaseZ/BL]
			
			
			# Compute perpendicular vectors for interactive and flanking particle
			
			if abs(BaseX) > 1e-3:
				I2 = r.random()
				I3 = r.random()
				I1 = -(I2*BaseY+I3*BaseZ)/(BaseX)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			elif abs(BaseZ) > 1e-3:
				I1 = r.random()
				I2 = r.random()
				I3 = -(I1*BaseX+I2*BaseY)/(BaseZ)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL  
			
			else:
				I1 = r.random()
				I3 = r.random()
				I2 = -(I1*BaseX+I3*BaseZ)/(BaseY)
			
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			I = [I1,I2,I3]
		
			F1 = np.cross(I,B)
			F2 = np.cross(B,I)
		
		
		
			for j in range(1,NspaceS+1):
				GAtomC = GAtomC + 1
				NewX = (j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (j*Pscale+1)*BaseZ + GScale*ZI[k]
				if j == 1:
					Bond.append([BT1,BaseA,GAtomC,' #C1Sspaj0'])
				else:
					Bond.append([BT1,GAtomC-1,GAtomC,' #C1Sspaj'])
					if j != 2:
						Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C1Sspa'])		
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
			IC = []
			FL1 = []
			FL2 = []    
			for j in range(1,NinterS+1):
				GAtomC = GAtomC + 1
				Link = GAtomC
				NewX = (NspaceS*Pscale+j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (NspaceS*Pscale+j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (NspaceS*Pscale+j*Pscale+1)*BaseZ + GScale*ZI[k]
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
				if j == 1:
					Bond.append([BT2,GAtomC-1,GAtomC,' #C1Sj0'])
					Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C1Sj0'])
				elif j == 2 :		
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Sj1'])
					Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C1Sj1'])
				else:
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Sj'])
					Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C1Sj'])			
				
				GAtomC = GAtomC + 1
				Inter = GAtomC	
				IC.append(GAtomC)
				NIX = D*I1 + NewX
				NIY = D*I2 + NewY
				NIZ = D*I3 + NewZ
				Atom.append([GAtomC,MolID,P1STA[j-1],NIX,NIY,NIZ])
				Bond.append([BT3,Link,Inter,' #C1SI'])
				
				GAtomC = GAtomC + 1
				Flank1 = GAtomC
				FL1.append(Flank1)		
				NF1X = D*F1[0] + NewX + D*I1
				NF1Y = D*F1[1] + NewY + D*I2
				NF1Z = D*F1[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
				Bond.append([BT3,Inter,Flank1,' #C1SF1'])      
				Bond.append([BT4,Link,Flank1,'#C1SLF1'])
				
				GAtomC = GAtomC + 1
				Flank2 = GAtomC
				FL2.append(Flank2)        
				NF2X = D*F2[0] + NewX + D*I1
				NF2Y = D*F2[1] + NewY + D*I2
				NF2Z = D*F2[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
				Bond.append([BT3,Inter,Flank2,' #C1SF2'])
				Bond.append([BT4,Link,Flank2,' #C1SLF2'])		
				Angl.append([AT3,Flank1,Inter,Flank2,' #C1SF'])
				
				if j == NinterS:
					GAtomC = GAtomC + 1
					FlankF = GAtomC			
					NFX = D*I1 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseX + GScale*XI[k]
					NFY = D*I2 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseY + GScale*YI[k]
					NFZ = D*I3 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseZ + GScale*ZI[k]
					Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
					Bond.append([BT3,Inter,FlankF,' #C1SFF'])
					Bond.append([BT4,Link,FlankF,' #C1SLFF'])			
					Angl.append([AT3,IC[-2],Inter,FlankF,' #C1SFF'])  
		
			Bond.append([BT3,IC[0],IC[1],' #C1SIB1'])		
		
			for n in range(0,NinterS-2):
				Angl.append([AT3,IC[n],IC[n+1],IC[n+2],' #C1SI'])
				Bond.append([BT3,IC[n+1],IC[n+2],' #C1SIB2'])		
		
				
		for i in C1L:
			BaseA = Atom[i][0]
			BaseX = Atom[i][3] - GScale*XI[k]
			BaseY = Atom[i][4] - GScale*YI[k]
			BaseZ = Atom[i][5] - GScale*ZI[k]
			
			BL = m.sqrt(BaseX**2+BaseY**2+BaseZ**2)
			B = [BaseX/BL,BaseY/BL,BaseZ/BL]
			
			
			# Compute perpendicular vectors for interactive and flanking particle
			
			if abs(BaseX) > 1e-3:
				I2 = r.random()
				I3 = r.random()
				I1 = -(I2*BaseY+I3*BaseZ)/(BaseX)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			elif abs(BaseZ) > 1e-3:
				I1 = r.random()
				I2 = r.random()
				I3 = -(I1*BaseX+I2*BaseY)/(BaseZ)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL  
			
			else:
				I1 = r.random()
				I3 = r.random()
				I2 = -(I1*BaseX+I3*BaseZ)/(BaseY)
			
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			I = [I1,I2,I3]
		
			F1 = np.cross(I,B)
			F2 = np.cross(B,I)    
			
			for j in range(1,NspaceL+1):
				GAtomC = GAtomC + 1
				NewX = (j*Pscale+1)*BaseX  + GScale*XI[k]
				NewY = (j*Pscale+1)*BaseY  + GScale*YI[k]
				NewZ = (j*Pscale+1)*BaseZ  + GScale*ZI[k]
				if j == 1:
					Bond.append([BT1,BaseA,GAtomC,' #C1Lspaj0'])
				else:
					Bond.append([BT1,GAtomC-1,GAtomC,' #C1Lspaj'])
					if j != 2:
						Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C1Lspa'])				
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
			IC = []
			FL1 = []
			FL2 = []     
			for j in range(1,NinterL+1):
				GAtomC = GAtomC + 1
				Link = GAtomC		
				NewX = (NspaceL*Pscale+j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (NspaceL*Pscale+j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (NspaceL*Pscale+j*Pscale+1)*BaseZ + GScale*ZI[k]
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
				if j == 1:
					Bond.append([BT2,GAtomC-1,GAtomC,' #C1Lj0'])
					Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C1Lj0'])
				elif j == 2 :		
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Lj1'])
					Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C1Lj1'])
				else:
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1LSj'])
					Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C1Lj'])		
		
				GAtomC = GAtomC + 1	
				Inter = GAtomC			
				IC.append(GAtomC)
				NIX = D*I1 + NewX
				NIY = D*I2 + NewY
				NIZ = D*I3 + NewZ
				Atom.append([GAtomC,MolID,P1LTA[j-1],NIX,NIY,NIZ])
				Bond.append([BT3,Link,Inter,' #C1LI'])
				
				GAtomC = GAtomC + 1    
				Flank1 = GAtomC
				FL1.append(Flank1)        
				NF1X = D*F1[0] + NewX + D*I1
				NF1Y = D*F1[1] + NewY + D*I2
				NF1Z = D*F1[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
				Bond.append([BT3,Inter,Flank1,' #C1LF1'])      
				Bond.append([BT4,Link,Flank1,'#C1LLF1'])  
				
				GAtomC = GAtomC + 1
				Flank2 = GAtomC
				FL2.append(Flank2) 		
				NF2X = D*F2[0] + NewX + D*I1
				NF2Y = D*F2[1] + NewY + D*I2
				NF2Z = D*F2[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
				Bond.append([BT3,Inter,Flank2,' #C1LF2'])
				Bond.append([BT4,Link,Flank2,' #C1LLF2'])		
				Angl.append([AT3,Flank1,Inter,Flank2,' #C1LF'])
				
				if j == NinterL:
					GAtomC = GAtomC + 1
					FlankF = GAtomC			
					NFX = D*I1 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseX + GScale*XI[k]
					NFY = D*I2 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseY + GScale*YI[k]
					NFZ = D*I3 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseZ + GScale*ZI[k]
					Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
					Bond.append([BT3,Inter,FlankF,' #C1LFF'])
					Bond.append([BT4,Link,FlankF,' #C1LLFF'])			
					Angl.append([AT3,IC[-2],Inter,FlankF,' #C1LFF'])  
		
			Bond.append([BT3,IC[0],IC[1],' #C1LIB1'])
		
			for n in range(0,NinterL-2):
				Angl.append([AT3,IC[n],IC[n+1],IC[n+2],' #C1LI'])
				Bond.append([BT3,IC[n+1],IC[n+2],' #C1LIB2'])		        
	
	
	elif TYP[k] == 2:
	
		for i in C1S:
		
			BaseA = Atom[i][0]
			BaseX = Atom[i][3] - GScale*XI[k]
			BaseY = Atom[i][4] - GScale*YI[k]
			BaseZ = Atom[i][5] - GScale*ZI[k]
			
			BL = m.sqrt(BaseX**2+BaseY**2+BaseZ**2)
			B = [BaseX/BL,BaseY/BL,BaseZ/BL]
			
			
			# Compute perpendicular vectors for interactive and flanking particle
			
			if abs(BaseX) > 1e-3:
				I2 = r.random()
				I3 = r.random()
				I1 = -(I2*BaseY+I3*BaseZ)/(BaseX)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			elif abs(BaseZ) > 1e-3:
				I1 = r.random()
				I2 = r.random()
				I3 = -(I1*BaseX+I2*BaseY)/(BaseZ)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL  
			
			else:
				I1 = r.random()
				I3 = r.random()
				I2 = -(I1*BaseX+I3*BaseZ)/(BaseY)
			
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			I = [I1,I2,I3]
		
			F1 = np.cross(I,B)
			F2 = np.cross(B,I)
		
		
		
			for j in range(1,NspaceS+1):
				GAtomC = GAtomC + 1
				NewX = (j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (j*Pscale+1)*BaseZ + GScale*ZI[k]
				if j == 1:
					Bond.append([BT1,BaseA,GAtomC,' #C1Sspaj0'])
				else:
					Bond.append([BT1,GAtomC-1,GAtomC,' #C1Sspaj'])
					if j != 2:
						Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C1Sspa'])		
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
			IC = []
			FL1 = []
			FL2 = []    
			for j in range(1,NinterS+1):
				GAtomC = GAtomC + 1
				Link = GAtomC
				NewX = (NspaceS*Pscale+j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (NspaceS*Pscale+j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (NspaceS*Pscale+j*Pscale+1)*BaseZ + GScale*ZI[k]
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
				if j == 1:
					Bond.append([BT2,GAtomC-1,GAtomC,' #C1Sj0'])
					Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C1Sj0'])
				elif j == 2 :		
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Sj1'])
					Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C1Sj1'])
				else:
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Sj'])
					Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C1Sj'])			
				
				GAtomC = GAtomC + 1
				Inter = GAtomC	
				IC.append(GAtomC)
				NIX = D*I1 + NewX
				NIY = D*I2 + NewY
				NIZ = D*I3 + NewZ
				Atom.append([GAtomC,MolID,P1STB[j-1],NIX,NIY,NIZ])
				Bond.append([BT3,Link,Inter,' #C1SI'])
				
				GAtomC = GAtomC + 1
				Flank1 = GAtomC
				FL1.append(Flank1)		
				NF1X = D*F1[0] + NewX + D*I1
				NF1Y = D*F1[1] + NewY + D*I2
				NF1Z = D*F1[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
				Bond.append([BT3,Inter,Flank1,' #C1SF1'])      
				Bond.append([BT4,Link,Flank1,'#C1SLF1'])
				
				GAtomC = GAtomC + 1
				Flank2 = GAtomC
				FL2.append(Flank2)        
				NF2X = D*F2[0] + NewX + D*I1
				NF2Y = D*F2[1] + NewY + D*I2
				NF2Z = D*F2[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
				Bond.append([BT3,Inter,Flank2,' #C1SF2'])
				Bond.append([BT4,Link,Flank2,' #C1SLF2'])		
				Angl.append([AT3,Flank1,Inter,Flank2,' #C1SF'])
				
				if j == NinterS:
					GAtomC = GAtomC + 1
					FlankF = GAtomC			
					NFX = D*I1 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseX + GScale*XI[k]
					NFY = D*I2 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseY + GScale*YI[k]
					NFZ = D*I3 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseZ + GScale*ZI[k]
					Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
					Bond.append([BT3,Inter,FlankF,' #C1SFF'])
					Bond.append([BT4,Link,FlankF,' #C1SLFF'])			
					Angl.append([AT3,IC[-2],Inter,FlankF,' #C1SFF'])  
		
			Bond.append([BT3,IC[0],IC[1],' #C1SIB1'])		
		
			for n in range(0,NinterS-2):
				Angl.append([AT3,IC[n],IC[n+1],IC[n+2],' #C1SI'])
				Bond.append([BT3,IC[n+1],IC[n+2],' #C1SIB2'])		
		
				
		for i in C1L:
			BaseA = Atom[i][0]
			BaseX = Atom[i][3] - GScale*XI[k]
			BaseY = Atom[i][4] - GScale*YI[k]
			BaseZ = Atom[i][5] - GScale*ZI[k]
			
			BL = m.sqrt(BaseX**2+BaseY**2+BaseZ**2)
			B = [BaseX/BL,BaseY/BL,BaseZ/BL]
			
			
			# Compute perpendicular vectors for interactive and flanking particle
			
			if abs(BaseX) > 1e-3:
				I2 = r.random()
				I3 = r.random()
				I1 = -(I2*BaseY+I3*BaseZ)/(BaseX)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			elif abs(BaseZ) > 1e-3:
				I1 = r.random()
				I2 = r.random()
				I3 = -(I1*BaseX+I2*BaseY)/(BaseZ)
				
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL  
			
			else:
				I1 = r.random()
				I3 = r.random()
				I2 = -(I1*BaseX+I3*BaseZ)/(BaseY)
			
				IL = m.sqrt(I1**2 + I2**2 + I3**2)
				I1 = I1/IL
				I2 = I2/IL
				I3 = I3/IL
				
			I = [I1,I2,I3]
		
			F1 = np.cross(I,B)
			F2 = np.cross(B,I)    
			
			for j in range(1,NspaceL+1):
				GAtomC = GAtomC + 1
				NewX = (j*Pscale+1)*BaseX  + GScale*XI[k]
				NewY = (j*Pscale+1)*BaseY  + GScale*YI[k]
				NewZ = (j*Pscale+1)*BaseZ  + GScale*ZI[k]
				if j == 1:
					Bond.append([BT1,BaseA,GAtomC,' #C1Lspaj0'])
				else:
					Bond.append([BT1,GAtomC-1,GAtomC,' #C1Lspaj'])
					if j != 2:
						Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C1Lspa'])				
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
			IC = []
			FL1 = []
			FL2 = []     
			for j in range(1,NinterL+1):
				GAtomC = GAtomC + 1
				Link = GAtomC		
				NewX = (NspaceL*Pscale+j*Pscale+1)*BaseX + GScale*XI[k]
				NewY = (NspaceL*Pscale+j*Pscale+1)*BaseY + GScale*YI[k]
				NewZ = (NspaceL*Pscale+j*Pscale+1)*BaseZ + GScale*ZI[k]
				Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
				if j == 1:
					Bond.append([BT2,GAtomC-1,GAtomC,' #C1Lj0'])
					Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C1Lj0'])
				elif j == 2 :		
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1Lj1'])
					Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C1Lj1'])
				else:
					Bond.append([BT2,GAtomC-4,GAtomC,' #C1LSj'])
					Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C1Lj'])		
		
				GAtomC = GAtomC + 1	
				Inter = GAtomC			
				IC.append(GAtomC)
				NIX = D*I1 + NewX
				NIY = D*I2 + NewY
				NIZ = D*I3 + NewZ
				Atom.append([GAtomC,MolID,P1LTB[j-1],NIX,NIY,NIZ])
				Bond.append([BT3,Link,Inter,' #C1LI'])
				
				GAtomC = GAtomC + 1    
				Flank1 = GAtomC
				FL1.append(Flank1)        
				NF1X = D*F1[0] + NewX + D*I1
				NF1Y = D*F1[1] + NewY + D*I2
				NF1Z = D*F1[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
				Bond.append([BT3,Inter,Flank1,' #C1LF1'])      
				Bond.append([BT4,Link,Flank1,'#C1LLF1'])  
				
				GAtomC = GAtomC + 1
				Flank2 = GAtomC
				FL2.append(Flank2) 		
				NF2X = D*F2[0] + NewX + D*I1
				NF2Y = D*F2[1] + NewY + D*I2
				NF2Z = D*F2[2] + NewZ + D*I3
				Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
				Bond.append([BT3,Inter,Flank2,' #C1LF2'])
				Bond.append([BT4,Link,Flank2,' #C1LLF2'])		
				Angl.append([AT3,Flank1,Inter,Flank2,' #C1LF'])
				
				if j == NinterL:
					GAtomC = GAtomC + 1
					FlankF = GAtomC			
					NFX = D*I1 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseX + GScale*XI[k]
					NFY = D*I2 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseY + GScale*YI[k]
					NFZ = D*I3 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseZ + GScale*ZI[k]
					Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
					Bond.append([BT3,Inter,FlankF,' #C1LFF'])
					Bond.append([BT4,Link,FlankF,' #C1LLFF'])			
					Angl.append([AT3,IC[-2],Inter,FlankF,' #C1LFF'])  
		
			Bond.append([BT3,IC[0],IC[1],' #C1LIB1'])
		
			for n in range(0,NinterL-2):
				Angl.append([AT3,IC[n],IC[n+1],IC[n+2],' #C1LI'])
				Bond.append([BT3,IC[n+1],IC[n+2],' #C1LIB2'])		        


# Fix image flags


AtomO = []


for i in Atom:
	if i[3] > GScale*BOXXH:
		XO = i[3] - GScale*(BOXXH-BOXXL)
		NX = 1
	elif i[3] < GScale*BOXXL:	
		XO = i[3] + GScale*(BOXXH-BOXXL)
		NX = -1
	else:
		XO = i[3]
		NX = 0

	if i[4] > GScale*BOXYH:
		YO = i[4] - GScale*(BOXYH-BOXYL)
		NY = 1
	elif i[4] < GScale*BOXYL:	
		YO = i[4] + GScale*(BOXYH-BOXYL)
		NY = -1
	else:
		YO = i[4]
		NY = 0
		
	if i[5] > GScale*BOXZH:
		ZO = i[5] - GScale*(BOXZH-BOXZL)
		NZ = 1
	elif i[5] < GScale*BOXZL:	
		ZO = i[5] + GScale*(BOXZH-BOXZL)
		NZ = -1
	else:
		ZO = i[5]
		NZ = 0		

	AtomO.append([i[0],i[1],i[2],XO,YO,ZO,NX,NY,NZ])


# Output to file
NAtom = len(AtomO)
NBond = len(Bond)
NAngl = len(Angl)

F = open(OUTNAME,'w')
F.write('LAMMPS Description \n')
F.write('\n')
F.write(str(NAtom) + ' atoms \n')
F.write(str(NBond) + ' bonds \n')
F.write(str(NAngl) + ' angles \n')
#F.write('0 dihedrals \n')
#F.write('0 impropers \n')
F.write('\n')
F.write(str(AtTN) + ' atom types \n')
F.write(str(BTN) + ' bond types \n')
F.write(str(ATN) + ' angle types \n')
#F.write('0 dihedral types \n')
#F.write('0 improper types \n')
F.write('\n')
F.write(str(GScale*BOXXL) + ' ' + str(GScale*BOXXH) + ' xlo xhi \n')
F.write(str(GScale*BOXYL) + ' ' + str(GScale*BOXYH) + ' ylo yhi \n')
F.write(str(GScale*BOXZL) + ' ' + str(GScale*BOXZH) + ' zlo zhi \n')
F.write('\n')
F.write('Masses \n')
F.write('\n')
F.write('1 1.0 \n')
F.write('2 1.0 \n')
F.write('3 1.0 \n')
F.write('4 1.0 \n')
F.write('5 1.0 \n')
F.write('6 1.0 \n')
F.write('7 1.0 \n')
F.write('8 1.0 \n')
F.write('9 1.0 \n')
F.write('10 1.0 \n')
F.write('\n')
F.write('Atoms \n')
F.write('\n')
for i in AtomO:
    F.write(str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' 0 ' + str(i[3])+ ' ' + str(i[4])+ ' ' + str(i[5])+ ' ' + str(i[6])+ ' ' + str(i[7]) + ' ' + str(i[8])+ '\n')
F.write('\n')
F.write('Bonds \n')
F.write('\n')
BCount = 0
for i in Bond:
    BCount = BCount + 1
    #F.write(str(BCount) + ' '  + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + i[3] +'\n') # Include Labels
    F.write(str(BCount) + ' '  + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) +'\n')
F.write('\n')
F.write('Angles \n')
F.write('\n')
ACount = 0
for i in Angl:
    ACount = ACount + 1
    #F.write(str(ACount) + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) +  ' ' + str(i[3]) + i[4] + '\n') # Include Labels
    F.write(str(ACount) + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) +  ' ' + str(i[3])+ '\n')

	

	
F.close()	



