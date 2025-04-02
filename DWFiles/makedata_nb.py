import math as m
import random as r
import sys as s
import numpy as np



NpolyS = int(s.argv[1]) # Number of long polymer chains in each particle
NpolyL = int(s.argv[2]) # Number of short polymer chains in each particle
NspaceL = int(s.argv[3]) # Number of spacer beads in each long polymer
NspaceS = int(s.argv[4]) # Number of spacer beads in each  short polymer


DX = 50

#Variables

phi = m.pi*(3-m.sqrt(5))

Nbig = 300  # Number of spheres in each core
#NpolyS = 20 # Number of long polymer chains in each particle
#NpolyL = 30 # Number of short polymer chains in each particle
#NspaceS = 1 # Number of spacer beads in each  short polymer
#NspaceL = 14  # Number of spacer beads in each long polymer
NinterS = 4 # Number of interactive beads in each short polymer 
NinterL = 4 # Number of interactive beads in each long polymer

D = 0.8 # Distance between linker and interactive/flanking bead
CRad = 5 # Radius of Core Particle
# Set types

AtTN = 8 # Number of Atom Types

Cor1 = 1 # Core Particle for Particle 1
SpaT = 2 # Spacer Particle
P1ST = [3, 3, 4, 4] # Particle 1 Short Polymer Interactive, Length equal NinterS
P1LT = [5, 5, 6, 6] # Particle 1 Long Polymer Interactive, Length equal NinterL
P2ST = [3, 3, 4, 4] # Particle 2 Short Polymer Interactive, Length equal NinterS
P2LT = [5, 5, 6, 6] # Particle 2 Long Polymer Interactive, Length equal NinterL
FlaT = 7 # Flanking Particle 
Cor2 = 1 # Core Particle for Particle 2

BTN = 4 # Number of Bond Types
BT1 = 1 # Spacer Bonds
BT2 = 2 # Base Bonds
BT3 = 3 # Flanking Bonds
BT4 = 4 # Flanking-Linker Bonds

ATN = 3 # Number of Angle Types
AT1 = 1 # Spacer Angles
AT2 = 2 # Base Angles
AT3 = 3 # Flanking Angles


BoxL = 100  # Box long length side, will be from -BoxL/2 to BoxL/2
BoxH = 60   # Box short length side, will be from -BoxH/2 to BoxH/2
omega = 1   # Size of beads

Cscale = 1 # Scaling on core
Pscale = 0.125 # Scaling on polymer



Mass = 1.0   # Mass of beads

Atom = []
Bond = []
Angl = []
#OUTNAME = TYPE + str(DX) + '.data'
OUTNAME = 'PNP.data'

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
print('Core diameter is ' + str(Dmax))



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
C1S = []
C1L = []
C2S = []
C2L = []


GAtomC = 0
for i in Core:
    GAtomC = GAtomC + 1
    Atom.append([GAtomC,1,Cor1,i[0]-DX/2,i[1],i[2]]) 

for i in CoreL:
    GAtomC = GAtomC + 1
    C1L.append(GAtomC-1)
    Atom.append([GAtomC,1,Cor1,i[0]-DX/2,i[1],i[2]]) 

for i in CoreS:
    GAtomC = GAtomC + 1
    C1S.append(GAtomC-1)
    Atom.append([GAtomC,1,Cor1,i[0]-DX/2,i[1],i[2]]) 

for i in Core:
    GAtomC = GAtomC + 1
    Atom.append([GAtomC,2,Cor2,i[0]+DX/2,i[1],i[2]])

for i in CoreL:
    GAtomC = GAtomC + 1
    C2L.append(GAtomC-1)		
    Atom.append([GAtomC,2,Cor2,i[0]+DX/2,i[1],i[2]])
	
for i in CoreS:
    GAtomC = GAtomC + 1
    C2S.append(GAtomC-1)		
    Atom.append([GAtomC,2,Cor2,i[0]+DX/2,i[1],i[2]])	

### Make polymer chains ###



# Make side chain atoms and bonds

# Choose Atom 
for i in C1S:
    BaseA = Atom[i][0]
    BaseX = Atom[i][3] + DX/2
    BaseY = Atom[i][4]
    BaseZ = Atom[i][5]
    
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
 
    MolID = 3
 
    for j in range(1,NspaceS+1):
        GAtomC = GAtomC + 1
        NewX = (j*Pscale+1)*BaseX - DX/2
        NewY = (j*Pscale+1)*BaseY
        NewZ = (j*Pscale+1)*BaseZ
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
        NewX = (NspaceS*Pscale+j*Pscale+1)*BaseX  - DX/2
        NewY = (NspaceS*Pscale+j*Pscale+1)*BaseY
        NewZ = (NspaceS*Pscale+j*Pscale+1)*BaseZ
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
        Atom.append([GAtomC,MolID,P1ST[j-1],NIX,NIY,NIZ])
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
            NFX = D*I1 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseX - DX/2
            NFY = D*I2 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseY
            NFZ = D*I3 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseZ
            Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
            Bond.append([BT3,Inter,FlankF,' #C1SFF'])
            Bond.append([BT4,Link,FlankF,' #C1SLFF'])			
            Angl.append([AT3,IC[-2],Inter,FlankF,' #C1SFF'])  

    Bond.append([BT3,IC[0],IC[1],' #C1SIB1'])		
#    Bond.append([BT3,FL1[0],FL1[1],' #C1SIBF1'])
#    Bond.append([BT3,FL2[0],FL2[1],' #C1SIBF2'])    
    for k in range(0,NinterS-2):
        Angl.append([AT3,IC[k],IC[k+1],IC[k+2],' #C1SI'])
        Bond.append([BT3,IC[k+1],IC[k+2],' #C1SIB2'])		
#        Bond.append([BT3,FL1[k+1],FL1[k+2],' #C1SIBFF1'])
#        Bond.append([BT3,FL2[k+1],FL2[k+2],' #C1SIBFF2'])         
        
for i in C1L:
    BaseA = Atom[i][0]
    BaseX = Atom[i][3] + DX/2
    BaseY = Atom[i][4]
    BaseZ = Atom[i][5]
    
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

    MolID = 4
    
    for j in range(1,NspaceL+1):
        GAtomC = GAtomC + 1
        NewX = (j*Pscale+1)*BaseX  - DX/2
        NewY = (j*Pscale+1)*BaseY
        NewZ = (j*Pscale+1)*BaseZ
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
        NewX = (NspaceL*Pscale+j*Pscale+1)*BaseX  - DX/2
        NewY = (NspaceL*Pscale+j*Pscale+1)*BaseY
        NewZ = (NspaceL*Pscale+j*Pscale+1)*BaseZ
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
        Atom.append([GAtomC,MolID,P1LT[j-1],NIX,NIY,NIZ])
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
            NFX = D*I1 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseX - DX/2
            NFY = D*I2 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseY
            NFZ = D*I3 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseZ
            Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
            Bond.append([BT3,Inter,FlankF,' #C1LFF'])
            Bond.append([BT4,Link,FlankF,' #C1LLFF'])			
            Angl.append([AT3,IC[-2],Inter,FlankF,' #C1LFF'])  

    Bond.append([BT3,IC[0],IC[1],' #C1LIB1'])
#    Bond.append([BT3,FL1[0],FL1[1],' #C1LIBF1'])
#    Bond.append([BT3,FL2[0],FL2[1],' #C1LIBF2'])       
    for k in range(0,NinterL-2):
        Angl.append([AT3,IC[k],IC[k+1],IC[k+2],' #C1LI'])
        Bond.append([BT3,IC[k+1],IC[k+2],' #C1LIB2'])		        
#        Bond.append([BT3,FL1[k+1],FL1[k+2],' #C1LIBFF1'])
#        Bond.append([BT3,FL2[k+1],FL2[k+2],' #C1LIBFF2'])  		
for i in C2S:
    BaseA = Atom[i][0]
    BaseX = Atom[i][3]  - DX/2
    BaseY = Atom[i][4]
    BaseZ = Atom[i][5]
    
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
      
    MolID = 5
    
    for j in range(1,NspaceS+1):
        GAtomC = GAtomC + 1
        NewX = (j*Pscale+1)*BaseX  + DX/2
        NewY = (j*Pscale+1)*BaseY
        NewZ = (j*Pscale+1)*BaseZ
        if j == 1:
            Bond.append([BT1,BaseA,GAtomC,' #C2Sspaj0'])
        else:
            Bond.append([BT1,GAtomC-1,GAtomC,' #C2Sspaj'])
            if j != 2:
                Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C2Sspa'])
        Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
    IC = []
    FL1 = []
    FL2 = []       
    for j in range(1,NinterS+1):
        GAtomC = GAtomC + 1
        Link = GAtomC		
        NewX = (NspaceS*Pscale+j*Pscale+1)*BaseX  + DX/2
        NewY = (NspaceS*Pscale+j*Pscale+1)*BaseY
        NewZ = (NspaceS*Pscale+j*Pscale+1)*BaseZ
        Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
        if j == 1:
        	Bond.append([BT2,GAtomC-1,GAtomC,' #C2Sj0'])
        	Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C2Sj0'])
        elif j == 2 :		
        	Bond.append([BT2,GAtomC-4,GAtomC,' #C2Sj1'])
        	Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C2Sj0'])
        else:
        	Bond.append([BT2,GAtomC-4,GAtomC,' #C2Sj'])
        	Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C2Sj0'])			

        GAtomC = GAtomC + 1	
        Inter = GAtomC			
        IC.append(GAtomC)
        NIX = D*I1 + NewX
        NIY = D*I2 + NewY
        NIZ = D*I3 + NewZ
        Atom.append([GAtomC,MolID,P2ST[j-1],NIX,NIY,NIZ])
        Bond.append([BT3,Link,Inter,' #C2SI'])
        
        GAtomC = GAtomC + 1  
        Flank1 = GAtomC
        FL1.append(Flank1)	        
        NF1X = D*F1[0] + NewX + D*I1
        NF1Y = D*F1[1] + NewY + D*I2
        NF1Z = D*F1[2] + NewZ + D*I3
        Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
        Bond.append([BT3,Inter,Flank1,' #C2SF1'])      
        Bond.append([BT4,Link,Flank1,'#C2SLF1'])    
        
        GAtomC = GAtomC + 1
        Flank2 = GAtomC	
        FL2.append(Flank2)        
        NF2X = D*F2[0] + NewX + D*I1
        NF2Y = D*F2[1] + NewY + D*I2
        NF2Z = D*F2[2] + NewZ + D*I3
        Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
        Bond.append([BT3,Inter,Flank2,' #C2SF2'])
        Bond.append([BT4,Link,Flank2,' #C2SLF2'])		
        Angl.append([AT3,Flank1,Inter,Flank2,' #C2SF'])

        if j == NinterS:
            GAtomC = GAtomC + 1
            FlankF = GAtomC		
            NFX = D*I1 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseX + DX/2
            NFY = D*I2 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseY
            NFZ = D*I3 + (NspaceS*Pscale+(NinterS+1)*Pscale+1)*BaseZ
            Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
            Bond.append([BT3,Inter,FlankF,' #C2SFF'])
            Bond.append([BT4,Link,FlankF,' #C2SLFF'])			
            Angl.append([AT3,IC[-2],Inter,FlankF,' #C2SFF'])  

    Bond.append([BT3,IC[0],IC[1],' #C2SIB1'])
#    Bond.append([BT3,FL1[0],FL1[1],' #C2SIBF1'])
#    Bond.append([BT3,FL2[0],FL2[1],' #C2SIBF2'])      
    for k in range(0,NinterS-2):
        Angl.append([AT3,IC[k],IC[k+1],IC[k+2],' #C2SI'])
        Bond.append([BT3,IC[k+1],IC[k+2],' #C2SIB2'])        
#        Bond.append([BT3,FL1[k+1],FL1[k+2],' #C2SIBFF1'])
#        Bond.append([BT3,FL2[k+1],FL2[k+2],' #C2SIBFF2'])
			
for i in C2L:
    BaseA = Atom[i][0]
    BaseX = Atom[i][3]  - DX/2
    BaseY = Atom[i][4]
    BaseZ = Atom[i][5]
    
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
    MolID = 6
	
    for j in range(1,NspaceL+1):
        GAtomC = GAtomC + 1
        NewX = (j*Pscale+1)*BaseX  + DX/2
        NewY = (j*Pscale+1)*BaseY
        NewZ = (j*Pscale+1)*BaseZ
        if j == 1:
            Bond.append([BT1,BaseA,GAtomC,' #C2Lspaj0'])
        else:
            Bond.append([BT1,GAtomC-1,GAtomC,' #C2Lspaj'])
            if j != 2:
                Angl.append([AT1,GAtomC-2,GAtomC-1,GAtomC, ' #C2Lspa'])			
        Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
    IC = []
    FL1 = []
    FL2 = []     
    for j in range(1,NinterL+1):
        GAtomC = GAtomC + 1
        Link = GAtomC			
        NewX = (NspaceL*Pscale+j*Pscale+1)*BaseX  + DX/2
        NewY = (NspaceL*Pscale+j*Pscale+1)*BaseY
        NewZ = (NspaceL*Pscale+j*Pscale+1)*BaseZ
        Atom.append([GAtomC,MolID,SpaT,NewX,NewY,NewZ])
        if j == 1:
        	Bond.append([BT2,GAtomC-1,GAtomC,' #C2Lj0'])
        	Angl.append([AT2,GAtomC-2,GAtomC-1,GAtomC,' #C2Lj0'])
        elif j == 2 :		
        	Bond.append([BT2,GAtomC-4,GAtomC,' #C2Lj1'])
        	Angl.append([AT2,GAtomC-5,GAtomC-4,GAtomC,' #C2Lj0'])
        else:
        	Bond.append([BT2,GAtomC-4,GAtomC,' #C2Lj'])
        	Angl.append([AT2,GAtomC-8,GAtomC-4,GAtomC,' #C2Lj0'])		

        GAtomC = GAtomC + 1
        Inter = GAtomC		
        IC.append(GAtomC)
        NIX = D*I1 + NewX
        NIY = D*I2 + NewY
        NIZ = D*I3 + NewZ
        Atom.append([GAtomC,MolID,P2LT[j-1],NIX,NIY,NIZ])
        Bond.append([BT3,Link,Inter,' #C2LI'])
        
        GAtomC = GAtomC + 1    
        Flank1 = GAtomC
        FL1.append(Flank1)        
        NF1X = D*F1[0] + NewX + D*I1
        NF1Y = D*F1[1] + NewY + D*I2
        NF1Z = D*F1[2] + NewZ + D*I3
        Atom.append([GAtomC,MolID,FlaT,NF1X,NF1Y,NF1Z])        
        Bond.append([BT3,Inter,Flank1,' #C2LF1'])      
        Bond.append([BT4,Link,Flank1,'#C2LLF1'])     
        
        GAtomC = GAtomC + 1
        Flank2 = GAtomC		
        FL2.append(Flank2)		
        NF2X = D*F2[0] + NewX + D*I1
        NF2Y = D*F2[1] + NewY + D*I2
        NF2Z = D*F2[2] + NewZ + D*I3
        Atom.append([GAtomC,MolID,FlaT,NF2X,NF2Y,NF2Z])        
        Bond.append([BT3,Inter,Flank2,' #C2LF2'])
        Bond.append([BT4,Link,Flank2,' #C2LLF2'])		
        Angl.append([AT3,Flank1,Inter,Flank2,' #C2LF'])

        if j == NinterL:
            GAtomC = GAtomC + 1
            FlankF = GAtomC				
            NFX = D*I1 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseX + DX/2
            NFY = D*I2 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseY
            NFZ = D*I3 + (NspaceL*Pscale+(NinterL+1)*Pscale+1)*BaseZ
            Atom.append([GAtomC,MolID,FlaT,NFX,NFY,NFZ])        
            Bond.append([BT3,Inter,FlankF,' #C2LFF'])
            Bond.append([BT4,Link,FlankF,' #C2LLFF'])			
            Angl.append([AT3,IC[-2],Inter,FlankF,' #C2LFF'])          
		
    Bond.append([BT3,IC[0],IC[1],' #C2LIB1'])
#    Bond.append([BT3,FL1[0],FL1[1],' #C2LIBF1'])
#    Bond.append([BT3,FL2[0],FL2[1],' #C2LIBF2'])     
    for k in range(0,NinterL-2):
        Angl.append([AT3,IC[k],IC[k+1],IC[k+2],' #C2LI'])
        Bond.append([BT3,IC[k+1],IC[k+2],' #C2LIB2'])        
#        Bond.append([BT3,FL1[k+1],FL1[k+2],' #C2LIBFF1'])
#        Bond.append([BT3,FL2[k+1],FL2[k+2],' #C2LIBFF2'])         
        

# Output to file
NAtom = len(Atom)
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
F.write(str(-BoxL/2) + ' ' + str(BoxL/2) + ' xlo xhi \n')
F.write(str(-BoxH/2) + ' ' + str(BoxH/2) + ' ylo yhi \n')
F.write(str(-BoxH/2) + ' ' + str(BoxH/2) + ' zlo zhi \n')
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
F.write('\n')
F.write('Atoms \n')
F.write('\n')
for i in Atom:
    F.write(str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' 0 ' + str(i[3])+ ' ' + str(i[4])+ ' ' + str(i[5]) + '\n')
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


### Compute properties:
#CRadius = m.sqrt(Core[0][0]**2 + Core[0][1]**2 + Core[0][2]**2)
#print(CRadius)
#print(SLength)
#print(LLength)
#




