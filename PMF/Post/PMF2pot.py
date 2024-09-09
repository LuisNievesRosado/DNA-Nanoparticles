import sys as s
import matplotlib.pyplot as plt
import numpy as np


INFILE1 = 'AA.out' 
INFILE2 = 'AB.out' 
INFILE3 = 'BB.out' 
Dia = 10.0
N = 39

F1 = open(INFILE1,'r')
F1.readline()

R1 = []
PMF1 = []

for i in range(0,N):
    A = F1.readline()
    spl = A.split()
    R1.append(float(spl[0]))
    PMF1.append(float(spl[1]))
    
F1.close()    

F2 = open(INFILE2,'r')
F2.readline()

R2 = []
PMF2 = []

for i in range(0,N):
    A = F2.readline()
    spl = A.split()
    R2.append(float(spl[0]))
    PMF2.append(float(spl[1]))
    
F2.close() 
   
F3 = open(INFILE3,'r')
F3.readline()

R3 = []
PMF3 = []

for i in range(0,N):
    A = F3.readline()
    spl = A.split()
    R3.append(float(spl[0]))
    PMF3.append(float(spl[1]))
F3.close() 


for i in range(0,N):
    PMF1[i] = PMF1[i] - PMF1[-1]
    PMF2[i] = PMF2[i] - PMF2[-1]
    PMF3[i] = PMF3[i] - PMF3[-1]


for i in range(0,N):
    R1[i] = R1[i]/Dia
    R2[i] = R2[i]/Dia
    R3[i] = R3[i]/Dia

DR1 = R1[1]-R1[0]
DR2 = R2[1]-R2[0]
DR3 = R3[1]-R3[0]

For1 = np.gradient(PMF1,DR1)
For2 = np.gradient(PMF2,DR2)
For3 = np.gradient(PMF3,DR3)

for i in range(0,N):
    For1[i] = -For1[i] 
    For2[i] = -For2[i] 
    For3[i] = -For3[i]

# Polymer extra for n = 5
#    
#PolyC1 = np.polyfit(R1[0:10],PMF1[0:10],5)
#PolyC2 = np.polyfit(R2[0:10],PMF2[0:10],5)
#PolyC3 = np.polyfit(R3[0:10],PMF3[0:10],5)
#
#Poly1 = [] 
#Poly2 = [] 
#Poly3 = [] 
#
#PolyF1 = [] 
#PolyF2 = [] 
#PolyF3 = [] 
#
#for i in range(0,10):
#    Poly1.append(PolyC1[0]*R1[i]**5+PolyC1[1]*R1[i]**4+PolyC1[2]*R1[i]**3+PolyC1[3]*R1[i]**2+PolyC1[4]*R1[i]+PolyC1[5])
#    Poly2.append(PolyC2[0]*R2[i]**5+PolyC2[1]*R2[i]**4+PolyC2[2]*R2[i]**3+PolyC2[3]*R2[i]**2+PolyC2[4]*R2[i]+PolyC2[5])
#    Poly3.append(PolyC3[0]*R3[i]**5+PolyC3[1]*R3[i]**4+PolyC3[2]*R3[i]**3+PolyC3[3]*R3[i]**2+PolyC3[4]*R3[i]+PolyC3[5])
#    PolyF1.append(-5*PolyC1[0]*R1[i]**4-4*PolyC1[1]*R1[i]**3-3*PolyC1[2]*R1[i]**2-2*PolyC1[3]*R1[i]-PolyC1[4])
#    PolyF2.append(-5*PolyC2[0]*R2[i]**4-4*PolyC2[1]*R2[i]**3-3*PolyC2[2]*R2[i]**2-2*PolyC2[3]*R2[i]-PolyC2[4])
#    PolyF3.append(-5*PolyC3[0]*R3[i]**4-4*PolyC3[1]*R3[i]**3-3*PolyC3[2]*R3[i]**2-2*PolyC3[3]*R3[i]-PolyC3[4])
#
#
## Add polymer extrapolation for small r:
#Rmin = 0.01
#Nsmall =10
#Rs1 = np.linspace(Rmin,R1[0],num=Nsmall,endpoint=False)
#Rs2 = np.linspace(Rmin,R2[0],num=Nsmall,endpoint=False)
#Rs3 = np.linspace(Rmin,R3[0],num=Nsmall,endpoint=False)
#
#RF1 = []
#RF2 = []
#RF3 = []
#
#FPMF1 = []
#FPMF2 = []
#FPMF3 = []
#
#FF1 = []
#FF2 = []
#FF3 = []
#
#for i in range(0,Nsmall):
#	RF1.append(Rs1[i])
#	RF2.append(Rs2[i])
#	RF3.append(Rs3[i])
#	FPMF1.append(PolyC1[0]*Rs1[i]**5+PolyC1[1]*Rs1[i]**4+PolyC1[2]*Rs1[i]**3+PolyC1[3]*Rs1[i]**2+PolyC1[4]*Rs1[i]+PolyC1[5])
#	FPMF2.append(PolyC2[0]*Rs2[i]**5+PolyC2[1]*Rs2[i]**4+PolyC2[2]*Rs2[i]**3+PolyC2[3]*Rs2[i]**2+PolyC2[4]*Rs2[i]+PolyC2[5])
#	FPMF3.append(PolyC3[0]*Rs3[i]**5+PolyC3[1]*Rs3[i]**4+PolyC3[2]*Rs3[i]**3+PolyC3[3]*Rs3[i]**2+PolyC3[4]*Rs3[i]+PolyC3[5])
#	FF1.append(-5*PolyC1[0]*Rs1[i]**4-4*PolyC1[1]*Rs1[i]**3-3*PolyC1[2]*Rs1[i]**2-2*PolyC1[3]*Rs1[i]-PolyC1[4]) 
#	FF2.append(-5*PolyC2[0]*Rs2[i]**4-4*PolyC2[1]*Rs2[i]**3-3*PolyC2[2]*Rs2[i]**2-2*PolyC2[3]*Rs2[i]-PolyC2[4]) 
#	FF3.append(-5*PolyC3[0]*Rs3[i]**4-4*PolyC3[1]*Rs3[i]**3-3*PolyC3[2]*Rs3[i]**2-2*PolyC3[3]*Rs3[i]-PolyC3[4]) 
#
#for i in range(0,N):
#	RF1.append(R1[i])
#	RF2.append(R2[i])
#	RF3.append(R3[i])
#	FPMF1.append(PMF1[i])
#	FPMF2.append(PMF2[i])
#	FPMF3.append(PMF3[i])
#	FF1.append(For1[i])
#	FF2.append(For2[i])
#	FF3.append(For3[i])
#

# Polymer extra for n = 4
    
#PolyC1 = np.polyfit(R1[1:10],PMF1[1:10],4)
#PolyC2 = np.polyfit(R2[1:10],PMF2[1:10],4)
#PolyC3 = np.polyfit(R3[1:10],PMF3[1:10],4)
#
#
#
#Poly1 = [] 
#Poly2 = [] 
#Poly3 = [] 
#
#PolyF1 = [] 
#PolyF2 = [] 
#PolyF3 = [] 
#
#for i in range(0,10):
#    Poly1.append(PolyC1[0]*R1[i]**4+PolyC1[1]*R1[i]**3+PolyC1[2]*R1[i]**2+PolyC1[3]*R1[i]+PolyC1[4])
#    Poly2.append(PolyC2[0]*R2[i]**4+PolyC2[1]*R2[i]**3+PolyC2[2]*R2[i]**2+PolyC2[3]*R2[i]+PolyC2[4])
#    Poly3.append(PolyC3[0]*R3[i]**4+PolyC3[1]*R3[i]**3+PolyC3[2]*R3[i]**2+PolyC3[3]*R3[i]+PolyC3[4])
#    PolyF1.append(-4*PolyC1[0]*R1[i]**4-3*PolyC1[1]*R1[i]**2-2*PolyC1[2]*R1[i]-PolyC1[3])
#    PolyF2.append(-4*PolyC2[0]*R2[i]**4-3*PolyC2[1]*R2[i]**2-2*PolyC2[2]*R2[i]-PolyC2[3])
#    PolyF3.append(-4*PolyC3[0]*R3[i]**4-3*PolyC3[1]*R3[i]**2-2*PolyC3[2]*R3[i]-PolyC3[3])
#
#
## Add polymer extrapolation for small r:
#Rmin = 0.01
#Nsmall =10
#Rs1 = np.linspace(Rmin,R1[0],num=Nsmall,endpoint=False)
#Rs2 = np.linspace(Rmin,R2[0],num=Nsmall,endpoint=False)
#Rs3 = np.linspace(Rmin,R3[0],num=Nsmall,endpoint=False)
#
#RF1 = []
#RF2 = []
#RF3 = []
#
#FPMF1 = []
#FPMF2 = []
#FPMF3 = []
#
#FF1 = []
#FF2 = []
#FF3 = []
#
#for i in range(0,Nsmall):
#	RF1.append(Rs1[i])
#	RF2.append(Rs2[i])
#	RF3.append(Rs3[i])
#	FPMF1.append(PolyC1[0]*Rs1[i]**4+PolyC1[1]*Rs1[i]**3+PolyC1[2]*Rs1[i]**2+PolyC1[3]*Rs1[i]+PolyC1[4])
#	FPMF2.append(PolyC2[0]*Rs2[i]**4+PolyC2[1]*Rs2[i]**3+PolyC2[2]*Rs2[i]**2+PolyC2[3]*Rs2[i]+PolyC2[4])
#	FPMF3.append(PolyC3[0]*Rs3[i]**4+PolyC3[1]*Rs3[i]**3+PolyC3[2]*Rs3[i]**2+PolyC3[3]*Rs3[i]+PolyC3[4])
#	FF1.append(-4*PolyC1[0]*Rs1[i]**3-3*PolyC1[1]*Rs1[i]**2-2*PolyC1[2]*Rs1[i]-PolyC1[3]) 
#	FF2.append(-4*PolyC2[0]*Rs2[i]**3-3*PolyC2[1]*Rs2[i]**2-2*PolyC2[2]*Rs2[i]-PolyC2[3]) 
#	FF3.append(-4*PolyC3[0]*Rs3[i]**3-3*PolyC3[1]*Rs3[i]**2-2*PolyC3[2]*Rs3[i]-PolyC3[3]) 
#
#for i in range(1,N):
#	RF1.append(R1[i])
#	RF2.append(R2[i])
#	RF3.append(R3[i])
#	FPMF1.append(PMF1[i])
#	FPMF2.append(PMF2[i])
#	FPMF3.append(PMF3[i])
#	FF1.append(For1[i])
#	FF2.append(For2[i])
#	FF3.append(For3[i])

# Polymer extra for n = 4
    
PolyC1 = np.polyfit(R1[1:10],PMF1[1:10],2)
PolyC2 = np.polyfit(R2[1:10],PMF2[1:10],2)
PolyC3 = np.polyfit(R3[1:10],PMF3[1:10],2)



Poly1 = [] 
Poly2 = [] 
Poly3 = [] 

PolyF1 = [] 
PolyF2 = [] 
PolyF3 = [] 

for i in range(0,10):
    Poly1.append(PolyC1[0]*R1[i]**2+PolyC1[1]*R1[i]+PolyC1[2])
    Poly2.append(PolyC2[0]*R2[i]**2+PolyC2[1]*R2[i]+PolyC2[2])
    Poly3.append(PolyC3[0]*R3[i]**2+PolyC3[1]*R3[i]+PolyC3[2])
    PolyF1.append(-2*PolyC1[0]*R1[i]-PolyC1[1])
    PolyF2.append(-2*PolyC2[0]*R2[i]-PolyC2[1])
    PolyF3.append(-2*PolyC3[0]*R3[i]-PolyC3[1])


# Add polymer extrapolation for small r:
Rmin = 0.00001
Nsmall = 10
Rs1 = np.linspace(Rmin,R1[2],num=Nsmall,endpoint=False)
Rs2 = np.linspace(Rmin,R2[2],num=Nsmall,endpoint=False)
Rs3 = np.linspace(Rmin,R3[2],num=Nsmall,endpoint=False)

RF1 = []
RF2 = []
RF3 = []

FPMF1 = []
FPMF2 = []
FPMF3 = []

FF1 = []
FF2 = []
FF3 = []

for i in range(0,Nsmall):
	RF1.append(Rs1[i])
	RF2.append(Rs2[i])
	RF3.append(Rs3[i])
	FPMF1.append(PolyC1[0]*Rs1[i]**2+PolyC1[1]*Rs1[i]+PolyC1[2])
	FPMF2.append(PolyC2[0]*Rs2[i]**2+PolyC2[1]*Rs2[i]+PolyC2[2])
	FPMF3.append(PolyC3[0]*Rs3[i]**2+PolyC3[1]*Rs3[i]+PolyC3[2])
	FF1.append(-2*PolyC1[0]*Rs1[i]-PolyC1[1]) 
	FF2.append(-2*PolyC2[0]*Rs2[i]-PolyC2[1]) 
	FF3.append(-2*PolyC3[0]*Rs3[i]-PolyC3[1]) 

for i in range(2,N):
	RF1.append(R1[i])
	RF2.append(R2[i])
	RF3.append(R3[i])
	FPMF1.append(PMF1[i])
	FPMF2.append(PMF2[i])
	FPMF3.append(PMF3[i])
	FF1.append(For1[i])
	FF2.append(For2[i])
	FF3.append(For3[i])













# Find minima


for i in range(0,len(PMF1)):
	if PMF1[i] > PMF1[i+1]:
		Min1 = PMF1[i+1]
	if PMF1[i] < PMF1[i+1]:
		break
for i in range(0,len(PMF2)):
	if PMF2[i] > PMF2[i+1]:
		Min2 = PMF2[i+1]
	if PMF2[i] < PMF2[i+1]:
		break	
		
for i in range(0,len(PMF3)):
	if PMF3[i] > PMF3[i+1]:
		Min3 = PMF3[i+1]
	if PMF3[i] < PMF3[i+1]:
		break		
		
		
	
MinI1 = PMF1.index(Min1)
MinI2 = PMF2.index(Min2)
MinI3 = PMF3.index(Min3)
MinR1 = R1[MinI1]
MinR2 = R2[MinI2]
MinR3 = R3[MinI3]
print('eps/eps = ' + str(Min2/Min1))
print('eps/eps = ' + str(Min2/Min3))
print('sig/sig = ' + str(MinR2/MinR1))
print('sig/sig = ' + str(MinR2/MinR3))





    
plt.figure()
plt.plot(R1,PMF1,label = 'AA',linewidth = 3,color='r',marker = 'o')
plt.plot(R3,PMF3,label = 'BB',linewidth = 3,color='b',marker = 'o')
plt.plot(R2,PMF2,label = 'AB',linewidth = 3,color='g',marker = 'o')
plt.axhline(y=0,linestyle='--',color='k')
plt.axvline(x=10/Dia,linestyle='--',color='k')
#plt.axvline(x=R1[MinI1],color='k')
#plt.axvline(x=R2[MinI2],color='k')
plt.ylim([-5,10])
plt.legend()
plt.ylabel('Potential of Mean Force')
plt.xlabel('R/Diameter')
plt.savefig('PMF.png')
#plt.show()


#plt.figure()
#plt.plot(R1,PMF1,label = 'AA',linewidth = 3,color='r')
#plt.plot(R3,PMF3,label = 'BB',linewidth = 3,color='b')
#plt.plot(R2,PMF2,label = 'AB',linewidth = 3,color='g')
#plt.plot(RF1,FPMF1,label = 'AA',linewidth = 3,color='r',linestyle='--')
#plt.plot(RF3,FPMF3,label = 'BB',linewidth = 3,color='b',linestyle='--')
#plt.plot(RF2,FPMF2,label = 'AB',linewidth = 3,color='g',linestyle='--')
#plt.axhline(y=0,linestyle='--',color='k')
#plt.axvline(x=10/Dia,linestyle='--',color='k')
#plt.legend()
#plt.ylabel('Potential of Mean Force')
#plt.xlabel('R/Diameter')
#plt.show()
#
#
#
#plt.figure()
#plt.plot(R1,For1,label = 'AA',linewidth = 3,color='r')
#plt.plot(R3,For3,label = 'BB',linewidth = 3,color='b')
#plt.plot(R2,For2,label = 'AB',linewidth = 3,color='g')
#plt.axhline(y=0,linestyle='--',color='k')
#plt.axvline(x=10/Dia,linestyle='--',color='k')
##plt.axvline(x=R1[MinI1],color='k')
##plt.axvline(x=R2[MinI2],color='k')
#plt.legend()
#plt.ylabel('Potential of Mean Force')
#plt.xlabel('R/Diameter')
#plt.show()
#
#plt.figure()
#plt.plot(R1,For1,label = 'AA',linewidth = 3,color='r')
#plt.plot(R3,For3,label = 'BB',linewidth = 3,color='b')
#plt.plot(R2,For2,label = 'AB',linewidth = 3,color='g')
#plt.plot(RF1,FF1,label = 'AA',linewidth = 3,color='r',linestyle='--')
#plt.plot(RF3,FF3,label = 'BB',linewidth = 3,color='b',linestyle='--')
#plt.plot(RF2,FF2,label = 'AB',linewidth = 3,color='g',linestyle='--')
#plt.axhline(y=0,linestyle='--',color='k')
#plt.axvline(x=10/Dia,linestyle='--',color='k')
#plt.legend()
#plt.ylabel('Potential of Mean Force')
#plt.xlabel('R/Diameter')
#plt.show()

    
    
G1 = open('AA.pot','w')
G1.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G1.write('# PMF Potential \n')
G1.write('\n')
G1.write('AA\n')
G1.write('N '+ str(N) +'\n')
G1.write('\n')
for i in range(0,N):
	G1.write(str(i+1) + ' '+ str(R1[i]) + ' ' + str(PMF1[i]) + ' ' +str(For1[i]) +'\n')
G1.close()

G2 = open('AB.pot','w')
G2.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G2.write('# PMF Potential \n')
G2.write('\n')
G2.write('AB\n')
G2.write('N '+ str(N) +'\n')
G2.write('\n')
for i in range(0,N):
	G2.write(str(i+1) + ' '+ str(R2[i]) + ' ' + str(FPMF2[i]) + ' ' +str(For2[i]) +'\n')
G2.close()

G3 = open('BB.pot','w')
G3.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G3.write('# PMF Potential \n')
G3.write('\n')
G3.write('BB\n')
G3.write('N '+ str(N) +'\n')
G3.write('\n')
for i in range(0,N):
	G3.write(str(i+1) + ' '+ str(R3[i]) + ' ' + str(PMF3[i]) + ' ' +str(For3[i]) +'\n')
G3.close()



GP1 = open('AAP.pot','w')
GP1.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
GP1.write('# PMF Potential \n')
GP1.write('\n')
GP1.write('AA\n')
GP1.write('N '+ str(len(RF1)) +'\n')
GP1.write('\n')
for i in range(0,len(RF1)):
	GP1.write(str(i+1) + ' '+ str(RF1[i]) + ' ' + str(FPMF1[i]) + ' ' +str(FF1[i]) +'\n')
GP1.close()

GP2 = open('ABP.pot','w')
GP2.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
GP2.write('# PMF Potential \n')
GP2.write('\n')
GP2.write('AB\n')
GP2.write('N '+ str(len(RF2)) +'\n')
GP2.write('\n')
for i in range(0,len(RF2)):
	GP2.write(str(i+1) + ' '+ str(RF2[i]) + ' ' + str(FPMF2[i]) + ' ' +str(FF2[i]) +'\n')
GP2.close()

GP3 = open('BBP.pot','w')
GP3.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
GP3.write('# PMF Potential \n')
GP3.write('\n')
GP3.write('BB\n')
GP3.write('N '+ str(len(RF3)) +'\n')
GP3.write('\n')
for i in range(0,len(RF3)):
	GP3.write(str(i+1) + ' '+ str(RF3[i]) + ' ' + str(FPMF3[i]) + ' ' +str(FF3[i]) +'\n')
GP3.close()