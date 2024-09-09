import sys as s
import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.exists('data'):
	os.makedirs('data')


K = 3

Nmax = 39
Dmin = 12
Dmax = 50.0




for n in ['AA','AB','BB']:
    R0DAT = []
    RDAT = []
    TIME = []
    for i in range(0,Nmax):
        DATA = []
        Infile = n + '/log.lammps.' + str(i)
        F = open(Infile,'r')

        InRun = 0 # Flag to set if the line is in a run or not.
    
        for line in F:
            spl = line.split()
            if InRun == 2:   # Line is in the run data
                if spl[0] == 'Loop': # Run is over, organize data
                    InRun = 0
                    NumData = []      # Will organize numerical data here
                    NumL = len(Label) # Number of columns of data
                    DC = {}           # Dictionary that will be added to DATA
                    for i in range(0,NumL):
                        NumData.append([])
                        for j in Num:
                            NumData[i].append(float(j[i]))
                        DC[Label[i]] = NumData[i]
                    DATA.append(DC)   # Save dictionary to DATA
                            
                else:
                    Num.append(spl)     
        
            if InRun == 1:   # Line is in the run labels
                InRun = 2 
                for i in spl:
                    Label.append(i)
        
            if spl != []:
                if spl[0] == 'Per': # Simulation data is about to start
                    InRun = 1
                    Label = []      # Clear labels previous run/ 
                    Num = []        # Clear data previous run     
        R0DAT.append(DATA[0]['v_for_r0'])
        RDAT.append(DATA[0]['c_for_temper'])
        TIME.append(DATA[0]['Step'])
        F.close()
    # Reorganize data into histograms    
    
    R0DATord = []
    RDAThist = []
    TIMEhist = []
    
    for i in np.linspace(Dmin,Dmax,Nmax):
        RDAThist.append([])
        TIMEhist.append([])
        iR = round(i,1)
        R0DATord.append(iR)
        
        for j in range(0,Nmax):
            for k in range(0,len(R0DAT[j])):
                if R0DAT[j][k] == iR:
                    RDAThist[-1].append(RDAT[j][k])
                    TIMEhist[-1].append(TIME[j][k])
    
    # Sort data by time
    
    RDATsort = []
    TIMEsort = []

    for i in range(0,len(TIMEhist)):
        RDATsort.append([])
        TIMEsort.append([])
        xy = zip(TIMEhist[i],RDAThist[i])
        xysort = sorted(xy)
        for j in xysort:
            TIMEsort[-1].append(j[0])
            RDATsort[-1].append(j[1])
    
    # Output files
    
    Mout = open('Meta' + n + '.wham','w')
    for i in range(0,len(R0DATord)):
        Mout.write('data/' + n + 'com_' + str(R0DATord[i]) + '.wham' + ' ' + str(R0DATord[i]) +  ' ' + str(K) + '\n' )
        Wout = open('data/' + n + 'com_' + str(R0DATord[i]) + '.wham','w')
        for j in range(0,len(RDATsort[i])):
            Wout.write(str(TIMEsort[i][j]) + ' ' + str(RDATsort[i][j]) + '\n')
        Wout.close()
    Mout.close()
    
    # Plotting
    
    plt.figure()
    for i in range(0,len(R0DATord)):
        plt.hist(RDATsort[i],20)
    plt.title(n + ' Histograms')    
    plt.savefig(n + 'Histo.png')
    #plt.show()
    


                


            
    
