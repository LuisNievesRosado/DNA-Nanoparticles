import sys as s
import matplotlib.pyplot as plt
import numpy as np
import os

Nmax = 39

for n in ['AA','AB','BB']:
    out_files = [f for f in os.listdir(n) if f.endswith('.out')]
    if len(out_files) != 1:
        raise ValueError('There should be only one .out file in the ' + n + ' directory')
    infile = out_files[0]
    

    TRY = np.zeros(Nmax-1)
    SUC = np.zeros(Nmax-1)
    PROB = np.zeros(Nmax-1)    

    F = open(n + '/' + infile,'r')
    for line in F:
        spl = line.split()
        if spl[0] == 'SWAP':
            ACCEPT = spl[6][0] # Remove trailing comma
            S1 = int(spl[8])
            S2 = int(spl[9].rstrip(spl[9][-1])) # Remove trailing comma
            S = min(S1,S2)
            TRY[S] = TRY[S] + 1
            if ACCEPT == '1':
                SUC[S] = SUC[S] + 1      
    F.close()
    
    for i in range(0,Nmax-1):
        if TRY[i] == 0:
            print('The ' + str(i) + ' ' + str(i+1) + ' swap was never attempted, setting probability negative')
            PROB[i] = -1
        else:
            PROB[i] = SUC[i]/TRY[i]*100
    
    
    plt.figure()
    plt.bar(range(0,Nmax-1),PROB)
    plt.title(n + ' Acceptance Probabilites')
    plt.savefig(n + 'Acp')
    #plt.show()
    
    
    
