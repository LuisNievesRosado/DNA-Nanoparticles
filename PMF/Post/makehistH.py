import sys as s
import matplotlib.pyplot as plt
import numpy as np
import os

#============================================================================#
# Simulation Variables 
#============================================================================#


K_spring = 3            # Spring constant between two particles


number_windows = 39     # Number of windows
min_dist = 12.0         # Smallest window distance 
max_dist = 50.0         # Larges window distance

#============================================================================#
# Create folder to organize data files
#============================================================================#

if not os.path.exists('data'):
    os.makedirs('data')

#============================================================================#
# Function to read data from a lammps log file
#============================================================================#

def read_lammps(filename):
    lammps_data = []    # Final data, will be a set of dictionaries
						# lammps_data[i] = data of ith+1 run accessible as a 
						# dictionary with keys equal to the column label 
						# and values the data of that column
						
    file = open(filename,'r') # Open file
    
    run_flag = 0	    # Flag to set where line is in output file 
					    # 0 = Out of run; 1 = In run labels; 2 = In run data
    
	# Go through file line by line and checking if we are in a run
    for line in file:
        spl = line.split()
        if run_flag == 2:                           # Line is in the run or at the end of run
            if spl[0] == 'Loop':                    # End of run is over, organize data of this run
                run_flag = 0                        # Reset flag
                final_run_data = []                 # Will organize numerical data here
                number_labels = len(data_label)     # Number of columns of data
                run_dict = {}                       # Dictionary that will have data of this run     
                for idx in range(number_labels):    # Loop over each column
                    final_run_data.append([])       # Separate list per column 
                    for jdx in data_values:         # Loop over all data in run
                        final_run_data[idx].append(float(jdx[idx])) # Add data
                    run_dict[data_label[idx]] = final_run_data[idx] # Populate dictionary
                lammps_data.append(run_dict)        # Save dictionary to output
            else:                                   
                data_values.append(spl)             # Save data of all columns   
        
        if run_flag == 1:                           # Line is in the run labels
            run_flag = 2                            
            for label in spl:                       
                data_label.append(label)            
        
        if spl:                                     # If line is not empty 
            if spl[0] == 'Per':                     # Simulation data is about to start
                run_flag = 1                        
                data_label = []                     # Clear labels previous run 
                data_values = []                    # Clear data previous run     	
    
    file.close()
    return lammps_data

#============================================================================#
# Function to reorganize data, bringing together all the data with same R0
#============================================================================#

def organize_R0(time,data_R0,dist_COM,min_dist,max_dist,number_windows):
    R0_order = []

    COM_out_of_order = []
    time_out_of_order = []

    COM_order = []
    time_order = []
    
	
    for idx in np.linspace(min_dist,max_dist,number_windows): 
        idx_int = round(idx,1) # In case of finite precision issues with linspace
		
		# We will now search all data for points with R0 = idxR and store it.
		# The data will not be ordered in time, we do that in a later step.
        R0_order.append(idx_int)	
		
        COM_out_of_order.append([])
        time_out_of_order.append([])
        
        for jdx in range(number_windows):               # Go through data, sim by sim
            for kdx in range(len(data_R0[jdx])):        # Check R0's of simulation
                if data_R0[jdx][kdx] == R0_order[-1]:   # If we have the right R0 store it 
                    COM_out_of_order[-1].append(dist_COM[jdx][kdx])
                    time_out_of_order[-1].append(time[jdx][kdx])
		
	# We now order the data in time
         
    for idx in range(len(time_out_of_order)):
        COM_order.append([])
        time_order.append([])
        temp_sort = sorted(zip(time_out_of_order[idx],COM_out_of_order[idx]))
        for data_point in temp_sort:
            time_order[-1].append(data_point[0])		
            COM_order[-1].append(data_point[1])

    return R0_order, time_order, COM_order

#============================================================================#
# For all three simulations, loop over output files and read and organize data
#============================================================================#


for sim_type in ['AA','AB','BB']:
    data_R0 = []          # Store spring goal distance
    dist_COM = []         # Store COM real distance 
    time = []             # Track simulation time
    for idx in range(number_windows):
	
		# Read data of logfile idx.
        infile = sim_type + '/log.lammps.' + str(idx)
        data_run = read_lammps(infile) 

		# Save simulation time, R0 and COM distances
        time.append(data_run[0]['Step'])
        data_R0.append(data_run[0]['v_for_r0'])
        dist_COM.append(data_run[0]['c_for_temper'])

    # Organize data in R0
    R0_order, time_order, COM_order = organize_R0(time, data_R0, dist_COM,
                                                  min_dist, max_dist,
                                                  number_windows)
    
	
    # Output files for wham post processing
    
    meta_file = open('Meta' + sim_type + '.wham','w')
    for idx in range(len(R0_order)):
        # Meta file has directions to all individual data files of that simulation 	
        meta_file.write('data/' + sim_type + 'com_' + str(R0_order[idx]) 
                        + '.wham' + ' ' + str(R0_order[idx]) +  ' ' 
                        + str(K_spring) + '\n' )
        
		
		# Data file has distance vs time for all R0's
        data_file = open('data/' + sim_type + 'com_' 
                         + str(R0_order[idx]) + '.wham','w')
						 
        for jdx in range(len(COM_order[idx])):
            data_file.write(str(time_order[idx][jdx]) + ' ' 
                            + str(COM_order[idx][jdx]) + '\n')
							
        data_file.close()
    meta_file.close()
    
    # Output histogram plots
    
    plt.figure()
    for idx in range(len(R0_order)):
        plt.hist(COM_order[idx],20)
    plt.title(sim_type + ' Histograms')    
    plt.savefig(sim_type + 'Histo.png')
    #plt.show()