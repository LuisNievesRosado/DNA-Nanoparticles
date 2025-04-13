import sys as s
import matplotlib.pyplot as plt
import numpy as np


#import torch # Getting some weird memory errors on my PC that are fixed with this import
			  # not actually used in code.


#============================================================================#
# Variables
#============================================================================#
diameter = 10.0
number_windows = 39

#============================================================================#
# Function to read output from wham
#============================================================================#

def wham_read(filename, number_windows, diameter):
    distance = np.zeros(number_windows)
    PMF = np.zeros(number_windows)
	
	# Read in PMF from wham file
    file = open(filename,'r')
    line = file.readline()
    for idx in range(number_windows):
        line = file.readline()
        spl = line.split()
        distance[idx] = float(spl[0])
        PMF[idx] = float(spl[1])	
    file.close()
    # Zero PMF at large distances
    PMF = PMF - PMF[-1]
	
	# Scale distance by diameter
    distance = distance/diameter

    # Compute force, -grad PMF, note that here Delta D is always constant
    force = -np.gradient(PMF,distance[1]-distance[0])

    return distance, PMF, force

#============================================================================#
# Function to extrapolate small distances
#============================================================================#

def polynomial_extra(n_degree, n_points, min_distance, extra_points 
                             ,distance, PMF, force):
    
	# Use first n_points to compute a polynomial of n_degree 

    constant = np.polyfit(distance[1:n_points], PMF[1:n_points], n_degree)

    new_dist = np.linspace(min_distance,distance[2],
                             num=extra_points,endpoint=False)

    new_PMF = np.zeros(extra_points)
    new_force = np.zeros(extra_points)
	
    for idx in range(n_degree + 1):
        new_PMF = new_PMF + constant[idx]*new_dist**(n_degree-idx)
        if idx != n_degree:
            new_force = (new_force - 
                (n_degree - idx)*constant[idx]*new_dist**(n_degree-idx-1))
    
    # Add extrapolation to initial data
    extra_distance = np.concatenate((new_dist, distance[2:]))
    extra_PMF = np.concatenate((new_PMF, PMF[2:]))
    extra_force = np.concatenate((new_force, force[2:]))

    return extra_distance, extra_PMF, extra_force

#============================================================================#
# Function to output potential file for LAMMPS
#============================================================================#

def output_potential(name, type, distance, PMF, force):
    file = open(name,'w')
    
    file.write('# DATE: 2021-11-17 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
    file.write('# PMF Potential \n')
    file.write('\n')
    file.write(type + '\n')
    file.write('N ' + str(len(distance)) + '\n')
    file.write('\n')
    for idx in range(len(distance)):
        file.write(str(idx + 1) + ' ' + str(distance[idx]) + ' ' 
                   + str(PMF[idx]) + ' ' + str(force[idx]) + '\n')	
    file.close() 	

#============================================================================#
# Process AA, AB and BB data
#============================================================================#

# Create figure of PMF data, add plots as data is processes
plt.figure()

for type in ['AA','AB','BB']:
    filename = type + '.out'
    
    # Read data from wham results    
    distance, PMF, force = wham_read(filename, number_windows, diameter)

    # Add polynomial extrapolation to results 
	# n_degree = 2 best from experimentation
    extra_distance, extra_PMF, extra_force = \
            polynomial_extra(2, 10, 0.00001, 10, distance, PMF, force) 
                                     
    

    # Output potential to file, both original and extrapolated	
    output_potential(type + '.pot', type, distance, PMF, force)
    output_potential(type + 'P.pot', type, extra_distance, 
                     extra_PMF, extra_force)

    # Add plot to figure
    if type == 'AA':
        plot_color = 'r'	
    elif type == 'AB':
        plot_color = 'g'		
    elif type == 'BB':
        plot_color = 'b'		
		
    plt.plot(distance, PMF, label = type, linewidth = 3,
            color=plot_color, marker = 'o')



# Finish plot and output
plt.axhline(y=0.0, linestyle='--', color='k')
plt.axvline(x=1.0, linestyle='--', color='k')
plt.ylim([-5,10])
plt.legend()
plt.ylabel('Potential of Mean Force')
plt.xlabel('R/Diameter')
plt.savefig('PMF.png')