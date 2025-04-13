import math as m
import random as r
import sys as s
import numpy as np

#============================================================================#
# File inputs
#============================================================================#
pair_type = s.argv[1]                # AA vs AB vs BB

number_chain_short = int(s.argv[2])  # Number of long DNA chains in each particle
number_chain_long = int(s.argv[3])   # Number of short DNA chains in each particle
number_spacer_long = int(s.argv[4])  # Number of spacer beads in each long DNA
number_spacer_short = int(s.argv[5]) # Number of spacer beads in each short DNA

#============================================================================#
# Hard coded variables.
#============================================================================#

# Name of output file
outname = 'PNP.data'     # PNP from Polymer NanoParticle, legacy name before
                         # we settled on DNA 
                         
# Box properties
box_length = 100  # Box long length side, will be from -BoxL/2 to BoxL/2
box_height = 60   # Box short length side, will be from -BoxH/2 to BoxH/2
dx = 50           # Initial Distance between two particles

# Core properties
number_core = 300  # Number of spheres in each core
core_radius = 5    # Radius of Core Particle

# Bead properties
omega = 1.0  # Size of beads
mass = 1.0   # Mass of beads

# Scaling constants
dist_linker = 0.8  # Distance between linker and interactive/flanking bead
site_scale = 1.125 # Scaling on DNA attachment site
DNA_scale = 0.125  # Scaling on DNA chains

# Number of DNA beads
number_DNA_short = 2 # Number of interactive beads in each short DNA
number_DNA_long = 4  # Number of interactive beads in each long DNA

# Type Dictionaries
atom_type = {}
bond_type = {}
angle_type = {}


# Atom Types 
atom_type['number'] = 8 # Number of Atom Types
if pair_type == 'AB':
    atom_type['core_1']  = 1 			# Core Particle for Particle 1
    atom_type['spacer']  = 2 			# Spacer Particle
    atom_type['short_1'] = [5, 3]       # Particle 1 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_1']  = [3, 5, 4, 6] # Particle 1 Long DNA Beads, Length equal number_DNA_long
    atom_type['short_2'] = [3, 5]       # Particle 2 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_2']  = [5, 3, 6, 4] # Particle 2 Long DNA Beads, Length equal number_DNA_long
    atom_type['flanker'] = 7            # Flanking Particle
    atom_type['core_2']  = 8            # Core Particle for Particle 2
elif pair_type == 'AA':
    atom_type['core_1']  = 1 			# Core Particle for Particle 1
    atom_type['spacer']  = 2 			# Spacer Particle
    atom_type['short_1'] = [5, 3]       # Particle 1 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_1']  = [3, 5, 4, 6] # Particle 1 Long DNA Beads, Length equal number_DNA_long
    atom_type['short_2'] = [5, 3]       # Particle 2 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_2']  = [3, 5, 4, 6] # Particle 2 Long DNA Beads, Length equal number_DNA_long
    atom_type['flanker'] = 7            # Flanking Particle
    atom_type['core_2']  = 1            # Core Particle for Particle 2
elif pair_type == 'BB':
    atom_type['core_1']  = 8 			# Core Particle for Particle 1
    atom_type['spacer']  = 2 			# Spacer Particle
    atom_type['short_1'] = [3, 5]       # Particle 1 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_1']  = [5, 3, 6, 4] # Particle 1 Long DNA Beads, Length equal number_DNA_long
    atom_type['short_2'] = [3, 5]       # Particle 2 Short DNA Beads, Length equal number_DNA_short
    atom_type['long_2']  = [5, 3, 6, 4] # Particle 2 Long DNA Beads, Length equal number_DNA_long
    atom_type['flanker'] = 7            # Flanking Particle
    atom_type['core_2']  = 8            # Core Particle for Particle 2

# Bond Types

bond_type['number'] = 4 # Number of Bond Types
bond_type['spacer'] = 1 # Spacer Bonds
bond_type['base']   = 2 # Base Bonds
bond_type['inter']  = 3 # Interactive-spacer Bonds
bond_type['linker'] = 4 # Flanking-Linker Bonds

# Angle Types
angle_type['number'] = 3 # Number of Angle Types
angle_type['spacer'] = 1 # Spacer Angles
angle_type['base']   = 2 # Base Angles
angle_type['flank']  = 3 # Flanking Angles
						 

#============================================================================#
# Create workhorse class
#============================================================================#

class DNA:
    """Class for all atom creation and final output"""
    
    def __init__(self,atom_type, bond_type, angle_type, site_scale, DNA_scale,
                 dist_linker, core_radius, number_core, number_chain_long,
                 number_chain_short):
		
		# Type Dictionaries
        self.atom_type = atom_type
        self.bond_type = bond_type
        self.angle_type = angle_type
		
		# Scaling parameters
        self.site_scale = site_scale		
        self.DNA_scale = DNA_scale
        self.dist_linker = dist_linker
		
		# Template Cores

        self.core_radius = core_radius
        
        self.number_core = number_core
        self.number_chain_long = number_chain_long
        self.number_chain_short= number_chain_short
		
        self.template_core = np.zeros((number_core,3))
        self.template_long = np.zeros((number_chain_long,3))
        self.template_short = np.zeros((number_chain_short,3))

        # Global Atom counter
        self.atom_count = 0
		
		# Atom = [atom_ID, molecule_ID, atom_type, x, y, z] charge = 0 for all
        self.atoms = []
		
        # Bond = [bond_type, atom_ID_1, atom_ID_2]
        self.bonds = []
		
        # Angle = [angle_type, atom_ID_1, atom_ID_2, atom_ID_3] 
        self.angles = []
	
    def make_template_core(self):
	    #"""Create template cores"""
        
        phi = m.pi*(3-m.sqrt(5)) # Constant for particles in a spherical surface
		
		# Make central core, with radius core_radius
        for idx in range(self.number_core):
        	# Creates points in unit sphere roughly uniformly, then scales 
            y = 1 - (idx/(float(self.number_core) - 1))*2
            
            scale_factor =  m.sqrt(1 - y**2)		
            
            theta = phi*idx		
            x = m.cos(theta)*scale_factor		
            z = m.sin(theta)*scale_factor
        	
            self.template_core[idx,:] = [self.core_radius*x,
                                         self.core_radius*y,
                                         self.core_radius*z]

		# Make long core, with radius core_radius*site_scale
        for idx in range(self.number_chain_long):
        	# Creates points in unit sphere roughly uniformly, then scales 
            y = 1 - (idx/(float(self.number_chain_long) - 1))*2
            
            scale_factor =  m.sqrt(1 - y**2)		
            
            theta = phi*idx		
            x = m.cos(theta)*scale_factor		
            z = m.sin(theta)*scale_factor
        	
            self.template_long[idx,:] = [self.core_radius*self.site_scale*x,
                                         self.core_radius*self.site_scale*y,
                                         self.core_radius*self.site_scale*z]    	
		
		# Make long core, with radius core_radius*site_scale, rotate until distance with 
		# the attachments points for long chains is maximized
		
        rotation = np.linspace(-m.pi, m.pi, 50)
        rotation_distance = []
        for angle in rotation:
			# Propose an angle and create points
            temp_template = np.zeros((self.number_chain_short,3))
            for idx in range(self.number_chain_short):
            	# Creates points in unit sphere roughly uniformly, then scales
                x = 1 - (idx/(float(self.number_chain_short) - 1))*2
                
                scale_factor =  m.sqrt(1 - x**2)		
                
                theta = phi*idx + angle		
                y = m.cos(theta)*scale_factor		
                z = m.sin(theta)*scale_factor
            	
                temp_template[idx,:] = [self.core_radius*self.site_scale*x,
                                        self.core_radius*self.site_scale*y,
                                        self.core_radius*self.site_scale*z] 	
		    
			# Check minimum distance between points
            minimum_distance = 10000
            for point_long in self.template_long:
                for point_short in temp_template:
                    distance = np.linalg.norm(point_long-point_short)
                    if distance < minimum_distance:
                        minimum_distance = distance
			
            rotation_distance.append(minimum_distance)			
			
		# Choose angle with maximum rotation distance	
        real_angle = rotation[rotation_distance.index(max(rotation_distance))]	
        print('Attachment sites distance is ' + str(max(rotation_distance)))
		
		# Recreate template	with optimal angle
        for idx in range(self.number_chain_short):
            # Creates points in unit sphere roughly uniformly, then scales 
            x = 1 - (idx/(float(self.number_chain_short) - 1))*2
            
            scale_factor =  m.sqrt(1 - x**2)		
            
            theta = phi*idx + real_angle		
            y = m.cos(theta)*scale_factor		
            z = m.sin(theta)*scale_factor
            
            self.template_short[idx,:] = [self.core_radius*self.site_scale*x,
                                          self.core_radius*self.site_scale*y,
                                          self.core_radius*self.site_scale*z]	
	
    def make_core(self, core_type, molecule_id, x, y, z):
        """Create a core, along with DNA chain attachment sites"""
		
        # Track python index of Core and DNA attachment sites
        core = []
        core_long = []
        core_short = []

		
        # Create core:
        for point in self.template_core:
            core.append(self.atom_count)
            self.atom_count += 1
            self.atoms.append([self.atom_count,molecule_id,core_type,
                               point[0] + x, point[1] + y, point[2] + z])
							
        # Create long core:
        for point in self.template_long:
            core_long.append(self.atom_count)
            self.atom_count += 1
            self.atoms.append([self.atom_count,molecule_id,core_type,
                               point[0] + x, point[1] + y, point[2] + z])
							   
        # Create short core:
        for point in self.template_short:
            core_short.append(self.atom_count)
            self.atom_count += 1
            self.atoms.append([self.atom_count,molecule_id,core_type,
                               point[0] + x, point[1] + y, point[2] + z])							
							
        return core, core_long, core_short
							
    def add_DNA(self, core, number_spacer, particle_type, molecule_id, 
                x, y, z):
        """Add DNA strands to the particle, choose between short and long"""
        interactive_ids = [] # Track id's of all interactive beads, outputted
		
        for site in core:
            site_atom = self.atoms[site][0]  # ID of attachment site
            site_x = self.atoms[site][3] - x # x of attachment site, centered to 0 to compute vectors		
            site_y = self.atoms[site][4] - y # y of attachment site, centered to 0 to compute vectors				
            site_z = self.atoms[site][5] - z # z of attachment site, centered to 0 to compute vectors		 		
						
            site_vec = np.array([site_x, site_y, site_z])
            site_norm = site_vec/np.linalg.norm(site_vec)			
			
			# Compute perpendicular vector to get direction of attached beads
				 
            if abs(site_x) > 1e-3: # Checking that base_x is non_zero
                perp_vec_y = r.random()
                perp_vec_z = r.random()
                perp_vec_x = -(perp_vec_y*site_y + perp_vec_z*site_z)/(site_x)
        
                perp_vec_1 = np.array([perp_vec_x, perp_vec_y, perp_vec_z])
                perp_norm_1 = perp_vec_1/np.linalg.norm(perp_vec_1)				
            
            elif abs(site_z) > 1e-3: # Checking that base_y is non_zero
                perp_vec_x = r.random()
                perp_vec_y = r.random()
                perp_vec_z = -(perp_vec_x*site_x + perp_vec_y*site_y)/(site_z)
        
                perp_vec_1 = np.array([perp_vec_x, perp_vec_y, perp_vec_z])
                perp_norm_1 = perp_vec_1/np.linalg.norm(perp_vec_1)	

            else:                   # They can't all be zero
                perp_vec_x = r.random()
                perp_vec_z = r.random()
                perp_vec_y = -(perp_vec_x*site_x + perp_vec_z*site_z)/(site_y)
        
                perp_vec_1 = np.array([perp_vec_x, perp_vec_y, perp_vec_z])
                perp_norm_1 = perp_vec_1/np.linalg.norm(perp_vec_1)	

			# Remaining perpendicular vector can be found with cross product
            perp_norm_2 = np.cross(perp_norm_1,site_norm)
			
			# Add spacer beads
            for idx in range(1,number_spacer + 1):
                self.atom_count += 1
                
                spacer_id = self.atom_count
                spacer_x = (idx*DNA_scale + 1)*site_x + x
                spacer_y = (idx*DNA_scale + 1)*site_y + y				
                spacer_z = (idx*DNA_scale + 1)*site_z + z				
				
                if idx == 1: # Attachment bond
                    self.bonds.append([self.bond_type['spacer'],
                                       site_atom, spacer_id])
									   
                else: # Spacer - Spacer bond
                    self.bonds.append([self.bond_type['spacer'],
                                       spacer_id - 1, spacer_id])
									   
                    if idx > 2: # No attachment angles
					    # Angle between spacer beads
                        self.angles.append([self.angle_type['spacer'],
                                           spacer_id - 2, spacer_id - 1,
                                           spacer_id])

                self.atoms.append([self.atom_count, molecule_id,
                                   self.atom_type['spacer'],
                                   spacer_x, spacer_y, spacer_z])						   

            # Add base spacers beads, interactive beads and flanking beads
            length_interactive = len(self.atom_type[particle_type])			
            
            inter_chain = []  # Tracking for bonds and angles
            flank_1 = []      # Tracking for bonds and angles 
            flank_2 = []      # Tracking for bonds and angles

            for idx in range(1,length_interactive + 1):
				# Base spacer bead
                self.atom_count += 1
                base_id = self.atom_count # Tracking for bonds and angles

                base_x = (((number_spacer + idx)*self.DNA_scale + 1)
                          *site_x + x)
                base_y = (((number_spacer + idx)*self.DNA_scale + 1)
                          *site_y + y)
                base_z = (((number_spacer + idx)*self.DNA_scale + 1)
                          *site_z + z)

                self.atoms.append([base_id, molecule_id, 
                                   self.atom_type['spacer'],
                                   base_x, base_y, base_z])
								   
                if idx == 1:
					# Bond between base spacer and previous spacer
                    self.bonds.append([self.bond_type['base'],
                                       base_id - 1, base_id])
									   
					# Angle between base spacer and two previous spacers
                    self.angles.append([self.angle_type['base'],
                                        base_id - 2, base_id - 1, base_id])	
										
                elif idx == 2: # Skip the interactive and two flanking
					# Bond between base spacer and previous spacer
                    self.bonds.append([self.bond_type['base'],
                                       base_id - 4, base_id])
									   
     				# Angle between base spacer and two previous spacers
                    self.angles.append([self.angle_type['base'],
                                        base_id - 5, base_id - 4, base_id])
										
                else: # Skip the interactive and two flanking twice
					# Bond between base spacer and previous spacer
                    self.bonds.append([self.bond_type['base'],
                                       base_id - 4, base_id])
									   
     				# Angle between base spacer and two previous spacers
                    self.angles.append([self.angle_type['base'],
                                        base_id - 8, base_id - 4, base_id])

				# Interactive bead - in the direction of perp_norm_1 away from base
                self.atom_count += 1
                inter_id = self.atom_count # Tracking for bonds and angles
                inter_chain.append(inter_id)

                inter_x = base_x + self.dist_linker*perp_norm_1[0] 
                inter_y = base_y + self.dist_linker*perp_norm_1[1] 
                inter_z = base_z + self.dist_linker*perp_norm_1[2] 

                self.atoms.append([inter_id, molecule_id,
                                   self.atom_type[particle_type][idx - 1],
                                   inter_x, inter_y, inter_z])	
								   
                # Bond between interactive and base spacer
                self.bonds.append([self.bond_type['inter'],
                                   base_id, inter_id])
								   
                if idx > 1: 				   
                    # Bond between interactive and previous interactive
                    self.bonds.append([self.bond_type['inter'],
                                       inter_chain[-2], inter_id]) 
									   
                if idx > 2:
					# Angle between interactive and previous interactives
                    self.angles.append([self.angle_type['flank'],
                                inter_chain[-3], inter_chain[-2], inter_id])
				
                # Add flanking bead 1 - perpendicular to spacers and interactive, co-planar with interactive
		
                self.atom_count += 1
                flank_1_id = self.atom_count
                flank_1.append(flank_1_id)
				
                flank_1_x = inter_x + self.dist_linker*perp_norm_2[0]  
                flank_1_y = inter_y + self.dist_linker*perp_norm_2[1]  		
                flank_1_z = inter_z + self.dist_linker*perp_norm_2[2]  		
		
                self.atoms.append([flank_1_id, molecule_id,
                                   self.atom_type['flanker'],
                                   flank_1_x, flank_1_y, flank_1_z])	
								   
                # Bond between flanker and interactive
                self.bonds.append([self.bond_type['inter'],
                                   inter_id, flank_1_id])
								   
                # Bond between flanker and base spacer		
                self.bonds.append([self.bond_type['linker'],
                                   base_id, flank_1_id])		

                # Add flanking bead 2 - perpendicular to spacers and interactive, co-planar with interactive
		
                self.atom_count += 1
                flank_2_id = self.atom_count
                flank_2.append(flank_2_id)
				
                flank_2_x = inter_x - self.dist_linker*perp_norm_2[0]  
                flank_2_y = inter_y - self.dist_linker*perp_norm_2[1]  		
                flank_2_z = inter_z - self.dist_linker*perp_norm_2[2]  		
		
                self.atoms.append([flank_2_id, molecule_id,
                                   self.atom_type['flanker'],
                                   flank_2_x, flank_2_y, flank_2_z])	
								   
                # Bond between flanker and interactive
                self.bonds.append([self.bond_type['inter'],
                                   inter_id, flank_2_id])
								   
                # Bond between flanker and base spacer		
                self.bonds.append([self.bond_type['linker'],
                                   base_id, flank_2_id])
								   
				# Angle between flanking beads and interactive bead				   
                self.angles.append([self.angle_type['flank'],
                                    flank_1_id, inter_id, flank_2_id])


				# Add final flanker to the end, colinear with interactives
                if idx == length_interactive:
                    self.atom_count += 1
                    flank_end_id = self.atom_count
					
                    flank_end_x = self.atoms[inter_chain[-1] - 1][3] \
                                  + self.dist_linker*site_norm[0] 							
									
                    flank_end_y = self.atoms[inter_chain[-1] - 1][4] \
                                  + self.dist_linker*site_norm[1] 								
									
                    flank_end_z = self.atoms[inter_chain[-1] - 1][5] \
                                  + self.dist_linker*site_norm[2]  								
        
                    self.atoms.append([flank_end_id, molecule_id,
                                       self.atom_type['flanker'],
                                       flank_end_x, flank_end_y, flank_end_z])	
					
					# Bond between interactive and final flanker				   
                    self.bonds.append([self.bond_type['inter'],
                                       inter_id, flank_end_id])
									   
					# Bond between interactive and final spacer				   
                    self.bonds.append([self.bond_type['linker'],
                                       inter_id, flank_end_id])		
									   
					# Angle between final flanking beads and previous interactive beads				   
                    self.angles.append([self.angle_type['flank'],
                                    inter_chain[-2], inter_id, flank_end_id])		
    
	
            interactive_ids.append(inter_chain)
	
        return interactive_ids
	
	
    def make_data_file(self,outname,omega,mass,box_length,box_height):
        """Make data file with all atoms, bonds and angles"""
        file = open(outname,'w')
        file.write('LAMMPS Description \n')
        file.write('\n')
        file.write(str(len(self.atoms)) + ' atoms \n')		
        file.write(str(len(self.bonds)) + ' bonds \n')			
        file.write(str(len(self.angles)) + ' angles \n')			
        file.write('\n')			
        file.write(str(self.atom_type['number']) + ' atom types \n')		
        file.write(str(self.bond_type['number']) + ' bond types \n')			
        file.write(str(self.angle_type['number']) + ' angle types \n')			
        file.write('\n')				
        file.write(str(-box_length/2) + ' ' + str(box_length/2) 
                  + ' xlo xhi \n')
        file.write(str(-box_height/2) + ' ' + str(box_height/2)
                  + ' ylo yhi \n')
        file.write(str(-box_height/2) + ' ' + str(box_height/2) 
                  + ' zlo zhi \n')			
        file.write('\n')			
        file.write('Masses \n')				
        file.write('\n')
        for idx in range(self.atom_type['number']):
            file.write(str(idx) + ' ' + str(mass) + '\n')
        file.write('\n')
        file.write('Atoms \n')
        file.write('\n')	
        for atom in self.atoms:
            file.write(str(atom[0]) + ' ' +  str(atom[1]) + ' '
                       + str(atom[2]) + ' 0 ' + str(round(atom[3],4))			
                       + ' ' + str(round(atom[4],4)) + ' ' 
                       + str(round(atom[5],4)) + '\n')			
	
        file.write('\n')
        file.write('Bonds \n')
        file.write('\n')			
        bond_count = 0
        for bond in self.bonds:
            bond_count += 1
            file.write(str(bond_count) + ' ' + str(bond[0]) + ' ' 
                       + str(bond[1]) + ' ' + str(bond[2]) + '\n')		
        file.write('\n')
        file.write('Angles \n')
        file.write('\n')				
        angle_count = 0
        for angle in self.angles:
            angle_count += 1
            file.write(str(angle_count) + ' ' + str(angle[0]) + ' ' 
                       + str(angle[1]) + ' ' + str(angle[2]) + ' '
                       + str(angle[3]) + '\n')				
        file.close()
		
		
#============================================================================#
# Make PMF pair
#============================================================================#
# Call class
PMF = DNA(atom_type,bond_type,angle_type,site_scale,DNA_scale,
          dist_linker,core_radius,number_core,number_chain_long,
          number_chain_short)

# Create core templates
PMF.make_template_core()

# Create particle 1 cores
core_1, core_1_long, core_1_short = PMF.make_core(atom_type['core_1'],
                                                  1,-dx/2, 0.0, 0.0)

# Create particle 1 short DNA chains
interactive_1S = PMF.add_DNA(core_1_short,number_spacer_short,'short_1',3,
                          -dx/2, 0.0, 0.0)

# Create particle 1 long DNA chains
interactive_1L = PMF.add_DNA(core_1_long,number_spacer_long,'long_1',4,
                          -dx/2, 0.0, 0.0)


# Create particle 2 cores
core_2, core_2_long, core_2_short = PMF.make_core(atom_type['core_2'],
                                                  2,dx/2, 0.0, 0.0)

# Create particle 2 short DNA chains
interactive_2S = PMF.add_DNA(core_2_short,number_spacer_short,'short_2',5,
                          dx/2, 0.0, 0.0)

# Create particle 2 long DNA chains
interactive_2L = PMF.add_DNA(core_2_long,number_spacer_long,'long_2',6,
                          dx/2, 0.0, 0.0)

# Make data file

PMF.make_data_file(outname,omega,mass,box_length,box_height)

# Output interactive bead ids for post-processing

file_1S = open('chain1S.txt','w')
for chain in interactive_1S:
    for bead in chain:
        file_1S.write(str(bead) + ' ')
    file_1S.write('\n')
file_1S.close()

file_1L = open('chain1L.txt','w')
for chain in interactive_1L:
    for bead in chain:
        file_1L.write(str(bead) + ' ')
    file_1L.write('\n')
file_1L.close()

file_2S = open('chain2S.txt','w')
for chain in interactive_2S:
    for bead in chain:
        file_2S.write(str(bead) + ' ')
    file_2S.write('\n')
file_2S.close()

file_2L = open('chain2L.txt','w')
for chain in interactive_2L:
    for bead in chain:
        file_2L.write(str(bead) + ' ')
    file_2L.write('\n')
file_2L.close()




















	





















