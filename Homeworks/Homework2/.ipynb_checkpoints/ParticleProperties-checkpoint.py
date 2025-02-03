import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, particle_type, particle_n):
    '''
    This function relays important information on the particle that is called by the Inputs
    Inputs: 
        filename: A string containing name of the file that contains the data you to be accessed
        particle_type: An int referring to the type of particle you want information on. This variable can assume
            three distinct values: 1, 2 and 3. 
                Type 1 = Dark Matter particle, 
                Type 2 = Disk particle,
                Type 3 = Halo particle
        particle_n: An int referring to the number of the specific particle you want information about

    Outputs:
        This function outputs three values
        1. D_distance: This is the 3D distance of the particle in regards to the center of the galaxy 
            (astropy unit of kpc)
        2. D_velocity: This is the 3D velocity of the particle in a cartesian coordinate system 
            centered on the center of the galaxy (astropy unit of km/s)
        3. mass: Total mass of the particle (astropy unit of solar masses) 
    '''
    #Calling the Read function and storing the third Output, data, in the variable "data"
    data = Read(filename)[2]
    #Creating an index that associates the column 'type' in data to the Input 'particle_type'
    #This makes it so 'particle_type' can be used to create three categories for the particles
    index = np.where(data['type'] == particle_type)
    
    #Creating a 'mass' variable that stores the value found in the column 'm' for a given 'index' and 'particle_n'
    mass = data['m'][index][particle_n]*10**10*u.Msun

    #Creating a 'position' variable for each direction that stores the value found in the column 'x','y' or 'z' for a given 
    #'index' and 'particle_n' in units of kpc, These values are related to the distance of the particle to the center of the
    #galaxy
    x_position = data['x'][index][particle_n-1]*u.kpc
    y_position = data['y'][index][particle_n-1]*u.kpc
    z_position = data['z'][index][particle_n-1]*u.kpc

    #Creating a 'velocity' variable for each direction that stores the value found in the column 'vx','vy' or 'vz' for a given 
    #'index' and 'particle_n' in units of km/s, These values are related to the velocity of the particle in different directions
    x_velocity = data['vx'][index][particle_n-1]*(u.km/u.s)
    y_velocity = data['vy'][index][particle_n-1]*(u.km/u.s)
    z_velocity = data['vz'][index][particle_n-1]*(u.km/u.s)

    #Calculating the distance of the particle to the center of the galaxy in 3D and storing the value in the variable 'D_distance'
    D_distance = np.around(np.sqrt(x_position**2 + y_position**2 + z_position**2), 3)

    #Calculating the velocity of the particle in 3D and storing the value in the variable 'D_velocity'
    D_velocity = np.around(np.sqrt(x_velocity**2 + y_velocity**2 + z_velocity**2), 3)

    #Returning 'D_distance', 'D_velocity' and 'mass'
    return D_distance, D_velocity, mass