import numpy as np
import astropy.units as u
from ReadFile import Read


def ComponentMass(filename, p_type):
    '''
    This function calculates the total mass of a certain particle type in a galaxy .txt file
    Inputs:
        filename (.txt file): name of the file referring to a galaxy simulation
        p_type (int raging from 1 to 3): integer that refers to a particle type
    Outputs:
        TotMass (float): float containing the total mass for a given particle 
            type in units of 10^12 solar masses
    '''
    
    #Using the Read function from the previous assignment and storing data onto the 'data' variable
    data = Read(filename)[2]
    
    #Indexing the data according to the particle type
    index = np.where(data['type'] == p_type)
    
    #Creating a variable for the total mass to be stored
    TotMass = 0
    
    #For loop iterating over every row and adding up the values
    for i in data['m'][index]:
        TotMass += i*10**10*u.Msun
    
    return np.around(TotMass/(10**12*u.Msun), 3)