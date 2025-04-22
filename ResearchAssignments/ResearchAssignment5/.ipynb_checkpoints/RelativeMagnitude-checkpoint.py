#Functions used for ReasearchAssignment 5
# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


def RelativeMagnitude(a, b): #From HW6
    """
    Function to calculate an array of the relative magnitudes between a pair of 3D arrays.
    Example: center of mass position (x,y,z) of galaxy or velocity (vx, vy, vz)
    
    INPUTS:
        array1: 3D numpy array
        array2: 3D numpy array
        
    OUTPUTS:
        relmag: 1D numpy array
            Array of the relative vector magnitudes (array1 - array2)
    """
    A = a[:,0] - b[:,0] # compute relative "x" values
    B = a[:,1] - b[:,1] # compute relative "y" values
    C = a[:,2] - b[:,2] # compute relative "z" values
    
    relmag = np.sqrt(A**2 + B**2 + C**2) # compute magnitudes
    
    return relmag