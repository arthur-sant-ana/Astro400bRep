import numpy as np
import astropy.units as u
#Defining a reading function that reads a file and stores certain values
def Read(filename):
    '''
    This function reads a file and returns some specific values
    Inputs:
        filename: takes a filename as input, this is the file you want to be read
    Outputs:
        time: Time in astropy Myrs, corresponds to the first line of the document
        N_part: # of particles in the system, second line of the document
        data: An array with the rest of the data in the document
    '''
    file = open(filename, 'r')

    #Reading the first line and storing the values
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    #Reading the second line and storing the values
    line2 = file.readline()
    label, value = line2.split()
    N_part = float(value)
    #Closing file
    file.close()
    #Storing the rest of the file into 'data'
    data = np.genfromtxt(filename,dtype=None, names=True, skip_header=3)
    return time, N_part, data