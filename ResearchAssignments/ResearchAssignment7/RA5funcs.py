#Functions used for ReasearchAssignment 5
# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


#import previously developed functions
#Allows us to get the mass of each seperate component for a galaxy
from ComponentMassFuncs import ComponentMass
#Reading the .txt files
from ReadFile import Read
#Allows us to get the COM position and velocity
from CenterOfMass import CenterOfMass
from RelativeMagnitude import RelativeMagnitude
from RotateFrame import RotateFrame


#Setting up the masses
M_host = ComponentMass('../../../M31/M31_000.txt', 1) + ComponentMass('../../../M31/M31_000.txt', 2) + ComponentMass('../../../M31/M31_000.txt', 3)+ ComponentMass('../../../MW/MW_000.txt', 1) + ComponentMass('../../../MW/MW_000.txt', 2)+ ComponentMass('../../../MW/MW_000.txt', 3)
#M_host has to be changing with time since MW and M31 are merging
M_sat = ComponentMass('../../../M33/M33_000.txt', 1) + ComponentMass('../../../M33/M33_000.txt', 2)
#M_sat is also changing with time since M33 is losing mass every time particles get beyond the Jacobi Radius
#Getting the orbital info of M31 and M33
M31_Orbit = np.genfromtxt("Orbit_M31.txt")
M33_Orbit = np.genfromtxt("Orbit_M33.txt")
M31_M33_Sep = RelativeMagnitude(M33_Orbit[:, 1:4], M31_Orbit[:, 1:4])
M31_M33_Vel = RelativeMagnitude(M33_Orbit[:, 4:], M31_Orbit[:, 4:])


def JacobiRadius(M_sat, M_host, r):
    '''
    This function will calculate the Jacobi radius for a given galaxy at a given moment

    inputs:
        M_sat: Mass of the satelite
        M_host: Mass of the host
        r: Separation between host and satelite

    outputs;
        R_j: Jacobi Radius
    '''
    R_j = r*(M_sat/(2*M_host))**(1/3)
    return R_j


def VelocityDisp(filenumber):
    '''
    This function uses code from Lab7 to plot the M33 disk particles color coded by velocity
    '''
    #Code from Lab 7 for contour and scatter plots
    COMD = CenterOfMass(f'../../../M33/M33_{filenumber}.txt',2)
    
    #This file is equivalent to the M33_Orbit[52] (because M33_Orbit has step=5)
    COMP = COMD.COM_P(0.1)
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 
    
    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)
    
    # Determine velocities of disk particles relative to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 
    
    # total velocity 
    vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)
    
    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T
    
    #Rotate frame
    rn, vn = RotateFrame(r, v)
    
    #Arrays for Jacobi radius
    magnitudes = np.sqrt(np.sum(rn**2, axis=1))
    indx = filenumber/5.0
    threshold = np.array([JacobiRadius(M_sat, M_host, M31_M33_Sep[int(indx)])])
    print(threshold)
    
    # Create a mask for positions that meet the condition
    mask = magnitudes > threshold
    
    # Apply the mask to filter positions and corresponding velocities
    r_J = rn[mask]
    v_J = vn[mask]

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.scatter(rn[:,1], rn[:,2], c=vn[:,0])
    ax.add_patch(plt.Circle((0, 0), threshold, color='r', fill=False))

    #colorbar
    cbar = plt.colorbar()
    cbar.set_label('Vx (km/s)', size=22)
    
    # Add axis labels
    plt.xlabel('y (kpc)', fontsize=22)
    plt.ylabel('z (kpc)', fontsize=22)
    
    
    
    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    
    plt.show()

def VelocityDispMask(filenumber):
    '''
    This function uses code from Lab7 to plot the M33 disk particles color coded by velocity
        plus it masks particles within the jacobi radius for any given snapshot
    
    Inputs:
        filenumber [int]: number of the snapshot M33 file that will be used
    
    Outputs:
        r_J [array]: Array containing the positions of particles beyond the Jacobi Radius
        v_J [array]: Array containing the velocity vectors of particles beyond the Jacobi Radius
        threshold [array]: Array containing the Jacobi Radius for that specific snapshot
        magnitudes [array]: Array containing the magnitudes of the positions 
    '''
    #Code from Lab 7 for contour and scatter plots
    COMD = CenterOfMass(f'../../../M33/M33_{filenumber}.txt',2)

    #COMD for M31 - Used for velocity dispersion of stellar streams compared to M31 COM motion
    COMD2 = CenterOfMass(f'../../../M31/M31_{filenumber}.txt',2)
    
    #This file is equivalent to the M33_Orbit[52] (because M33_Orbit has step=5)
    COMP = COMD.COM_P(0.1)
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

    COMP2 = COMD2.COM_P(0.1)
    COMV2 = COMD2.COM_V(COMP2[0], COMP2[1], COMP2[2])
    
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 
    
    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)
    
    # Determine velocities of disk particles relative to COM motion
    vxD = COMD.vx - COMV2[0].value 
    vyD = COMD.vy - COMV2[1].value 
    vzD = COMD.vz - COMV2[2].value 
    
    # total velocity 
    vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)
    
    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T
    
    #Rotate frame
    rn, vn = RotateFrame(r, v)
    
    #Arrays for Jacobi radius
    magnitudes = np.sqrt(np.sum(rn**2, axis=1))
    indx = filenumber/5.0
    threshold = np.array([JacobiRadius(M_sat, M_host, M31_M33_Sep[int(indx)])])
    
    # Create a mask for positions that meet the condition
    mask = magnitudes > threshold
    
    # Apply the mask to filter positions and corresponding velocities
    r_J = rn[mask]
    v_J = vn[mask]

    return r_J, v_J, threshold, magnitudes

def PlotVelDisp(filenumber):
    '''
    Function that plots the velocity coded scatter plot for a given snapshot

    Inputs:
        filenumber [int]: Snapshot that will have its velocity graphed

    Outputs:
        Scatter plot coded by velocity with all particles within the Jacobi Radius blocked out, 
        and with a 2D circle representing the Jacobi Radius for that snapshot
    
    '''
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.scatter(VelocityDispMask(filenumber)[0][:,1], VelocityDispMask(filenumber)[0][:,2], c=VelocityDispMask(filenumber)[1][:,0])
    ax.add_patch(plt.Circle((0, 0), VelocityDispMask(filenumber)[2], color='r', fill=False))
    
    #colorbar
    cbar = plt.colorbar()
    cbar.set_label('Vx (km/s)', size=22)
    
    # Add axis labels
    plt.xlabel('y (kpc)', fontsize=22)
    plt.ylabel('z (kpc)', fontsize=22)
    
    
    
    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    
    plt.show()

def PositionvsDispersion(filenumber, window_size):
    '''
    Function that calculates the Dispersion for particles beyond the jacobi Radius, for any given distance from the center of M33
    
    Inputs:
        filenumber [int]: Snapshot that will be used
        window_size [int]: Window_size to be used for the dispersion
    
    Outputs:
        r_dispersion_centers [array]: Array containing the sorted distances of particles from the center of M33
        dispersion [array]: Array containing the velocity dispersion 
    '''
    # Compute magnitude of velocity at each time step
    v_magnitude = np.linalg.norm(VelocityDispMask(filenumber)[1], axis=1)
    r_magnitude = np.linalg.norm(VelocityDispMask(filenumber)[0], axis=1)
    
    # Use a moving window over radius to compute local dispersion
    # Sort by radius
    sorted_indices = np.argsort(r_magnitude)
    r_sorted = r_magnitude[sorted_indices]
    v_sorted = v_magnitude[sorted_indices]
    
    # Define window size in terms of number of points
    dispersion = np.array([
        np.std(v_sorted[i:i + window_size])
        for i in range(len(r_sorted) - window_size)
    ])
    r_dispersion_centers = np.array([
        np.mean(r_sorted[i:i + window_size])
        for i in range(len(r_sorted) - window_size)
    ])


    return r_dispersion_centers, dispersion

def AvgDisp(start, end, window_size):
    '''
    Function that computes the average dispersion for multiple snapshots
    
    Inputs:
        start [int]: snapshot to start 
        end [int]: snapshot to end
        window_size [int]: window size to be used for the dispersion claculation
    
    Outputs:
        dispersions [array]: An array containing the mean velocity dispersion for each given snapshot
            (in implements of 5 so that it matches the "M33_Orbit" framework)
    '''
    num_steps = (end - start) // 5
    dispersions = np.zeros(num_steps)

    for i, index in enumerate(range(start, end, 5)):
        dispersions[i] = np.mean(PositionvsDispersion(index, window_size))

    return dispersions
