#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tunneling simulation code
"""
import numpy as np
from scipy.sparse import diags
from numpy.linalg import inv
from scipy.signal import gaussian
from tqdm import tqdm

global tmax, dt, xmax, dx, sigma0, x0, k0, pot_h, pot_w, pot_shape, pi, datadir

def set_sim_globals(tmax_in, dt_in, xmax_in, dx_in, sigma0_in, x0_in, k0_in, pot_h_in, pot_w_in, pot_shape_in, datadir_in):
    """To set the global variables as entered in the wrapper code
    Args: 
        tmax_in (scalar): the number of timesteps in the simulation
        dt_in (scalar): the size of timesteps
        xmax_in (scalar): the size of the x-axis to be simulated
        dx_in (scalar): the spacing of the x-grid
        sigma0_in (scalar): the width of your Gaussian wavepacket
        x0_in (scalar): the (starting) center of your wavepacket
        k0_in (scalar): the momentum of the wavepacket
        pot_h_in (list of scalars): the height of the potential barrier
        pot_w_in (list of scalars): the width of the potential barrier
        pot_shape_in (string): the shape of the potential barrier (square, triangular, Gaussian, step, double_barrier or custom)
        datadir_in (string): the directory to save the files"""
    global tmax, dt, xmax, dx, sigma0, x0, k0, pot_h, pot_w, pot_shape, pi, datadir
    tmax, dt, xmax, dx, sigma0, x0, k0, pot_h, pot_w, pot_shape = tmax_in, dt_in, xmax_in, dx_in, sigma0_in, x0_in, k0_in, pot_h_in, pot_w_in, pot_shape_in
    datadir = datadir_in
	
    pi = np.pi

    return


def set_up_system(height, width):
    """To initialize the size, timesteps, potential, hamiltonian and wavefunction of the system.
    Returns:
        x_pos (array): grid of x-positions
        potential (array): potential barrier
        psi_o (array): the initial wavefunction
        hamiltonian (array): the hamiltonian in matrix form"""
    x_pos = np.linspace(-xmax/2, xmax/2, int(xmax/dx))
    psi_o = get_psi_init(x_pos)
    potential = create_pot(x_pos, height, width)
    hamiltonian = create_ham(xmax, dx, potential)

    return x_pos, potential, psi_o, hamiltonian

def create_pot(x_pos, height, width):
    """To create the potential barrier. The barrier is always created at the center of the grid.
    Args:
        x_pos (array): grid of x-positions
        shape (string): the shape of the barrier. Either square, triangle or Gaussian.
    Returns:
        pot (array): the potential barrier"""
    pot = np.zeros_like(x_pos)
    center_val = int((xmax/dx)/2)
    start = int(center_val - (width/2))
    stop = int(center_val + (width/2))	
    if pot_shape == 'square':  
        pot = np.zeros_like(x_pos)
        pot[start:stop] = height
        
    elif pot_shape == 'triangle':
        pot[start:stop] = height/2 - (height/width)*x_pos[start:stop]
        
    elif pot_shape == 'gaussian':
        pot = height*gaussian(xmax/dx, width)
        
    elif pot_shape == 'step':
        pot[center_val+1:] = height
    
    elif pot_shape == 'double_barrier':
        pot[start-2:start] = height
        pot[stop:stop+2] = height
    
    #=========CUSTOM POTENTIAL BARRIER=========    
    elif pot_shape == 'custom':
    #IF THE USER WOULD LIKE TO USE THEIR OWN CUSTOM POTENTIAL, IT MAY BE HARD-CODED HERE
    	pot=height*np.sin(x_pos* (2*pi/width))
    #=========================================
    
    return pot

def create_ham(xmax, dx, pot):
    """To calculate the Hamiltonian (in matrix form).
    Args:
        xmax (scalar): the length of the position array
        dx (scalar): the spacing between two x-values
        pot (array): the potential barrier
    Returns: the Hamiltonian (H = p^2/2m + V)."""
    A1 = np.tile(2, int(xmax/dx))
    A2 = np.tile(-1, int(xmax/dx))
    Ekin = diags([A1, A2, A2], [0, -1, 1]).toarray()
    Epot = diags([pot], [0]).toarray()
    return Ekin/(2*dx**2) + Epot

def get_psi_init(x_pos):
    """To generate the intital wavefunction, which is assumed to be a Gaussian wavepacket.
    Args:
        x_pos (array): grid of x-positions
    Returns:
        psi_o (array): the wavefunction at t = 0
        """
    psi_o = (2*pi*sigma0**2)**(-1/4)*np.exp(-(x_pos-x0)**2/(4*sigma0**2))*np.exp(1j*k0*x_pos)

    return psi_o

def update_system(psi_old, ham):
    """To calculate the wavefunction at the next timestep according to the Schroedinger equation.
    Args:
        psi_old (array): the wavefunction of the previous timestep
        ham (array): the hamiltonian
    Returns:
        psi_new_norm (array): the (normalized) new wavefunction"""
    CN_mat = np.matmul(inv(np.identity(int(xmax/dx)) - dt*ham/(2j)),(np.identity(int(xmax/dx)) + dt*ham/(2j)))
    psi_new = np.dot(CN_mat, psi_old)
    psi_new_norm = psi_new/np.sqrt(np.trapz(abs(psi_new)**2, dx = dx)) #renormalize the wavefunction
    
    return psi_new_norm

def coefficients(xmax, psi_new, dx, height, width):
    """To calculate the transmission and reflection coefficients after the interaction with the potential barrier.
    Args:
        xmax (scalar): the length of the position array
        psi_new (array): the wavefunction after interaction with the barrier
        dx (scalar): the spacing between two x-values
    Returns:
        R_coeff (scalar): the reflection coefficient
        T_coeff (scalar): the transmission coefficient"""
    T_coeff = np.trapz(np.abs(psi_new[int(xmax/dx/2+width/2)::])**2, dx = dx)
    R_coeff = np.trapz(np.abs(psi_new[:int(xmax/dx/2-width/2):])**2, dx = dx)
    return R_coeff, T_coeff


def run_system():
    """To run the simulation."""
    print('Running Sim')
    
    coeff_file = open(datadir+'coeffs_{}_#h{}_#w{}.csv'.format(pot_shape, len(pot_h), len(pot_w)), 'w' )
    coeff_file.write('#Height,Width,Reflection,Transmission\n')
    coeff_file.close()
    coeff_file = open(datadir+'coeffs_{}_#h{}_#w{}.csv'.format(pot_shape, len(pot_h), len(pot_w)).format(pot_shape), 'a' )
    
    for height in pot_h: #iterate over input potential heights and widths
        for width in pot_w:
        	print('Height: '+str(height))
        	print('Width: '+str(width))
        	x_pos, potential, psi_o, hamiltonian = set_up_system(height,width) #establish potential and initial wavefunction
        	psi_new = np.zeros((int(tmax/dt), int(xmax/dx)), dtype = np.clongdouble) # to save the wavefunction
        	psi_new[0] = np.copy(psi_o) #copy initial wavefunction
        	for i in tqdm(range(1, int(tmax/dt))):
        		psi_new[i] = update_system(psi_new[i-1], hamiltonian) #step system, save new wavefunction
        		if i == int(tmax/dt)-1:
        			R_coeff, T_coeff = coefficients(xmax, psi_new[i], dx, height, width) #calculates reflection and transmission
        	print("Reflection: "+str(R_coeff))
        	print("Transmission: "+str(T_coeff))
        	
        	#coeffs_save = np.asarray([R_coeff, T_coeff])
        	
        	coeff_file.write('{},{},{},{}\n'.format(height, width, R_coeff, T_coeff))
        	np.savetxt(datadir+'psi_h{}_w{}_{}.txt'.format(height, width, pot_shape), psi_new) #saves data
        	np.savetxt(datadir+'pot_h{}_w{}_{}.txt'.format(height, width, pot_shape), potential)
        	np.savetxt(datadir+'xaxis_{}.txt'.format(xmax), x_pos)
        	print('Data saved')
        	
    coeff_file.close()


