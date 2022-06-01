import numpy as np
from QM_wavepacket_model import *
from QM_animation import *
from QM_plot_coeffs import *
import os

global tmax, dt, xmax, dx, sigma0, x0, k0, pot_h, pot_w, pot_shape, wavepacket, real, imaginary, datadir

#=================USER INPUT====================
#simulation properties
tmax = 200
dt = 1
xmax = 200
dx = 1

#wavepacket properties
sigma0 = 5.0
k0 = 1.0
x0 = -80

#potential properties
pot_h = [.5]#np.arange(.5,5.5,.5) #must be iterable
pot_w = [10] #must be iterable
pot_shape = 'gaussian' #options are 'square', 'triangle', 'gaussian', 'step', 'double_barrier', and 'custom'

#animation properties
wavepacket = True
real = True
imaginary = False

#Wrapper options:
run_sim = True #runs the simulation with the given parameters
animate = True #creates an animation based on output files with given parameters
plot_coeffs = False #plots the reflection and transmission coefficients
datadir = './output/'

#================END USER INPUT====================

shape_options = ['square', 'triangle', 'gaussian', 'step', 'double_barrier', 'custom']

if not os.path.exists(datadir): #create output directory, if it doesn't exist
	print('Making output directory')
	os.mkdir(datadir)


if pot_shape in shape_options: #check that input potential shape is a valid option
	if run_sim: #runs the simulation
		set_sim_globals(tmax, dt, xmax, dx, sigma0, x0, k0, pot_h, pot_w, pot_shape, datadir)
		run_system()

	if animate: #creates the animation
		set_anim_globals(pot_w, pot_h, pot_shape, xmax, wavepacket, real, imaginary, datadir)
		run_animation()
	
	if plot_coeffs: #plots T and R coefficients
		set_plotting_globals(datadir, pot_h, pot_w, pot_shape)
		coeffs_plot()
		
else: #if input potential shape is not valid
	print("Invalid potential shape. Options are 'square', 'triangle', 'gaussian', 'step', 'double_well', or 'custom'")