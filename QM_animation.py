"""
==================
Animated line plot
==================

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os


global pot_w, pot_h, pot_shape, xmax, wavepacket, real, imaginary, wp, psi_real, psi_imag, datadir

def set_anim_globals(pot_w_in, pot_h_in, pot_shape_in, xmax_in, wavepacket_in, real_in, imaginary_in, datadir_in):
	'''To set the global variables as entered in the wrapper code
    	Args: 
    		pot_w_in (list of scalars): the width of the potential barrier 
    		pot_h_in (list of scalars): the height of the potential barrier 
    		pot_shape_in (string): the shape of the potential barrier (square, triangular, Gaussian, step, double_barrier or custom)
    		xmax_in (scalar): the size of the x-axis to be simulated
    		wavepacket_in (boolean): whether the user would like the wavepacket plotted in the animation
    		real_in (boolean): whether the user would like the real part of psi plotted in the animation
    		imaginary_in (boolean): whether the user would like the imaginary part of psi plotted in the animation
    		datadir_in (string): path to the directory of data save files'''
	
	global pot_w, pot_h, pot_shape, xmax, wavepacket, real, imaginary, datadir
	pot_w, pot_h, pot_shape, xmax, wavepacket, real, imaginary = pot_w_in, pot_h_in, pot_shape_in, xmax_in, wavepacket_in, real_in, imaginary_in
	datadir = datadir_in
	
	return

def read_in_data(height, width):
	'''reads in data from output files for plotting. Psi, potential, and x axis data are read in'''
	global psi_array, pot_array, x_array
	
	psi_array = np.loadtxt(datadir+'psi_h{}_w{}_{}.txt'.format(height, width, pot_shape), dtype = 'cdouble')
	pot_array = np.loadtxt(datadir+'pot_h{}_w{}_{}.txt'.format(height, width, pot_shape))
	x_array = np.loadtxt(datadir+'xaxis_{}.txt'.format(xmax))
	
	return
	

def set_up_fig():
	'''sets up plotting animation figure, creates ax objects for the wavepacket, real, and imaginary parts'''
	
	global wp, psi_real, psi_imag, legend

	fig, ax = plt.subplots()
	ax.set_ylim(-1.2,1.2)
	ax.set_xlim(int(.75*x_array.min()), int(.75*x_array.max()))
	ax.plot(x_array,pot_array, c='grey', linestyle = 'dashed', label = 'Potential') #plots potential on all frames
	
	nan_array = np.empty_like(x_array)
	nan_array[:] = np.nan
	#print(nan_array)
	
	if wavepacket:
		wp, = ax.plot(np.copy(nan_array))
		wp.set_color('blue')
		wp.set_label('Wavepacket')
	
	if real:
		psi_real, = ax.plot(np.copy(nan_array))
		psi_real.set_color('orange')
		psi_real.set_label('Real Part')
	
	if imaginary:
		psi_imag, = ax.plot(np.copy(nan_array))
		psi_imag.set_color('green')
		psi_imag.set_label('Imaginary Part')
	
	plt.legend()
	return fig, ax

def animate(frame):
	'''checks what parts of the wave the user wants to see and returns the data to be animated'''
	
	returns = []
	if wavepacket:
		wp.set_ydata(abs(psi_array[frame]))  # update the data.
		wp.set_xdata(x_array)
		returns.append(wp) #adds to returned data
	
	if real:
		psi_real.set_ydata(np.real(psi_array[frame]))  # update the data.
		psi_real.set_xdata(x_array)
		returns.append(psi_real)
	
	if imaginary:
		psi_imag.set_ydata(np.imag(psi_array[frame]))  # update the data.
		psi_imag.set_xdata(x_array)
		returns.append(psi_imag)

	#print(line.get_ydata())
	return returns


def run_animation():
	'''to be called by wrapper code, creates, runs, and saves the animation'''
	print('Creating Animation')
	writer = animation.PillowWriter() #pillowWriter is good for saving Gifs
	anidir = './animations/' #directory to save animations
	if not os.path.exists(anidir): #creating directory if it doesn't exist
		print('Making animation directory')
		os.mkdir(anidir)
	
	for height in pot_h: #iterates over input potential heights and widths
		for width in pot_w:
			print('Height: '+str(height))
			print('Width: '+str(width))
			
			read_in_data(height, width)
			print('Read in data')
			fig, ax = set_up_fig()
			
			print('Animating')
			ani = animation.FuncAnimation(fig, animate, frames = psi_array.shape[0], interval=30, blit=True, repeat = True)
			
			ani.save(anidir+'ani_h{}_w{}_{}.gif'.format(height, width, pot_shape)) #saves animation as gif
			print('Animation Saved')
			plt.show()
