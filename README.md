# COP Assignment 3
### Fall 2022
### Marieke Visscher and Catherine Slaughter

## Description

	A numerical code to simulate quantum wavepacket propagation through various potential barriers

## Dependencies

	This code uses a few relatively basic packages that should be easily installable via your preferred python package management method (conda, pip, etc.)
	- numpy
	- matplotlib
	- scipy 
	- tqdm
	- os

## Usage

	The easiest way to run this simulation code is through the provided wrapper (run_simulation.py). This piece of code allows the user to change input parameter values, and pick which pieces of code to be run. 
	The available input parameters are:

- Simulation properties:
	- tmax (scalar): the number of timesteps in the simulation
    - dt (scalar): the size of timesteps
    - xmax (scalar): the size of the x-axis to be simulated, the code centers the x-axis at zero, so xmax = 200 results in x = [-100,100]
    - dx (scalar): the spacing of the x-grid

- Wavepacket properties:
	- sigma0 (scalar): the width of your Gaussian wavepacket
    - x0 (scalar): the initial center of your wavepacket
    - k0 (scalar): the momentum of the wavepacket

- Potential properties:
	- pot_h (list of scalars): the height(s) of the potential barrier
    - pot_w (list of scalars): the width(s) of the potential barrier, in units of dx
    - pot_shape (string): the shape of the potential barrier

- Animation properties:
	- wavepacket (boolean): whether the user would like the wavepacket plotted in the animation
    - real (boolean): whether the user would like the real part of psi plotted in the animation
    - imaginary (boolean): whether the user would like the imaginary part of psi plotted in the animation

- Wrapper options:
	- run_sim (boolean): to run the simulation with the given parameters
	- animate (boolean): to create an animation based on output files with given parameters
	- datadir (str): the directory that simulation data is saved to and pulled from, directory is created by the code if it doesn't already exist


### Potential barrier options:

Valid values for pot_shape are 'square', 'triangle', 'gaussian', 'step', 'double_well', and 'custom.' 
Each may interpret the input height and width slightly differently:

#### Square: 

A simple square barrier centered at x = 0 with width = pot_w*dx and height = pot_h. If pot_h < 0, a potential well is created

#### Triangle: 

A right-triangle barrier centered at x = 0. The right angle is along with x-axis, at the "left hand" side. The base has length pot_w*dx and the height is pot_h.

#### Gaussian:

A gaussian barrier centered at x = 0 with standard deviation pot_w*dx and height pot_h

#### Step:
	
A step barrier at x = 0 of height pot_h. Note that the step option does not make use of pot_w

#### Double well: 
	
Two "square" barriers centered around x = 0, each of width 2*dx and height pot_h. The inner edges of the barriers are a distance pot_w*dx apart

#### Custom:
	
Should the user want to study the effects of a different potential shape on the movement of a quantum wavepacket, the option exists to hard-code their own.
Where marked in QM_wavepacket_model.py (around line 80), simply set the variable pot equal to some array of equal length to x_pos. By default and as an example, we have input a sine wave potential here
