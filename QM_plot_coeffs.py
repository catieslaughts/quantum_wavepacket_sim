import matplotlib.pyplot as plt
import pandas as pd

global datadir, pot_h, pot_w, pot_shape, bg_file, plt_bg, ind_var

#For if the user would like some comparison data plotted in the background:
plt_bg = True
bg_file = './output/coeffs_square_#h11_#w1.csv' #hard-coded in by user

ind_var = 'Height' #the variable that serves as the independent variable, options are "Height" or "Width"


def set_plotting_globals(datadir_in, pot_h_in, pot_w_in, pot_shape_in):
	global datadir, pot_h, pot_w, pot_shape
	
	datadir, pot_h, pot_w, pot_shape = datadir_in, pot_h_in, pot_w_in, pot_shape_in
	
	return

def coeffs_plot():
	data = pd.read_csv(datadir+'coeffs_{}_#h{}_#w{}.csv'.format(pot_shape, len(pot_h), len(pot_w)))
	data.columns =['Height', 'Width', 'Reflection', 'Transmission']
	
	if plt_bg:
		bg_data = pd.read_csv(bg_file)
		bg_data.columns =['Height', 'Width', 'Reflection', 'Transmission']
	
	if ind_var == 'Height':
		for width in pot_w:
			plt_data = data[data['Width'] == width]
			
			if plt_bg:
				bg_plt_data = bg_data[bg_data['Width'] == width]
				plt.plot(bg_plt_data[ind_var], bg_plt_data['Reflection'], c='grey', linestyle = 'dashed', label = 'R (square barrier)')
				plt.plot(bg_plt_data[ind_var], bg_plt_data['Transmission'], c='darkgrey', linestyle = 'dashed', label = 'T (square barrier)')
			
			plt.plot(plt_data[ind_var], plt_data['Reflection'], c='Green', label = 'R ({} barrier)'.format(pot_shape))
			plt.plot(plt_data[ind_var], plt_data['Transmission'], c='Blue', label = 'T ({} barrier)'.format(pot_shape))	
			plt.title('T and R for a {} barrier with width={}'.format(pot_shape, width))
			plt.xlabel('{}'.format(ind_var))
			plt.ylabel('Coefficient Value')
			plt.legend()
		
			plt.savefig(datadir+'coeff_plot_{}_{}_{}.png'.format(pot_shape,ind_var, width), dpi=600)
			plt.show()
		
	if ind_var == 'Width':
		for height in pot_h:
			plt_data = data[data['Height'] == height]
			
			if plt_bg:
				bg_plt_data = bg_data[bg_data['Height'] == height]
				plt.plot(bg_plt_data[ind_var], bg_plt_data['Reflection'], c='grey', linestyle = 'dashed', label = 'R (square barrier)')
				plt.plot(bg_plt_data[ind_var], bg_plt_data['Transmission'], c='darkgrey', linestyle = 'dashed', label = 'T (square barrier)')
			
			plt.plot(plt_data[ind_var], plt_data['Reflection'], c='Green', label = 'R, {}'.format(pot_shape))
			plt.plot(plt_data[ind_var], plt_data['Transmission'], c='Blue', label = 'T, {}'.format(pot_shape))
		
			plt.title('T and R for a {} barrier with height={}'.format(pot_shape, height))
			plt.xlabel('{}'.format(ind_var))
			plt.ylabel('Coefficient Value')
			plt.legend()
		
			plt.savefig(datadir+'coeff_plot_{}_{}_{}.png'.format(pot_shape,ind_var, height), dpi=600)
			plt.show()
	
	