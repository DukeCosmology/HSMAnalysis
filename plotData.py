# Input Format: HDF5_File_Name Mode(0 - Single Var/ 1- Two Var/ 2- Two Var With Truth Value) 
#                                    N_Bins Entry1 <Entry2> <Truth> X-Label

import sys
import matplotlib.pyplot as plt
from pandas import read_hdf
import pandas
from scipy.stats import norm
from scipy.optimize import curve_fit
import numpy as np
import galsim

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x-x0)**2/(2*sigma**2))


def main(argv):

	filename = argv[1]
	mode = int(argv[2])
	nbins = int(argv[3])
	x_label = ""
	data = []
	
	# Reading the HDF5 File
	hdf = read_hdf(filename,'shapeData')
	
	# Plotting
	plt.figure(1)
	
	if mode == 0:
		# Single Variable
		data = hdf[argv[4]]
		x_label = argv[5]
		
	if mode == 1:
		# Double Variable
		data = (hdf[argv[4]] - hdf[argv[5]])#/((hdf[argv[5]] + hdf[argv[4]])/2)
		#data = hdf[argv[4]] - hdf[argv[5]]
		#var1 = hdf[argv[4]]
		#var2 = hdf[argv[5]]
		
		# Calculating fractional error for each pair
		#for index in range(len(var1)):
			#data.append((var1[index] - var2[index])/((var1[index] + var2[index])/2))
		
		x_label = argv[6]
		print max(data)
		#plt.xlim(0.02,0.04)
		
	if mode == 2:
		#plt.xlim(-0.25,0.25)
		# Comparison against a truth value
		diff = (hdf[argv[4]] - hdf[argv[5]]).values.tolist()
		truth = float(argv[6]) # Single Value
		# Find fractional error using this single 'truth' value for all values in 'est'
		outliers = []
		for val in diff:
			# Calculating fractional error and appending
			frac_error = (val - truth)/((val + truth)/2)
			data.append(frac_error)
			
			# Some extra checks for tail
			if frac_error <= -0.1:
				outliers.append(diff.index(val))
			
			#if frac_error > 1 or frac_error < -1:
				
				#index = diff.index(val)
				#print index
				#print "g2 frac error:"
				#print frac_error
				#print hdf[argv[4]][index]
				#print "g1 frac error:"
				#print (hdf['Gaussian_LSST_g1'][index] - (-0.027))/((hdf['Gaussian_LSST_g1'][index] + (-0.027))/2)
				#print hdf[argv[5]][index]
				#print frac_error
				#print "\n"
				
			#else:
				#data.append(frac_error)
			
		# Writing the images of all the outliers into a separate data cube
		datacube = galsim.fits.readCube('gaussian_lsst.fits')
		img_outliers = []
		for item in outliers:
			img_outliers.append(datacube[item])
		
		galsim.fits.writeCube(img_outliers, 'outliers_g2.fits')
			
		# Writing the data for all the outliers into a separate dataframe
		hdf_outlier = hdf.iloc[outliers]
		df = pandas.HDFStore('outliers_g2.h5', mode='w')
		df.put('shapeData', hdf_outlier, format='table', data_columns=True)
		df.close()
			
		
				
			
		x_label = argv[7]
		
	# Normal fit to the data
	print "Normal Fit :::"
	(mu, sigma) = norm.fit(data)
	print('Mean (Mu): %f' % mu)
	print('Standard Deviation (Sigma): %f' % sigma)
	
	# Normal fit to the data
	#xmin = min(data)
	#xmax = max(data)
	#x_axis = np.arange(xmin, xmax, 0.001)
	#plt.plot(x_axis, norm.pdf(x_axis, mu, sigma))
	
	# Histogram
	(n, bins, patches) = plt.hist(data, color='k', alpha=0.5, bins=nbins, fill = False, histtype='step')
	plt.xlabel(x_label)
	
	# Correcting the bins
	correction = (bins[1] - bins[0])/2
	bins = [(bin_val + correction) for bin_val in bins]
	bins = bins[:nbins]
	
	# Fitting a Gaussian
	try:
		popt, pcov = curve_fit(gaussian, bins, n) #, p0 = [1,0.00016,1]) #, p0 = [200, 25000, 10000])
		#plt.plot(x_axis, gaussian(x_axis, *popt), label='fit')
		error = np.sqrt(np.diag(pcov))
		print "Gaussian Fit :::"
		print popt
		print error
		
		# Generating the fractional difference plot for the fit
		#true_mean = 5.e4
		#true_sigma = 2.e4
		#frac_diff = []
		#for index in range(len(bins)):
			#truth = gaussian(bins[index], 5000, true_mean, true_sigma)
			#recon = n[index]
			#frac_diff.append((recon - truth)/((recon + truth)/2))
			
		#plt.figure(2)
		#plt.hist(frac_diff, color='k', alpha=0.5, bins=nbins, fill = False, histtype='step')
	
	except RuntimeError:
		print "Gaussian Fit Failed!"
	
	# Display
	plt.show()
	
	# Print histogram data
	#print n
	#print bins

	
	
if __name__ == "__main__":
    main(sys.argv)