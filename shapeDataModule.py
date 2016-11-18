import sys
import os
import math
import numpy
import logging
import galsim
import pyfits
import pandas

# Function to process the COSMOS Catalog header
def procHeaderCosmos(base, index):

	# Determining the COSMOS catalog directory and file name
	# Directory
	if 'dir' in base['input']['cosmos_catalog'][0]:
		dir = base['input']['cosmos_catalog'][0]['dir']
	else:
		dir = os.path.join(galsim.meta_data.share_dir,'COSMOS_25.2_training_sample')
		
	# File Name
	if 'file_name' in base['input']['cosmos_catalog'][0]:
		file_name = base['input']['cosmos_catalog'][0]['file_name']
	else:
		file_name = 'real_galaxy_catalog_25.2.fits'
		
	# Accessing the appropriate HDU from the COSMOS Catalog
	hdu_cat = pyfits.open(os.path.join(dir,file_name))
	
	# Reading the data from the HDU
	cat_data = hdu_cat[1].data
	# COSMOS identifier from the Leauthaud et al. (2010) catalog
	ident = cat_data[index][0]
	# Right ascension and declination (J2000, degrees)
	ra = cat_data[index][1]
	dec = cat_data[index][2]
	# F814W magnitude
	mag = cat_data[index][3]
	
	return (ident, ra, dec, mag)
    #return (0,0,0,0)

class HSMShapeData(galsim.config.ExtraOutputBuilder):

	def initialize(self, data, scratch, config, base, logger):
		
		#print "initialize()"
		
		# Using the base class function
		super(HSMShapeData,self).initialize(data,scratch,config,base,logger)
		
		# Initializing some additional arrays in scratch that will be used to hold the
		# shape parameters and other related information
		self.scratch['g1'] = []
		self.scratch['g2'] = []
		self.scratch['ident'] = []
		self.scratch['dec'] = []
		self.scratch['ra'] = []
		self.scratch['mag'] = []
		self.scratch['flux'] = []
		self.scratch['sigma'] = []
		self.scratch['shape_g1'] = []
		self.scratch['shape_g2'] = []
		self.scratch['shear_in_g1'] = []
		self.scratch['shear_in_g2'] = []
		#self.scratch['g'] = []
		#self.scratch['beta'] = []
		
		# Removing any old output file if present
		#file_name = os.path.join(config['dir'],config['file_name'])
		#if os.path.isfile(file_name):
		#	os.remove(file_name)
		#print config['gal_order']
		#print config['psf_order']
	
	def processImage(self, index, obj_nums, config, base, logger):
	
		#print "processImage()"
	
		# Printing status of run
		if index % 100 == 0:
			print index 
			
		
		# Finding the type of galaxy used for constructing the image
		gal_list = base['gal']['items']
		gal_index = base['gal']['index']
		gal_type = config['gal_tags'][gal_index]
		# Pushing this information into the scratch dictionary for later use
		self.scratch['gal_type'] = gal_type
		
		#print base['gal'] #['index']['current_val']
		#print base['gal']['items'][0]['index']['current_val']
		
		# Extracting Input Shear
		try:
			shear_in =  gal_list[gal_index]['shear']['current_val']
		except:
			shear_in = galsim.Shear(g1 = 0, g2 = 0)
			
		#print shear_in.g1
		#print shear_in.g2
			
		# Storing shear in scratch dictionary
		#self.scratch['g'].append(shear_in.g)
		#self.scratch['beta'].append(shear_in.beta) 
		
		# Extracting value of psf type from scratch to name fields in dataframe 
		psf_list =  base['psf']['items']
		psf_index = base['psf']['index']
		psf_type = config['psf_tags'][psf_index]
		self.scratch['psf_type'] = psf_type
		
		# Finding the type of psf used for constructing the image
		#psf_index = base['psf']['index']
		#psf_type =  psf_list[type_index]['type']
		#self.scratch['psf_type'] = psf_index
		
		# Get the current image from data
		image = base['current_image']
		
		# Drawing the current PSF being used to generate the image 
		if 'psf_image' not in self.scratch:
			#print "Making PSF Image"
			psf = base['psf']['current_val']
			size = base['image']['size']
			pixel_scale = base['image']['pixel_scale']
			psf_image = galsim.Image(size, size)
			psf.drawImage(image=psf_image, scale=pixel_scale)
			self.scratch['psf_image'] = psf_image
				
		# Altering HSM Parameters
		big_hsm = galsim.hsm.HSMParams(max_mom2_iter = 800000)
		
		# Run HSM over the image
		corrected_shear = galsim.hsm.EstimateShear(image, self.scratch['psf_image'], strict=False, shear_est = 'KSB',
							hsmparams = big_hsm) 
		
		# Printing the error message
		#if corrected_shear.error_message != '':
			#print index
			#print corrected_shear.error_message
			
		# Adding the Input Shear Information
		self.scratch['shear_in_g1'].append(shear_in.g1)
		self.scratch['shear_in_g2'].append(shear_in.g2)
			
				
		# Processing according to the image source type
		if gal_type.upper().startswith('COSMOS'):
		
			# COSMOS Galaxy
			# Appending shapeData into the scratch dictionary
			self.scratch['g1'].append(corrected_shear.corrected_g1) #observed_shape.g1)
			self.scratch['g2'].append(corrected_shear.corrected_g2) #observed_shape.g2)
			#print corrected_shear.meas_type
			# Get the general galaxy info from procHeader()
			
			
			# If non-sheared, then the extra information is to be pushed into scratch
			if not gal_type.upper().endswith('SH'):
				cosmos_index = base['gal']['items'][0]['index']['current_val']
				headerInfo = procHeaderCosmos(base, cosmos_index)
				self.scratch['ident'].append(headerInfo[0])
				self.scratch['ra'].append(headerInfo[1])
				self.scratch['dec'].append(headerInfo[2])
				self.scratch['mag'].append(headerInfo[3])
				#print headerInfo[0]
			# If sheared, then pop out the unnecessary information from scratch
			else:
				#print procHeaderCosmos(base, base['image_num'])[0]
				self.scratch.pop('ident', None)
				self.scratch.pop('ra', None)
				self.scratch.pop('dec', None)
				self.scratch.pop('mag', None)
				
			# Popping out the keys that are not needed in this case
			self.scratch.pop('flux', None)
			self.scratch.pop('sigma', None)
			self.scratch.pop('shape_g1', None)
			self.scratch.pop('shape_g2', None)
			
		if gal_type.upper().startswith('GAUSSIAN') :
			# Sheared Gaussian Galaxy
			# Finding the inherent galaxy (Gaussian, in this case) shape
			gal_shape = gal_list[gal_index]['ellip']['current_val']
			gal_flux = gal_list[gal_index]['flux']['current_val']
			gal_sigma = gal_list[gal_index]['sigma']['current_val']
			# Estimated shear
			#estimated_shear = corrected_shear.observed_shape - gal_shape
			# Appending the estimated shear into the scratch dictionary
			#self.scratch['g1'].append(estimated_shear.g1)
			#self.scratch['g2'].append(estimated_shear.g2)
			
			# Pushing all the keys that are needed for a Gaussian galaxy
			self.scratch['flux'].append(gal_flux)
			self.scratch['sigma'].append(gal_sigma)
			self.scratch['shape_g1'].append(gal_shape.g1)
			self.scratch['shape_g2'].append(gal_shape.g2)
			self.scratch['g1'].append(corrected_shear.corrected_g1)
			self.scratch['g2'].append(corrected_shear.corrected_g2)
			#print corrected_shear.meas_type
			#print corrected_shear.correction_status
			
			# Popping out keys that are not required in this case
			self.scratch.pop('ident', None)
			self.scratch.pop('ra', None)
			self.scratch.pop('dec', None)
			self.scratch.pop('mag', None)
				
 		
	def finalize(self, config, base, main_data, logger):
	
		#print "finalize()"
	
		# Extracting the gal_type and psf_type from scratch
		gal_type = self.scratch['gal_type']
		psf_type = self.scratch['psf_type']
		
		# Popping the 'psf_image' from the scratch dictionary as not needed anymore
		self.scratch.pop('psf_image', None)
		self.scratch.pop('gal_type', None)
		self.scratch.pop('psf_type', None)
		
		
		# Extracting value of galaxy type from scratch to name fields in dataframe
		#gal_index = base['gal']['index']
		#gal_type = config['gal_tags'][gal_index]
		
		# Extracting value of psf type from scratch to name fields in dataframe
		#psf_index = base['psf']['index']
		#psf_type = config['psf_tags'][psf_index]

		# Converting the scratch dictionary into a Pandas Dataframe
		# with appropriate naming convention
		
		# Extra Check (Need to be removed later)
		#print "------"
		#print self.scratch['g1'][0:10]
		#print "------"
		#print self.scratch['g2'][0:10]
		#print "------"
		#print self.scratch['shape_g1'][0:10]
		#print "------"
		#print self.scratch['shape_g2'][0:10]
		#print "------"
		
		if gal_type.upper() == 'COSMOS':
			# COSMOS
			columns = ('DEC', gal_type + '_' + psf_type + '_g1', gal_type + '_' + psf_type + 
						'_g2', gal_type +'_IDENT', gal_type + '_MAG', 'RA', 
						gal_type + '_INPUT_g1', gal_type + '_INPUT_g2' )
		
		elif gal_type.upper() == 'COSMOS_SH':
			# COSMOS Sheared
			columns = (gal_type + '_' + psf_type + '_g1', gal_type + '_' + psf_type + '_g2',
						gal_type + '_INPUT_g1', gal_type + '_INPUT_g2')
						
		else:
			# For all other types
			#print self.scratch.keys()
			#sys.exit()
			#columns = (gal_type + '_' + psf_type + '_g1', gal_type + '_' + psf_type + '_g2')
			columns = (gal_type + '_' + psf_type + '_flux', gal_type + '_' + psf_type + '_g1',
						gal_type + '_' + psf_type + '_g2', gal_type + '_' + psf_type + '_shape_g1', 
						gal_type + '_' + psf_type + '_shape_g2', gal_type + '_INPUT_g1',
						gal_type + '_INPUT_g2', gal_type + '_' + psf_type + '_sigma')
			
		shapeData = pandas.DataFrame(self.scratch)
		#print self.scratch['shear_in_g1']
		#print shapeData.keys()
		shapeData.columns = columns
		#print shapeData
		
		return shapeData
		
	def writeFile(self, file_name, config, base, logger):
	
		#print "writeFile()"
		
		# Writing/Appending the shapeData dataframe into a HDF5 file
		# Making the output directory if not present
		if not os.path.isdir(file_name.split('/')[0]):
			os.mkdir(file_name.split('/')[0])
		
		# First Run (New File)
		if not os.path.isfile(file_name):
			hdf = pandas.HDFStore(file_name, mode='w')
			hdf.put('shapeData', self.final_data, format='table', data_columns=True)
			hdf.close()
		# If file already exists in the output directory (Append to File)
		else:
			# Read out the current dataframe from the HDF5 file
			hdf = pandas.HDFStore(file_name)
			dataInit = hdf.get('shapeData')
			# Merge the current dataframe with the dataframe read from file
			cols_to_use = self.final_data.columns.difference(dataInit.columns)
			dataFin = pandas.concat([dataInit, self.final_data[cols_to_use]], axis=1)
			#print dataFin
			hdf.put('shapeData', dataFin, format='table', data_columns=True)
			hdf.close()

		
# Registering the extra output type
galsim.config.RegisterExtraOutput('HSMShapeData', HSMShapeData())