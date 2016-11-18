import sys
import os
import shutil
import datetime
import logging
import galsim

class FileCleaner(galsim.config.OutputBuilder):
		
	def writeFile(self, data, file_name):
		# No files need to be written here
		# This function cleans up the output area by transferring all the previously
		# written files into a time coded backup folder
		
		# Making the output directory if not present
		if os.path.isdir(os.path.join(os.getcwd(), file_name.split('/')[0])):
			
			# Checking if the output directory is empty
			if not os.listdir(os.path.join(os.getcwd(), file_name.split('/')[0])) == "":
				
				# Output directory is not empty
				# Making a new time coded backup folder in the output directory
				backup_dir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    			#os.makedirs(backup_dir)
    			
    			# Moving all the files from the output directory to the backup directory
    			shutil.move(os.path.join(os.getcwd(),file_name.split('/')[0]), backup_dir)
		
		
	def canAddHdus(self):
		# Don't add extra HDUs
		return False

		
# Registering the extra output type
galsim.config.RegisterOutputType('FileCleaner', FileCleaner())













# 
# 
# 
# 
# def dataCatalog(cat_filename, cat_dir, pos):
#     # Accessing the appropriate HDU from the COSMOS Catalog
#     hdu_cat = pyfits.open(os.path.join(cat_dir,cat_filename))
#     # Reading the data from the HDU
#     cat_data = hdu_cat[1].data
#     # COSMOS identifier from the Leauthaud et al. (2010) catalog
#     ident = cat_data[pos][0]
#     # Right ascension and declination (J2000, degrees)
#     ra = cat_data[pos][1]
#     dec = cat_data[pos][2]
#     # F814W magnitude
#     mag = cat_data[pos][3]
# 
#     return (ident, ra, dec, mag)
# 
# 
# def HSMShapeData(config, base):
# 
#     # Parsing the configuration
# 
#     logging.basicConfig(
#         format="%(message)s",
#         level=logging.INFO,
#         stream=sys.stdout)
#     logger = logging.getLogger("HSM Estimation")
#     
#     # Reading the configuration file using PyYAML
#     stream = file(argv[1], 'r')
#     config = yaml.load(stream)
# 
#     # Defining parameters used below
# 
#     # Big FFT parameters (Adjusting Max FFT Size)
#     big_fft = galsim.GSParams(maximum_fft_size=24 * 4096)
#     #big_hsm = galsim.hsm.HSMParams(max_mom2_iter=400000,failed_moments=-22.0)
# 
#     cat_filename = 'real_galaxy_catalog_25.2.fits'
#     # dir = 'data'
# 
#     # Make output directory if not already present.
#     if not os.path.isdir('output'):
#         os.mkdir('output')
# 
#     # Galaxy Shear Parameters
#     # Using e1-e2 shear
#     gal_e1 = 0.027
#     gal_e2 = 0.031
# 
#     # Number of galaxies to process
#     try:
#     	ngal = int(config['HSM_Comparison']['ngal'])
#     	
#     except (ValueError, TypeError, KeyError):
#     	logger.info('Argument Missing/Wrong Type!')
#     	logger.info('Need to specify an integer number of galaxies.')
#     	sys.exit()
#     
#     # LSST Filter Used
#     if 'filter' in config['HSM_Comparison']:
#     	filter = config['HSM_Comparison']['filter']
#     else:
#     	filter = 'r'
#     
#     # HSM Non-Convergence Error Code
#     hsm_error_code = -1000.
# 
#     # Output Storage
#     # Galaxy Information
#     ident = []
#     ra = []
#     dec = []
#     mag = []
#     # Images
#     im_hst_list = []
#     im_hst_sh_list = []
#     im_lsst_list = []
#     im_lsst_sh_list = []
#     # HST Data
#     e1_hst = []
#     e2_hst = []
#     e1_hst_sh = []
#     e2_hst_sh = []
#     amp_hst = []
#     amp_hst_sh = []
#     # LSST Data
#     e1_lsst = []
#     e2_lsst = []
#     e1_lsst_sh = []
#     e2_lsst_sh = []
#     amp_lsst = []
#     amp_lsst_sh = []
# 
#     # Image Parameters
#     im_size = 64
#     pixel_scale = 0.2      # arcsec
#     sky_level = 1.e6        # ADU / arcsec^2
# 
#     logger.info('Starting HSM-based Shear Testing using the following parameters:')
#     logger.info('    - real galaxies from catalog %r', cat_filename)
#     logger.info('    - pixel scale = %.2f', pixel_scale)
#     logger.info('    - filter = %s', filter)
#     logger.info('    - Applied gravitational shear = (%.3f,%.3f)',gal_e1,gal_e2)
# 
#     # Read in galaxy catalog
#     if 'catalog_directory' in config['HSM_Comparison']:
#     	cat_dir = config['HSM_Comparison']['catalog_directory']
#         cosmos_cat = galsim.COSMOSCatalog(file_name = cat_filename, dir = cat_dir)
#     else:
#         cat_dir = os.path.join(galsim.meta_data.share_dir,'COSMOS_25.2_training_sample')
#         cosmos_cat = galsim.COSMOSCatalog()
# 
#     logger.info('Reading in real galaxies from catalog')
# 
#     # Making the PSFs
#     # HST PSF
#     # bigger than the actual HST F814W PSF/ diam(m) and lam(nm)
#     psf_hst = galsim.ChromaticOpticalPSF(lam=1000.0, diam=2.4, aberrations = numpy.zeros(12))
# 
#     # LSST PSF (Parameters based on the lsst.yaml)
#     # Atmospheric Component (Type: Kolmogorov) (Chromatic)
#     fwhm = 0.67 # Value depends on wavelength of monochromatic PSF
#     # Constructing atmospheric PSF
#     base_atmos = galsim.Kolmogorov(fwhm=fwhm) # 500 nm
#     lsst_atmos = galsim.ChromaticAtmosphere(base_atmos, 500, zenith_angle=0 * galsim.degrees)
# 
#     # Optical Component
#     lam = 700 # nanometers
#     diam = 8.4 # meters
#     lam_over_diam = 0.017 # (700 nm / 8.4 m) * 206265 = 0.017 arcsec
#     obscuration = 0.4  # linear obscuration (3.4m/8.4m = 0.4)
#     nstruts = 4  # number of supports
#     strut_thick = 0.03  # relative to dimeter (25cm in this case)
#     # Arbitrary, intentionally not aligned with the axis of the image
#     strut_angle = 10 * galsim.degrees
#     # Plausible Values of aberrations
#     # Using Noll Convention
#     aberrations = numpy.zeros(12)
#     aberrations[4] = 0.06 # index 4 (defocus)
#     aberrations[5] = 0.02 # index 5 (astig1)
#     aberrations[6] = -0.03 # index 6 (astig2)
#     aberrations[7] = 0.02 # index 7 (coma1)
#     aberrations[8] = -0.04 # index 8 (coma2)
#     # Constructing optical PSF
#    #  lsst_optics = galsim.OpticalPSF(
# #         lam_over_diam=lam_over_diam,
# #         defocus=defocus,
# #         astig1=astig1,
# #         astig2=astig2,
# #         coma1=coma1,
# #         coma2=coma2,
# #         obscuration=obscuration,
# #         nstruts=nstruts,
# #         strut_thick=strut_thick,
# #         strut_angle=strut_angle)
# 
#     # Defining a chromatic LSST PSF
#     lsst_optics = galsim.ChromaticOpticalPSF(lam = lam, diam = diam, aberrations = aberrations)
# 
#     # Constructing the composite LSST PSF based on the optical and atmospheric
#     # components
#     psf_lsst = galsim.Convolve([lsst_atmos, lsst_optics])
#     
#     # HST Filter (F814W)
#     bp_file = os.path.join(galsim.meta_data.share_dir, 'wfc_F814W.dat.gz')
#     filter_hst = galsim.Bandpass(bp_file, wave_type='ang').thin().withZeropoint(25.94)
#     
#     # Defining the LSST Filters
#     try:
#     	filter_names = 'ugrizy'
#     	filter_lsst = {}
#     	for filter_name in filter_names:
#         	filter_filename = os.path.join(config['HSM_Comparison']['filter_data'], 'LSST_{}.dat'.format(filter_name))
#         	filter_lsst[filter_name] = galsim.Bandpass(filter_filename, 'nm')
#         	filter_lsst[filter_name] = filter_lsst[filter_name].thin(rel_err=1e-4)
#     
#     except (ValueError, TypeError, KeyError):
#     	logger.info('Argument Missing/Wrong Type!')
#     	logger.info('Need to specify the location of LSST Filters.')
#     	sys.exit()
#     
#     logger.info('Read in filters')
# 
# 
#     # Build the images from the catalog
#     for k in numpy.arange(ngal):
#     
#         logger.info('Estimating Parameters For Galaxy #%d', k)
# 
#         # Initializing the random number generator
#         rng = galsim.UniformDeviate()
# 
#         # Constructing the galaxy
#         # Since, no index is specified, galaxies are chosen at random from the
#         # catalog
#         gal_hst = cosmos_cat.makeGalaxy(k, gal_type='parametric', chromatic = 'True', gsparams=big_fft)
# 
#         # Flux Scaling for LSST
#         # ratio of (LSST diam)^2/(HST diam)^2
#         diam_sq_ratio = (8.4**2) / (2.4**2)
#         # exposure time of LSST in secs compared to 1 sec for HST
#         exp_time = 16
#         gal_lsst = gal_hst.withScaledFlux(flux_ratio=diam_sq_ratio * exp_time)
# 
#         # Applying transformations to galaxies
#         # Rotation, Shear, Magnification
#         # Rotating by a random angle
#         # theta = 2.*math.pi * rng() * galsim.radians
#         # gal_hst = gal_hst.rotate(theta)
#         # gal_lsst = gal_lsst.rotate(theta)
# 
#         # Applying the desired shear
#         gal_hst_sh = gal_hst.shear(e1=gal_e1, e2=gal_e2)
#         gal_lsst_sh = gal_lsst.shear(e1=gal_e1, e2=gal_e2)
# 
#         # Make the combined profile
#         # Non-sheared
#         final_hst = galsim.Convolve([psf_hst, gal_hst], gsparams=big_fft)
#         final_lsst = galsim.Convolve([psf_lsst, gal_lsst], gsparams=big_fft)
#         # Sheared
#         final_hst_sh = galsim.Convolve([psf_hst, gal_hst_sh], gsparams=big_fft)
#         final_lsst_sh = galsim.Convolve(
#             [psf_lsst, gal_lsst_sh], gsparams=big_fft)
#         # pdb.set_trace()
# 
#         logger.info('\tGalaxy Construction Done...')
#         #logger.info(time.time())
# 
#         # Offset by up to 1/2 pixel in each direction
#         # We had previously (in demo4 and demo5) used shift(dx,dy) as a way to
#         # shift the center of the image.  Since that is applied to the galaxy,
#         # the units are arcsec (since the galaxy profile itself doesn't know
#         # about the pixel scale).  Here, the offset applies to the drawn image,
#         # which does know about the pixel scale, so the units of offset are
#         # pixels, not arcsec.  Here, we apply an offset of up to half a pixel
#         # in each direction.
#         # dx = rng() - 0.5
#         # dy = rng() - 0.5
# 
#         # Drawing the galaxy profiles
#         # HST
#         # Non-sheared
#         im_hst = galsim.Image(im_size, im_size)
#         final_hst.drawImage(filter_hst, image=im_hst, scale=pixel_scale)
#         # Sheared
#         im_hst_sh = galsim.Image(im_size, im_size)
#         final_hst_sh.drawImage(filter_hst, image=im_hst_sh, scale=pixel_scale)
# 
#         # LSST
#         # Non-sheared
#         im_lsst = galsim.Image(im_size, im_size)
#         final_lsst.drawImage(filter_lsst[filter], image=im_lsst, scale=pixel_scale)
#         # Sheared
#         im_lsst_sh = galsim.Image(im_size, im_size)
#         final_lsst_sh.drawImage(filter_lsst[filter], image=im_lsst_sh, scale=pixel_scale)
# 
#         logger.info('\tImages Rendered...')
#         #logger.info(time.time())
# 
#         # Add a constant background level
#         # background = sky_level * pixel_scale**2
#         # im += background
# 
#         # Add Poisson noise.  This time, we don't give a sky_level,
#         # since we have already added it to the image, so we don't want
#         # any more added.  The sky_level parameter really defines how much
#         # _extra_ sky should be added above what is already in the image.
#         # im.addNoise(galsim.PoissonNoise(rng))
# 
#         # Store that into the list of all images
#         # all_images += [im_real]
#         # t5 = time.time()
# 
#         # Creating the output directory if not available
#         if not os.path.isdir('output'):
#             os.mkdir('output')
# 
#         # Writing the HST and LSST Images
#         # Storing a list of output images
#         im_hst_list += [im_hst]
#         im_hst_sh_list += [im_hst_sh]
#         im_lsst_list += [im_lsst]
#         im_lsst_sh_list += [im_lsst_sh]
# 
#         # Running HSM estimation
#         shape_hst = im_hst.FindAdaptiveMom(strict=False)
#         shape_hst_sh = im_hst_sh.FindAdaptiveMom(strict=False)
#         shape_lsst = im_lsst.FindAdaptiveMom(strict=False)
#         shape_lsst_sh = im_lsst_sh.FindAdaptiveMom(strict=False)
#         # pdb.set_trace()
#         
#         
#         # Storing results in appropriate arrays
#         if shape_hst.error_message != "":
#         	# Error in HSM
#         	# Assigning an error code = hsm_error_code.
#         	# HST
#         	e1_hst.append(hsm_error_code)
#         	e2_hst.append(hsm_error_code)
#         	amp_hst.append(hsm_error_code)
#         	e1_hst_sh.append(hsm_error_code)
#         	e2_hst_sh.append(hsm_error_code)
#         	amp_hst_sh.append(hsm_error_code)
#         	# LSST
#         	e1_lsst.append(hsm_error_code)
#         	e2_lsst.append(hsm_error_code)
#         	amp_lsst.append(hsm_error_code)
#         	e1_lsst_sh.append(hsm_error_code)
#         	e2_lsst_sh.append(hsm_error_code)
#         	amp_lsst_sh.append(hsm_error_code)
#         	# Debug info
#         	logger.info('\tHSM Estimation Failed ...')
#         	logger.info('\t%s',shape_hst.error_message)
#         else:
# 			# HSM Converged
# 			# HST
# 			e1_hst.append(shape_hst.observed_shape.e1)
# 			e2_hst.append(shape_hst.observed_shape.e2)
# 			amp_hst.append(shape_hst.moments_amp)
# 			e1_hst_sh.append(shape_hst_sh.observed_shape.e1)
# 			e2_hst_sh.append(shape_hst_sh.observed_shape.e2)
# 			amp_hst_sh.append(shape_hst_sh.moments_amp)
# 			# LSST
# 			e1_lsst.append(shape_lsst.observed_shape.e1)
# 			e2_lsst.append(shape_lsst.observed_shape.e2)
# 			amp_lsst.append(shape_lsst.moments_amp)
# 			e1_lsst_sh.append(shape_lsst_sh.observed_shape.e1)
# 			e2_lsst_sh.append(shape_lsst_sh.observed_shape.e2)
# 			amp_lsst_sh.append(shape_lsst_sh.moments_amp)
# 			# Debug info
# 			
# 			logger.info('\tHSM Estimation Done ...\n')
# 			#logger.info(time.time())
# 
#         # Additional data to write into the data frame
#         # Ident, RA, Dec and Mag
#         cat_data = dataCatalog(cat_filename, cat_dir, k)
#         ident.append(cat_data[0])
#         ra.append(cat_data[1])
#         dec.append(cat_data[2])
#         mag.append(cat_data[3])
# 
#     # Writing image files to 4 related data cubes
#     galsim.fits.writeCube(im_hst_list, os.path.join('output', 'im_hst.fits'))
#     galsim.fits.writeCube(
#         im_hst_sh_list, os.path.join(
#             'output', 'im_hst_sh.fits'))
#     galsim.fits.writeCube(im_lsst_list, os.path.join('output', 'im_lsst.fits'))
#     galsim.fits.writeCube(
#         im_lsst_sh_list, os.path.join(
#             'output', 'im_lsst_sh.fits'))
# 
#     logger.info('\nImages Written to File...')
# 
#     # Creating a dataframe out of LST and HST data
#     data = {
#         'IDENT': ident,
#         'RA': ra,
#         'DEC': dec,
#         'MAG': mag,
#         'E1_HST': e1_hst,
#         'E2_HST': e2_hst,
#         'AMP_HST': amp_hst,
#         'E1_SH_HST': e1_hst_sh,
#         'E2_SH_HST': e2_hst_sh,
#         'AMP_HST_SH': amp_hst_sh,
#         'E1_LSST': e1_lsst,
#         'E2_LSST': e2_lsst,
#         'AMP_LSST': amp_lsst,
#         'E1_SH_LSST': e1_lsst_sh,
#         'E2_SH_LSST': e2_lsst_sh,
#         'AMP_LSST_SH': amp_lsst_sh}
#         
#     shapedata = pd.DataFrame(data, index=numpy.arange(ngal) + 1)
#     # Can also append data to the Data Frame like: "shapedata['A'] = A" where
#     # A is some column array
# 
#     # Make a Panda hdf5 data store and frame.
#     if 'output_filename' in config['HSM_Comparison']:
#         filename = os.path.join('output', config['HSM_Comparison']['output_filename'] )
#     else:
#         filename = os.path.join('output', 'shapeData.h5')
#     hdf = pd.HDFStore(filename, mode='w')
#     hdf.put('shapeData', shapedata, format='table', data_columns=True)
#     hdf.close()
# 
