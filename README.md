## The HSMShapeData module in GalSim:

The **HSMShapeData** module is defined to be used as an external module as a part of the YAML based configuration in GalSim. The module provides the user with a framework to generate shear data from multiple sources and store the resulting information in a HDF file.

The implementation can be summarized as follows:

* Defining an independent module, written in Python, that can be imported into the YAML configuration file.
* Read the image data at source, run HSM over it and store the shear data as a Pandas data-frame into a HDF5 file, which can be used as an input to the analysis code.
* An independent analysis code, one of whose functions is to compare any 2 variables out of the above descibed data-frame.


## Usage:

* Import the module into the configuration (YAML) file (using the 'modules' field) and define other generic parameters (using fields such as 'input', 'image', 'stamp', etc.) that are common throughout the configuration file.

* Define a list of galaxies and PSFs that will be referenced later to contruct images.

* Construct the images by referencing the list of galaxies and PSFs defined before and generate shear data from them using the *HSMShapeData* module.  

For example, the following is a code snippet showing the module being used as an extra output type in GalSim to generate shear data from a particular combination of the galaxy and the PSF.

```
 HSMShapeData: 
 		# Output directory
        dir: output_yaml 
        # Output filename
        file_name: shapeData.h5 
        # Need to specify the order of entries in the 'psf' and 'gal' list 
        defined above. The names in the list would also be used as tags for 
        the keys in the output data-frame.
        # The tags that are valid currently are:
        # Galaxies:
        # 'COSMOS': for COSMOS galaxies, 'Gaussian': for Gaussian Galaxies
        # An extra "_SH" is added if the galaxy is sheared
        # PSFs:
        # 'HST' for HST PSF, 'LSST' for LSST PSF
        # This would probably be changed later to something that is more 
        flexible/versatile.
        gal_tags: ['COSMOS', 'COSMOS_SH', 'Gaussian']
        psf_tags: ['HST', 'LSST']


```