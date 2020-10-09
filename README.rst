Index
-----
  0.- Installation
  I.- How to Calculate Photometric Redshifts/K-Corrections
 II.- Setting the Luminosity Priors for Photometric Redshifts
III.- Setting the Cosmology
 IV.- Summary of Available Functions
  V.- Description of Available Functions

===================================================================

0.- Installation
    ------------

For installation instructions, see INSTALL file.

In a nutshell, first define the environment variable LRTPATH to the
folder where the files will be installed. Then run

make
make install

It is recommended that you run the test codes in the test folder and
compare the outputs to the zphot.bak and kcorr.bak files.


===================================================================


I.- How to Calculate Photometric Redshifts/K-Corrections
    ----------------------------------------------------

The LRT libraries use the low resolution spectral templates of Assef
et al. (2007) to calculate useful quantities like photometric
redshifts, K-corrections and bolometric luminosities.

To install the libraries, please refer to the INSTALL file included
with the package. 

REMEMBER TO SET LRTPATH BEFORE RUNNING THE CODE.

The LRT libraries provide functions for your Fortran-77 codes. To use
the functions, a few steps have to be followed. In here we show how to
construct a basic code that will calculate photometric redshifts and
K-corrections.

The order is the following: 1) Create the photometry description
file. 2) Initialize the libraries from your code. 3) Call the
functions for doing the respective calculations.

1) Create the Photometry Description File
   --------------------------------------

First, put the filters you will use in the Filters folder were the
code have been installed. This file has two columns: the first one is
wavelength in Angstroms and the second one is the response. Lines
starting with a # will be considered comments. THE EXTENSION OF THE
FILE MUST BE .FILTER .

Second, create the photometry description file. This file has three
columns. The first is the name of the filter file (without the .filter
termination). The second is the normalization flag: 1 if Vega
magnitude/fluxes are being used, 2 if IRAC and 3 if AB. The third and
last column is the mean flux of your calibration source in Jy through
the corresponding filter. A list of Vega mean fluxes can be found in
http://nsted.ipac.caltech.edu/NStED/docs/parhelp/Photometry.html and
AB ones are always 3631 Jy.

An example of this file can be found in the bandmag.dat file included
with the test files. Feel free to modify and use all the test files.

See webpage for a list of included filters.




2) Initialize the Libraries from your Code
   ---------------------------------------

Before calling the library functions, the initialization functions
must be called. For photometric redshifts, 

	call pzinit(filtname,nspec,zmin,zmax,dz,verbose)

where

	filtname: string	The name of the photometry
                  		description file ('bandmag.dat' 
				in the test) 

	nspec   : integer       The number of templates to use 
		  		(3 and 4 are supported. See Assef 
				et al. 2007  for details)

	zmin    : real*8 	Minimum redshift to search

	zmax	: real*8	Maximum redshift to search

	dz	: real*8	Redshift grid separation

	verbose : integer 	If 1, information will be print to the
				STDOUT. If 0, nothing will be printed
				except for possible errors.

Notice that the photoz_grid.dat file will be created. This file is for
internal use of the pza function. For K-corrections and all other
functions,

	call kcinit(filtname,nspec,verbose)

where filtname, nspec and verbose are the same as above.




3) Call the functions for doing the respective calculations.
   --------------------------------------------------------

For photometric redshifts, having declared 'nchan' filters: 

	call pza(mag,emag,maguse,zphot,chigal,chinop,op,chi2zop)

where

	mag(nchan): real*8 	Input array of size nchan with the 
				Magnitudes or Fluxes (in Jy) in each
				band arranged in the same order as in
				the photometry description file.

	emag(nchan): real*8 	Input array of size nchan with the 
				magnitude or fluxes (in Jy) errors in 
				each band

	maguse(nchan):integer 	Input array of size nchan with the use
				flags. If maguse[j] = 1, the band will
				be used for estimating photometric
				redshifts. If maguse[j] = 0, it will not.
				Still, a modeled magnitude or flux
				will be	returned for this filter.

	zphot	: real*8	Output. Photometric redshift output 
				variable

	chigal	: real*8	Output. Chi2 + Prior of the best fit

	chinop	: real*8	Output. Chi2 of the best fit

	op 	: integer 	Input. 0 if input is in Fluxes (in Jy) 
				and 1 if in magnitudes.

	chi2zop : integer 	Input. If >0, the chi-squared distribution 
				of the object for which a photometric
				redshift is calculated will be written
				to the fort.90 file. The first line
				has two columns, the chi2zop flag
				value (should be the galaxy ID for
				easier recognition) and the number of
				lines in the redshift grid, N. The
				following N lines show the
				distribution and have three columns:
				redshift, chi2, chi2+prior. Following
				objects on the same run are appended
				to the file.


For K-Corrections, 

	call kca(mag,emag,maguse,z,z0,magmod,magcorr,comp,op)

where

	mag(nchan): real*8 	Input array of size nchan with the 
				Magnitudes or Fluxes (in Jy) in each
				band arranged in the same order as in
				the photometry description file.

	emag(nchan): real*8 	Input array of size nchan with the 
				magnitude or fluxes (in Jy) errors in 
				each band

	maguse(nchan):integer 	Input array of size nchan with the use
				flags. If maguse[j] = 1, the band will
				be used for estimating photometric
				redshifts. If maguse[j] = 0, it will not.
				Still, a modeled magnitude or flux
				will be	returned for this filter.

	z	: real*8	Input Redshift of the object

	z0	: real*8	Input Redshift to which K-correct. 
				Typically 0.d0

	magmod(nchan): real*8 	Output array of size nchan on which the
				best fit model fluxes or magnitudes
				are returned

	magcorr(nchan):real*8	Output array of size nchan on which the 
				K-corrections will be returned. 
				mag - magcorr = K-corrected magnitude. 
				Flux/Flux_corr= K-corrected flux.

	comp(nspec): real*8	Output array of size nspec (3 or 4) on 
				which the best fit specific
				luminosities of each component are
				returned.

	op 	: integer 	Input. 0 if input is in Fluxes (in Jy) 
				and 1 if in magnitudes.	

=======================================================================


II.- Setting the Luminosity Priors for Photometric Redshifts
     -------------------------------------------------------

	The LRT libraries can use luminosity function priors for
estimating photometric redshifts. As discussed in Assef et al. (2007),
the priors improve the estimates by 5 - 10% in general.

	The priors are set by the "pzinit" function which in turn
calls the "setlumprior" internal function. The latter will look for a
file called 'prior.dat' in the folder where you are running your
program. This file must have 1 column and 4 lines with 1 number each
(everything after the first number in a line is considered
comments). Lines starting with a '#' are considered comments.

The first number indicates whether to use or not the luminosity
priors. If 1, the priors will be used and if 0 they will not.

The second number indicates the band to which this will be applied. The
number must correspond to the photometry description file location of
the filter from top to bottom.

The third is the value of M*.

The fourth is the value of alpha.


If no file is provided, the function will not use a luminosity prior.

=====================================================================

III.- Setting the Cosmology
      ---------------------

	This is very similar to setting the luminosity priors. The
cosmology is set by the "pzinit" and/or "kcinit" function, which in
turns call the "setdist" internal function. The latter will search for
a file called cosmo.dat in the folder where your program is run. This
file must have 1 column and 4 lines with 1 number each. Lines starting
with a '#' are considered comments.

The first number is Omega Matter.
The second is Omega_Lambda.
The third is Omega_k.
The fourth is H0 in units of km/s/Mpc .

If no file is provided a cosmology of (Omega_Matter, Omega_Lambda,
Omega_k, H0) = (0.3, 0.7, 0.0, 70.0) is assumed.

=====================================================================


IV.- Summary of Available Functions
     ------------------------------


The LRT main functions are pza and kca, which calculate photometric
redshifts and K corrections respectively, but we also provided other
potentially useful functions. A summary in alphabetical order follows:


bol_lum	: 	Subroutine that returns the bolometric luminosity of a
		galaxy given the specific luminosities returned by kca.

get_mags: 	Subroutine that given specific luminosities for each 
		component and a redshift, returns the associated
		magnitudes in the set of filters of the Photometry
		Description File. This Subroutine can be used to
		create mock galaxies.

DL 	: 	Function that returns the luminosity distance at a certain
		redshift

fitzero : 	Subroutine that estimates corrections to the input
		photometric zero points. This can significantly improve
		photometric redshifts and K-corrections.

kca	:	Subroutine to calculate K-Corrections.

pza	:	Subroutine to calculate Photometric Redshifts

kcinit 	:	Subroutine to initialize the K-corrections (and most other)
		routines.

pzinit	: 	Subroutine to initialize the photometric redshifts routines.

vc	:	Function that returns the co-moving volume at redshift z.

vmax 	:	Function that returns the maximum redshift at to which a
		galaxy could be found given a limiting magnitude.

==========================================================================

V.- Description of Available Functions
    ----------------------------------

For pza and kca see section I-3. For pzinit and kcinit see section I-2.


- real*8 function bol_lum(comp)
  -----------------------------

	Parameters:	

	  comp(nspec):	real*8		INPUT. Specific luminosities 
					of each template
					component. This vector is an
					output of kca.

	Description: 	Given the output vector of KCA with the specific
			bolometric luminosities, returns the
			bolometric luminosity of the galaxy in units
			of 10^10 solar luminosities. In practice, it
			just adds the components of vector comp. 

	Requisites:	call kcinit
			call kca


- subroutine get_mags(comp,z,jymodtot)
  ------------------------------------

	Parameters: 	

	  comp(nspec):	real*8		INPUT. Specific luminosities 
					of each template
					component. This vector is an
					output of kca.

	  z:		real*8		INPUT. Redshift of the object.

	  jymodtot(nchan): real*8	OUTPUT. Vector that holds the 
					expected fluxes in each filter.

	Description: 	Subroutine that given specific luminosities for each 
			component and a redshift, returns the associated
			fluxes in the set of filters of the Photometry
			Description File. This Subroutine can be used to
			create mock galaxies.

	Requisites:	call kcinit or pzinit
	


- real*8 function DL(z)
  ---------------------

	Parameters:	

	  z:		real*8		INPUT. Redshift of the object.

	Description:	Returns the luminosity distance to redshift z.

	Requisites:	call pzinit or call kcinit


- subroutine fitzero(filename, niter, chifrac, corr, op)
  ------------------------------------------------

	Parameters:
	  
	  filename: 	character*100 	INPUT. Name of the training set
					file. See below for the
					description of the file
					format.

	  niter:	integer		INPUT. Number of iterations for 
					estimating the
					corrections. Convergence is
					somewhat slow, but
					10 or 15 iterations are
					usually good enough.

	  chifrac:	real*8		INPUT. Fraction of objects to use, 
					arranged by chi-squared.

	  corr(nchan):	real*8		OUTPUT. Array of size nchan with 
					the zero point corrections,
					such that zeropoint*corr is
					the corrected zero point.

	  op :		 integer 	INPUT. 0 if input is in fluxes 
					(in Jy) and 1 if in magnitudes.

	Description:	Subroutine that estimates corrections to the input
			photometric zero points using a user-supplied
			training set. This corrections can
			significantly improve photometric redshifts and
			K-corrections. The format of the training set
			file must be redshift, flux (or magnitude) in
			each band, flux (or magnitude) error in each
			band, and flags for each band (1 if its to be
			used and 0 if not). Lines starting with # are
			considered comments. An example file (not
			usable but with the correct format) is
			provided in the test folder as zero.cat . Once
			the corrections have been estimated, they are
			automatically incorporated to the pza and kca
			functions.

	Requisites:	call pzinit or kcinit

- real*8 function vc(z)
  ---------------------

	Parameters:	

	  z:		real*8		INPUT. Redshift of the object.

	Description:	Returns the co-moving volume from redshift 0 
			to redshift z

	Requisites:	call pzinit or kcinit


- real*8 function vmax(comp,jchan,mlim,zlim,area) 	
  -----------------------------------------------

	Parameters: 	

	  comp(nspec):	real*8		INPUT. Specific luminosities 
					of each template
					component. This vector is an
					output of kca.

	  jchan:	integer		INPUT. Band in which the 
					magnitude limits applies.

	  mlim:		real*8		INPUT. Faint apparent magnitude 
					limit 

	  zlim:		real*8 		OUTPUT. Maximum redshift to which 
					galaxy can be found.

	  area:		real*8 		INPUT. Area of the survey in square 
					degrees.


	Description: 	Function that returns the co-moving volume on
			which a galaxy could have been detected given
			the magnitude limit of the survey and the area
			in square degrees. This function will also
			return the maximum redshift to which the
			galaxy could have been observed.
