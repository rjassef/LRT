Low Resolution Templates K Corrections and Photometric Redshifts
----------------------------------------------------------------

Copyrights: R.J. Assef & C.S. Kochanek, 2009

These libraries make use of the low resolution spectral templates of
Assef et al., 2009, ApJ, submitted . Please consult the paper for
information about the templates or the algorithms.

PLEASE REPORT ANY BUGS TO: rjassef[at]astronomy.ohio-state.edu

If you have any questions or suggestions please send them to the
address above.

This file describes the installation and testing of the libraries. The
README file describes the subroutines of the libraries.


======================================================================

INDEX
-----

  I.- REQUIREMENTS
 II.- COMPILATION AND INSTALLATION
III.- TESTING THE LIBRARIES
IV .- PYTHON MODULE

======================================================================


I.- REQUIREMENTS
    ------------

  1.- A Fortran 77 compiler. By default, the GNU compiler g77 is
      used. To use a different compiler, edit the file src/Makefile
      and modify the variables FC anf FCFLAGS accordingly.

  2.- LAPACK must be installed. Notice that if using LAPACK version
      3.2 or later, the libraries cannot be compiled with
      g77. Instead, you must use gfortran (or the compiler used to
      compile LAPACK if built from source).



I.- COMPILATION AND INSTALLATION
    ----------------------------

We have provided two methods for compiling and installing the
libraries. The first one (prefered) uses an evironmental variable
called LRTPATH pointing to an installation directory so you can call
the libraries without draging files around. The second one omits the
variable but you must copy certain files to every folder were you will
compile your codes.

After compiling anf installing, it is recomended you run the test
programs (see TESTING LIBRARIES section) to assure everything went ok.


1) Prefered Method

The programs provided here are libraries for your fortran 77 programs,
not stand alone programs. A template program that will work as stand
alone has also been included, but keep in mind that these algorithms
are intended as libraries.

A Makefile is provided to compile and install the libraries. First
define the environmental variable LRTPATH to point to the installation
directory, for example:

	setenv LRTPATH '/usr/local/lrt'       -----> tcsh

	LRTPATH '/usr/local/lrt'              -----> bash
	export LRTPATH


LRTPATH IS NECESARY TO RUN THE CODE, NOT JUST FOR INSTALLATION. You
should add the previous lines to your .tcshrc or .bashrc . If you wish
to use the codes without setting the LRTPATH, go to the end of this
section.

Once you have set LRTPATH, type

make 
make install

in the main folder. Remember that the environmental variable LRTPATH
should be set everytime the program is run. If the program is to be
installed to the folder where the code has been unpacked (most likely
where this INSTALL file is located and from where you will be runing
the Makefile), set the LRTPATH to such folder and omit "make install"
from above.

2) Alternative Method

Alternatevely, you can omit setting the LRTPATH. Just type make in the
main folder and then copy to the folder where your program will be
running the file 'liblrt.a' and the folders 'Filters' and 'specs'.


=========================================================================

II.- TESTING THE LIBRARIES
     ------- --- ---------

Two sample programs have been included to test that the code has been
properly installed and everything is working. One calculates
photometric redshifts for a sample of galaxies and the other estimates
K-Corrections.

Feel free to modify these codes.

These instructions assume that you used the prefered installation
method. If you didn't, copy the file 'liblrt.a' and the folders
'Filters' and 'specs' to the test folder, and change '$LRTPATH' for
'.' everywhere.


*Photometric Redshifts Test
 --------------------------

To compile and run the photometric redshifts test, first install the
libraries as described in the installation section. 

Next, go to the directory where you unpacked the original tar file
(most likely the folder where this INSTALL file is located) and enter
the test directory. Compile the code and run it:

	g77 -o zphot zphot.f -L$LRTPATH -llrt -llapack
	./zphot

Notice that this is the way you will compile every code. In case you
want to use another (more efficient) compiler, you have to recompile
and install the libraries using such compiler. For this, modify the
Makefile in src folder.

The output should be in a file called zphot.dat and should be equal to
the zphot.bak file. When running this code, the file photoz_grid.dat
will be created. This file is internal to the pzinit and pza
functions. The fort.90 file will also be created, containing the full
chi-squared distribution of each object.

*K-Corrections Test
 ------------------

To compile and run, type: 

	g77 -o kcorr kcorr.f -L$LRTPATH -llrt -llapack
	./kcorr

The output will be on the file kcorr.dat and should be equal to the
kcorr.bak file.

================================================================================

IV.- PYTHON MODULE
     -------------

     The python module is an f2py compiled version of the LRT
     libraries such that they can be called through python.

     For using the module see the instructions in the python_module
     folder.

     To install it, use commands

     make lrtpy
     make installpy

     in the root folder. Additional to creating the LRTPATH
     environment variable, you will have to add it to the PYTHONPATH
     environment variable you use. This would look like

     setenv PYTHONPATH $PYTHONPATH":"$LRTPATH   -----> tcsh

     PYTHONPATH=$PYTHONPATH:$LRTPATH            -----> bash
     export PYTHONPATH

     To test the installation worked properly, go to the test folder
     and type

     python test_python_module

**NOTE**

     After upgrading to astroconda from Ureka, this stopped
     working. The problem seems to be with the version of ld being
     used by Centos6, and nothing due to astroconda itself. To solver
     it, I installed a newer version of binutils from
     devtoolset-6. Specifically I did: 

     sudo yum install centos-release-scl

     sudo yum install devtoolset-6

     scl enable devtoolset-6 bash

     and then compiled the module. In different Linux distributions
     upgrading ld might work differently.

