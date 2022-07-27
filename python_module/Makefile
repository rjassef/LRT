#Fortran 77 Compiler
FC = gfortran

#Compiler Flags. Use -fbounds-check only for debugging.
FC_FLAGS = --fixed-line-length-none -std=legacy
#FC_FLAGS = '--fixed-line-length-none -std=legacy'
#FC_FLAGS='--fixed-line-length-none -std=legacy -fbounds-check'

#Removal Tool
RM = /bin/rm -f

#Set here the lapack and blas path if not standard.
LAPACKPATH = ''
#LAPACKPATH='/usr/lib64/atlas/'

#Source Files
LIB_SOURCE = ../src/kc.f ../src/pz.f ../src/common_lib.f ../src/kca_util.f ../src/util.f ../src/min.f ../src/fitzero.f ../src/heapsort.f ../src/star.f 

#Set here the location of libpython2.7 if needed.
#LIBPYPATH = "/home/rjassef/Software/Ureka/Ureka/python/lib/python2.7/config/"
#LIBPYPATH = "/home/rjassef/anaconda2/envs/iraf27/lib/python2.7/config/"
#LIBPYPATH = "/home/rjassef/anaconda2/envs/iraf27/lib/python2.7/config/"
LIBPYPATH = ''

#Name of the shared lrt library
LFILE = lrt

all: $(LFILE).so

$(LFILE).so: $(LIB_SOURCE)
	f2py --f77exec=$(FC) --f77flags="$(FC_FLAGS)" -L$(LIBPYPATH) -L$(LAPACKPATH) -llapack -lblas -m $(LFILE) -c $(LIB_SOURCE)

clean:
	$(RM) $(LFILE).so
