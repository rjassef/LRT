#Fortran 77 Compiler
FC = gfortran
#FC = g77
#FC = ifort

#Compiler Flags
ifeq ("$(shell uname -m)", "arm64")
FCFLAGS = -ffixed-line-length-none -std=legacy #-fbounds-check
else
FCFLAGS = -ffixed-line-length-none -std=legacy -mcmodel=medium #-fbounds-check
endif
#FCFLAGS = -ffixed-line-length-none
#FCFLAGS = -extend_source -O #-check

#AR and Ranlib
AR = /usr/bin/ar
RANLIB = /usr/bin/ranlib

#Removal Tool
RM = /bin/rm -f

#Object Files
LIB_OBJS = kc.o pz.o common_lib.o kca_util.o util.o min.o fitzero.o heapsort.o star.o

#Library Name
LFILE = liblrt.a

all: $(LFILE)

$(LFILE): $(LIB_OBJS)
	$(AR) rc $(LFILE) $(LIB_OBJS)
	$(RANLIB) $(LFILE)

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

clean: 
	$(RM) $(LIB_OBJS)

distclean: 
	$(RM) $(LIB_OBJS) $(LFILE)
