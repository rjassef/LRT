#!/usr/bin/env python

#This code is intended for testing the python module. It will grab the
#first object in the photometry file, calculate the photometric
#redshift, then do the SED modeling and make a plot of the best-fit
#SED with the photometry. 

#Remember that you need to have numpy and matpoltlib installed. Ureka
#should have all you need to run these codes.

import numpy as np
import SED_Model

#First, lets read the photometry of the first object from the file.
cat = open("gal_phot.dat")
ntarg = -1
nchan = -1
while True:
    line = cat.readline()
    if line[0]=="#":
        continue
    ntarg,nchan = [int(float(ix)) for ix in line.split()]
    break
z = float(cat.readline().split()[0])
jy = list()
for i in range(2,6):
    jy.extend([float(ix) for ix in cat.readline().split()])
jy = np.array(jy)
ejy = list()
for i in range(6,10):
    ejy.extend([float(ix) for ix in cat.readline().split()])
ejy = np.array(ejy)
jyuse = list()
for i in range(10,12):
    jyuse.extend([int(float(ix)) for ix in cat.readline().split()])
jyuse = np.array(jyuse).astype(np.int32)



#Now, lets create the object
gal = SED_Model.lrt_model()
#gal.iverbose = 1
#Lets fill up the info
gal.zspec = z
gal.jy = jy
gal.ejy = ejy
gal.jyuse = jyuse

#Get the photo z
gal.pz_fit()
print("Photo z = ",gal.zphot)
print("Should be 0.44")
print 

#Get the SED model
gal.kc_fit()

#Get the stellar mass.
print("Mstar/M_sun = {0:.3e}".format(gal.Mstar))
print("Should be 8.990e+10")

#Plot
print("In the plot you should see the chi_squared value of the fit.")
print("If everything went well, you should have obtained 1.61")
gal.plot()
#gal.plot_to_file("first_galaxy.eps")


#Get the best-fit model fluxes for a fake galaxy only composed from
#the E galaxy of the previous object.
#gal2 = SED_Model.lrt_model()
#gal2.zspec = gal.zspec
#gal2.comp = [0.,gal.comp[1],0.,0.]
#gal2.ebv = gal.ebv
#gal2.igm = gal.igm
#gal2.get_model_fluxes()
#print gal2.jymod



