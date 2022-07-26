import lrt
import numpy as np

class StarModel(object):

    def __init__(self, st="MS"):

        #Initialize the stellar fits.
        star_type = {
            "MS":1,
            "GS":2,
            "SGS":3,
        } 
        lrt.starinit('bandmag.dat',star_type[st],0)
        self.star_type = star_type

        #Initiliaze some parameters needed to pass to the fortran functions.
        self.comp = np.zeros(2)
        self.ns_best = None
        self.chi2 = None
        self.jymod = None

        return

    def fit(self):
        #Make sure that the fluxes have been loaded.
        try:
            lrt.star_fit(self.jy,self.ejy,self.jyuse,self.jymod,self.comp,self.ns_best,self.chi2,0)
        except AttributeError:
            print("Could not run lrt.star_fit. Make sure all fluxes have been set before calling.")

        return
