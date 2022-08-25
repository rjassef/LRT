import numpy as np
import matplotlib.pyplot as plt

import lrt

#The LRT object class
class lrt_model(object):

    #Initialize
    def __init__(self,jy=None,ejy=None,jyuse=None,z=None):

        #Values needed to consider that the object is declared.
        self.jy    = jy
        self.ejy   = ejy
        self.jyuse = jyuse
        self.z     = z

        self.comp = None
        self.zphot = None

        return

    @property
    def ahat(self):
        if self.comp is None:
            return None
        _ahat = np.zeros(self.nobj)
        for nobj in range(self.nobj):
            _ahat[nobj] = self.comp[nobj,0]/np.sum(self.comp[nobj])
        return _ahat

    @property
    def L6um(self,zuse=None):
        if self.comp is None:
            return None
        if zuse is None:
            if self.zbset is not None:
                zuse = self.zbest
            return None
        _L6um = np.zeros(self.nobj)
        for nobj in range(self.nobj):
            _L6um[nobj] = lrt.get_lnu(self.comp[nobj],zuse[nobj],6.)*3.e14/6.
        return _L6um

    ######

    @property
    def nobj(self):
        if not self.ready_to_fit:
            return None
        return jy.shape[0]

    #Prefer the provided redshifts. If any are <0, then use the photo-z if available. If neither are available, give an error message.
    @property
    def zbest(self):

        #First check if we will need to check the photo-zs.
        if len(self.z[self.z>0])==len(self.z):
            return self.z

        #If any are 0 or negative, we'll need photo-zs. Check if they exist.
        if self.zphot is None:
            return None

        #If they exist, replace the non-positive zs with zphots.
        zbest = np.where(self.z>0,self.z,self.zphot)
        return zbest
