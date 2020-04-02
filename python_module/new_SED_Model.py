import numpy as np
import matplotlib.pyplot as plt

import lrt

#The LRT object class
class lrt_model(object):

    #Initialize
    def __init__(self,_jy=None,_ejy=None,_jyuse=None,_z=None):

        #Values needed to consider that the object is declared.
        self.jy    = _jy
        self.ejy   = _ejy
        self.jyuse = _jyuse

        return

    @property
    def ahat(self):
        try:
            if self.comp is None:
                return None
        except NameError:
            return None
        _ahat = np.zeros(self.nobj)
        for nobj in range(self.nobj):
            _ahat[nobj] = self.comp[nobj,0]/np.sum(self.comp[nobj])
        return _ahat

    @property
    def L6um(self):
        try:
            if self.z is None or self.comp is None:
                return None
        except NameError:
            return None
        _L6um = np.zeros(self.nobj)
        for nobj in range(self.nobj):
            _L6um[nobj] = lrt.get_lnu(self.comp[nobj],self.z[nobj],6.)*3.e14/6.
        return _L6um

    ######

    @property
    def ready_to_fit(self):
        if self.jy is None or self.ejy is None or self.jyuse is None:
            print("Need to declare object.jy, object.ejy and object.jyuse in order to run anything.")
            return False
        return True

    @property
    def nobj(self):
        if not self.ready_to_fit:
            return None
        return jy.shape[0]

    @property
    def z(self):
        _z = np.zeros(self.nobj)
        for nobj in range(nobj):

        #Always use zspec if available and positive.
        try:
            _z = self.zspec
        except NameError:
            _z = -1.0*np.ones(self.nobj)

        if _z<0. or _z is None:
            try:
                _z = self.zphot
            except NameError:
                _z = None

        if _z is None:
            print("Must provide a redshift or run pz_fit.")

        return _z
