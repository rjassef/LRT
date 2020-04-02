import numpy as np

from .new_SED_Model import lrt_model

class lrt_fit(lrt_model):

    def __init__(self,_jy,_ejy,_jyuse):

        #If single object is loaded, add extra dimension for compatibility.
        if len(_jy.shape)==1:
            _jy    = np.expand_array(_jy   , axis=0)
            _ejy   = np.expand_array(_ejy  , axis=0)
            _jyuse = np.expand_array(_jyuse, axis=0)

        #Initiate the main object.
        super().__init__(_jy=_jy,_ejy=_ejy,_jyuse=_jyuse)


    def kc_fit(self):
        '''Calls kca to do a full SED fit using the kca function.
           If no redshift is provided, then uses the photo-z.
        '''

        #Check data is loaded.
        if self.jy is None or self.ejy is None or self.jyuse is None:
            print("Attempted fit without loading data")
            return

        #Check some redshift is provided.
        self.z = self.zspec
        if self.z==None or self.z<0.:
            self.z = self.zspec
            if self.zphot==None:
                print("Must provide a redshift or run pz_fit.")
                return
            else:
                self.z = self.zphot

        #Check kc has been initalized. Otherwise do it.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False

        #Check if the K-correction redshift is provided. Otherwise, set 0.
        if self.z0 == None:
            self.z0 = 0.

        #Initialize needed values.
        if self.ebv is None:
            self.ebv = np.array(-1.0)
        if self.igm == None:
            self.igm = np.array(-1.0)
        if self.chi2 == None:
            self.chi2 = np.array(-1.0)

        #Set the templates to use
        lrt.ivary.ivaryobj = np.array(self.tempuse).astype(np.int32)

        #Set whether to use the prior on reddening or not
        lrt.red_igm_prior.use_red_igm_prior = self.redigmprior

        #Initialize output arrays.
        self.jymod  = np.zeros(self.nchan)
        self.jycorr = np.zeros(self.nchan)
        self.cov    = np.asfortranarray(np.zeros((self.nspec,self.nspec)))

        #Initialize the template amplitudes if not initialized yet.
        if self.comp is None:
            self.comp = np.zeros(self.nspec)

        #Must make sure jyuse is really composed of integers.
        self.jyuse = self.jyuse.astype(int)

        #Run the fit.
        lrt.kca(self.jy,self.ejy,
                self.jyuse,self.z,
                self.z0,self.jymod,
                self.jycorr,self.comp,
                self.cov,self.ebv,self.igm,
                self.chi2,self.op)
