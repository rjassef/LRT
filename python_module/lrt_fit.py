import numpy as np

from .new_SED_Model import lrt_model

class lrt_fit(lrt_model):

    #Initialize
    def __init__(self,jy,ejy,jyuse=None,z=None):

        if len(shape(jy))==1:
            jy  = np.expand_dims(jy,axis=0)
            ejy = np.expand_dims(ejy,axis=0)

        if jyuse is None:
            jyuse = np.ones(jy.shape, dtype=np.int32)
        else:
            self.jyuse = jyuse
            #Must make sure jyuse is really composed of integers.
            self.jyuse = self.jyuse.astype(int)

        if z is None:
            self.z = np.ones(jy.shape[0])*-1
        else:
            self.z = z

        self.ebv  = -np.ones(jy.shape[0])
        self.igm  = -np.ones(jy.shape[0])
        self.chi2 = -np.ones(jy.shape[0])

        self.redigmprior = 1

        #Initiate the main object.
        super().__init__(jy=jy,ejy=ejy,jyuse=jyuse,z=z)
        return


    def kc_fit(self, zuse=None, z0=0.):
        '''Calls kca to do a full SED fit using the kca function.
           If no redshift is provided, then uses the photo-z.

           zuse = Redshift estimate. If not provided will use zbest.

           z0   = Redshift to which one wants to K-correct to.
        '''

        #Check redshifts are avilable.
        if zuse is None:
            if self.zbest is not None:
                zuse = self.zbest
            else:
                return

        #Check kc has been initalized. Otherwise do it.
        try:
            self._kcinit
        except NameError:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self.nspec = lrt.specmod1.nspec
            self._kcinit = True

        #Set the templates to use
        lrt.ivary.ivaryobj = np.array(self.tempuse).astype(np.int32)

        #Set whether to use the prior on reddening or not
        lrt.red_igm_prior.use_red_igm_prior = self.redigmprior

        #Initialize output arrays.
        self.jymod  = np.zeros(self.jy.shape)
        self.jycorr = np.zeros(self.jy.shape)
        self.cov    = np.asfortranarray(np.zeros((self.nobj,self.nspec,self.nspec)))
        self.comp = np.zeros((self.nobj,self.nspec))

        #For speed, we should later on implement a multi-core approach.
        for k in self.nobj:
            jy    = self.jy[nchan*k:nchan*(k+1)]
            ejy   = self.ejy[nchan*k:nchan*(k+1)]
            jyuse = self.jyuse[nchan*k:nchan*(k+1)]
            z     = zuse[k]

            jymod  = np.zeros(self.nchan)
            jycorr = np.zeros(self.nchan)
            comp   = np.zeros(self.nspec)
            cov    = np.zeros((self.nspec,self.nspec))
            ebv    = np.zeros(1)
            igm    = np.zeros(1)
            chi2   = np.zeros(1)

            #Run the fit.
            lrt.kca(jy, ejy, jyuse, z, z0, jymod, jycorr, comp, cov, ebv, self.igm, self.chi2, 0)

            self.jymod[k]  = jymod
            self.jycorr[k] = jycorr
            self.comp[k]   = comp
            self.cov[k]    = cov
            self.ebv[k]    = ebv[0]
            self.igm[k]    = igm[0]
            self.chi2[k]   = chi2[0] 
