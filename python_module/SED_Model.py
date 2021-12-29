import lrt
import numpy as np

import matplotlib.pyplot as plt
from pylab import *

#######

def readobj(file,n,data):

    ''' 
    Function to read from our standard format file.
    '''

    try:
        cat = open(file)
    except IOError:
        print("Cannot open file {0}".format(file))
        return

    #First line of file has the number of targets, channels and apertures.
    (ntarg,nchan,nap) = [float(x) for x in cat.readline().split()]

    #Doubles are three per line, integers are 6 per line.
    nljy = int(math.ceil(nchan/3.))
    nluse = int(math.ceil(nchan/6.))

    #Advance to the object position in the file.
    ntot = 2+nljy*2+nluse
    for l in range(1,n):
        for k in range(ntot):
            cat.readline()

    #Read the first two lines.
    (z,nline,ra)    = [float(x) for x in cat.readline().split()]
    (dec,aux1,aux2) = [float(x) for x in cat.readline().split()]

    #Read in the fluxes, errors and use flags.
    jy = list()
    for j in range(nljy):
        jy.extend([float(x) for x in cat.readline().split()])

    ejy = list()
    for j in range(nljy):
        ejy.extend([float(x) for x in cat.readline().split()])

    jyuse = list()
    for j in range(nluse):
        jyuse.extend([int(x) for x in cat.readline().split()])

    #Set a minimum error of 5%.
    for j in range(len(jy)):
        if jyuse[j]==1 and ejy[j]<0.05*jy[j]:
            ejy[j] = 0.05*jy[j]

    #Finally, assign read values to data object.
    data.zspec  = z
    data.ra     = ra
    data.dec    = dec
    data.jy     = np.array(jy)*1e-3
    data.ejy    = np.array(ejy)*1e-3
    data.jyuse  = np.array(jyuse)

    cat.close()

######

#Function simply to tell if python is being run interactively or not.
def in_ipython():
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True

######

#Main class.
class lrt_model(object):

    #Initialize.
    def __init__(self):
        #Data
        self.zspec = None
        self.jy    = None
        self.ejy   = None
        self.jyuse = None
        self.ra    = None
        self.dec   = None
        #Model
        self.z     = None
        self.jymod = None
        self.ebv = None
        self.igm = None
        self.chi2 = None
        self.z0 = 0.
        self.jycorr = None
        self.nchan = None
        self.comp = None
        self.cov  = None
        self.op   = 0
        self.zmin = 0.0
        self.zmax = 6.0
        self.dz   = 0.01
        self.chigal = None
        self.chinop = None
        self.chi2zop = 0
        self.zphot = None
        self._vec = None
        self._ahat = None
        self.nspec = 4
        self._kcinit = False
        self._pzinit = False
        self.iverbose = 0
        self._L6um = None
        self._L5100 = None
        self._Mstar = None
        self.tempuse = [1,1,1,1]
        self.redigmprior = 1
        self.name = None
        self._abs_mag = False

    @property
    def ahat(self):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
            return
        #Calculate ahat.
        self._ahat = self.comp[0]/np.sum(self.comp)
        return self._ahat

    @property
    def L6um(self):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
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
        #Calculate L6um.
        self._L6um = lrt.get_AGN_lnu(self.comp,self.z,6.)*3.e14/6.
        return self._L6um

    @property
    def L5100(self):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
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
        #Calculate L5100. Note this is lam*L_lam
        self._L5100 = lrt.get_AGN_lnu(self.comp,self.z,0.51)*3.e14/0.51
        return self._L5100

    def L_at_lam(self, lam):
        """lam must be in microns."""
        #Check kc has been initialized.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
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
        #Calculate L_at_lam. Note this is lam*L_lam
        self._L_at_lam = lrt.get_lnu(self.comp,self.z,lam)*3.e14/lam
        return self._L_at_lam

    @property
    def Mstar(self):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
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
        #Calculate Mstar.
        self._Mstar = lrt.ml_bell03_ii(self.comp,self.z)
        return self._Mstar

    @property
    def vec(self):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
            return
        vecfac = (lrt.dl(self.z))**2 * 1e10*3e-9/(1.+self.z)
        self._vec = self.comp*lrt.alphanorm.alpha_norm/vecfac
        return self._vec

    @property
    def abs_mag(self):

        if self._abs_mag is True:
            return self._abs_mag
        
        '''
        From loaded model parameters, get the expected magnitudes.
        '''
        
        #Check model parameters are loaded.
        if self.comp is None or self.z == None or self.ebv == None or \
                self.igm == None:
            print("Attempted to create model without all parameters set.")
            return
        
        #Check kc has been initialized.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
  
        jymod_absmag = np.zeros(self.nchan)

        #Run function.
        lrt.get_mags(self.comp,self.ebv,self.igm,-1.,jymod_absmag)

        self._abs_mag = -2.5*np.log10(jymod_absmag/lrt.cal1.jyzero[:self.nchan])

        return self._abs_mag


    ###

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


    ###

    def pz_fit(self):

        '''
        Runs pza to get the photo-z.
        '''

        #Check data has been loaded.
        if self.jy is None or self.ejy is None or self.jyuse is None:
            print("Attempted fit without loading data")
            exit(1)

        #Check pz has been initialized.
        if self._pzinit == False:
            lrt.pzinit("bandmag.dat",1,1,0,self.zmin,self.zmax,self.dz,
                       self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._pzinit = True
            self._kcinit = False
        
        #Initialize output arrays/scalars.
        if self.zphot == None:
            self.zphot = np.array(-1.0) 
        if self.chigal == None:
            self.chigal = np.array(-1.0)
        if self.chinop == None:
            self.chinop = np.array(-1.0)

        #Run fit.
        lrt.pza(self.jy,self.ejy,self.jyuse,self.zphot,self.chigal,
                self.chinop,self.op,self.chi2zop)

    ###

    def get_model_fluxes(self):
        '''
        From loaded model parameters, get the expected magnitudes.
        '''
        
        #Check some redshift is provided.
        self.z = self.zspec
        if self.z==None or self.z<0.: 
            self.z = self.zspec
            if self.zphot==None:
                print("Must provide a redshift or run pz_fit.")
                return
            else:
                self.z = self.zphot

        #Check model parameters are loaded.
        if self.comp is None or self.z == None or self.ebv == None or \
                self.igm == None:
            print("Attempted to create model without all parameters set.")
            return
        
        #Check kc has been initialized.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
  
        self.jymod = np.zeros(self.nchan)

        #Run function.
        lrt.get_mags(self.comp,self.ebv,self.igm,self.z,self.jymod)

    ###

    def plot(self):

        '''
        Plots the model and the data.
        '''

        #Check kc is initialized.Need it to calculate vec from comp.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
 
        #Check fit exists.
        if self.comp is None:
            print("Need a fit before plotting.")
            return

        #Check redshift to use.
        self.z = self.zspec
        if self.z==None or self.z<0.: 
            self.z = self.zspec
            if self.zphot==None:
                print("Must provide a redshift or run pz_fit.")
                return
            else:
                self.z = self.zphot
 
        #Get the vec vectors. This are like comp, but in the natural
        #units of the templates.
        #vecfac = (lrt.dl(self.z))**2 * 1e10*3e-9/(1.+self.z)
        #self.vec = self.comp*lrt.alphanorm.alpha_norm/vecfac

        #Create the arrays with the templates in the correct scale.
        nw = lrt.wavegrid.nwave
        s = np.ndarray((4,nw))
        for l in range(self.nspec):
            s[l] = self.vec[l]*lrt.specmod1.spec[l][0:nw]
        
        #Reddening over AGN component.
        s[0] *= 10**(-0.4*lrt.dust.tau[0:nw]*self.ebv)

        #IGM absorption
        tigm = np.array([lrt.transmit(lrt.wavegrid.bcen[k],self.z,self.igm) \
                             for k in range(nw)])
        for l in range(self.nspec):
            s[l] = s[l]/lrt.wavegrid.bcen[0:nw] * tigm

        #Best fit model.
        smod = copy(s[0])
        for l in range(1,self.nspec):
            smod += s[l]
        

        lam = np.array([lrt.cal1.lbar[j]/(1.+self.z) \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])
        jy  = np.array([(1.+self.z)*self.jy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])
        ejy = np.array([(1.+self.z)*self.ejy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])

        lamu = np.array([lrt.cal1.lbar[j]/(1.+self.z) \
                             for j in range(self.nchan) if(self.jyuse[j]==2)])
        ejyu = np.array([(1.+self.z)*0.5*self.ejy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==2)])

        jym  = (1.+self.z)*self.jymod/lrt.cal1.lbar[0:self.nchan]
        lamm = lrt.cal1.lbar[0:self.nchan]/(1+self.z)

        plt.autoscale(enable=True)
        plt.clf()
        
        plt.plot(lamm,jym,'cs',markerfacecolor='None')
        plt.errorbar(lam,jy,xerr=None,yerr=ejy,fmt='ro')
        if len(lamu)>0:
            plt.errorbar(lamu,0.9*ejyu,yerr=0.1*ejyu,lolims=True,fmt=None)
        plt.plot(lrt.wavegrid.bcen[0:nw],smod,'k-')
        plt.xscale('log')
        plt.yscale('log')

        plt.autoscale(enable=False)
        plt.plot(lrt.wavegrid.bcen[0:nw],s[0],'b-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[1],'r-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[2],'g-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[3],'c-')
        plt.plot(lrt.wavegrid.bcen[0:nw],smod,'k-')

        plt.plot(lamm,jym,'cs',markerfacecolor='None')
        plt.errorbar(lam,jy,xerr=None,yerr=ejy,fmt='ro')
        if len(lamu)>0:
        #    print lamu
        #    raw_input()
            plt.errorbar(lamu,0.9*ejyu,yerr=0.1*ejyu,lolims=True,fmt=None)

        #for j in range(len(ejyu)):
        #    plt.arrow(lamu[j], ejyu[j], 0., -0.5*ejyu[j],width=0.1*lamu[j])

        plt.xlabel(r'Rest-Frame $\lambda (\mu m)$')
        plt.ylabel(r'$\propto \nu F_{\nu}$')

        ax = gca()
        xloc = 0.03
        yloc = 1.00
        dyloc = 0.05
        if self.name is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,"%s" % (self.name),
                     transform = ax.transAxes)
        if self.chi2 is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$\chi^2$'+"= %.2f" % (self.chi2),
                     transform = ax.transAxes)
        if self.zspec!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$z_{spec}$'+"= %.2f" % (self.zspec),
                     transform = ax.transAxes)
        if self.zphot!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$z_{phot}$'+"= %.2f" % (self.zphot),
                     transform = ax.transAxes)
        if self.ahat!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$\hat{a}$'+" = %.2f" % (self.ahat),
                     transform = ax.transAxes)
        if self.ebv is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$E(B-V)_{AGN}$'+" = %.2f" % (self.ebv),
                     transform = ax.transAxes)

        if not in_ipython():
            plt.show(block=True)

        
    ###

    def plot_to_file(self,fig_name):

        '''
        Plots the model and the data.
        '''

        #Check kc is initialized.Need it to calculate vec from comp.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
 
        #Check fit exists.
        if self.comp is None:
            print("Need a fit before plotting.")
            return

        #Check redshift to use.
        self.z = self.zspec
        if self.z==None or self.z<0.: 
            self.z = self.zspec
            if self.zphot==None:
                print("Must provide a redshift or run pz_fit.")
                return
            else:
                self.z = self.zphot
 
        #Get the vec vectors. This are like comp, but in the natural
        #units of the templates.
        #vecfac = (lrt.dl(self.z))**2 * 1e10*3e-9/(1.+self.z)
        #self.vec = self.comp*lrt.alphanorm.alpha_norm/vecfac

        #Create the arrays with the templates in the correct scale.
        nw = lrt.wavegrid.nwave
        s = np.ndarray((4,nw))
        for l in range(self.nspec):
            s[l] = self.vec[l]*lrt.specmod1.spec[l][0:nw]
        
        #Reddening over AGN component.
        s[0] *= 10**(-0.4*lrt.dust.tau[0:nw]*self.ebv)

        #IGM absorption
        tigm = np.array([lrt.transmit(lrt.wavegrid.bcen[k],self.z,self.igm) \
                             for k in range(nw)])
        for l in range(self.nspec):
            s[l] = s[l]/lrt.wavegrid.bcen[0:nw] * tigm

        #Best fit model.
        smod = copy(s[0])
        for l in range(1,self.nspec):
            smod += s[l]
        

        lam = np.array([lrt.cal1.lbar[j]/(1.+self.z) \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])
        jy  = np.array([(1.+self.z)*self.jy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])
        ejy = np.array([(1.+self.z)*self.ejy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==1)])

        lamu = np.array([lrt.cal1.lbar[j]/(1.+self.z) \
                             for j in range(self.nchan) if(self.jyuse[j]==2)])
        ejyu = np.array([(1.+self.z)*0.5*self.ejy[j]/lrt.cal1.lbar[j] \
                             for j in range(self.nchan) if(self.jyuse[j]==2)])

        jym  = (1.+self.z)*self.jymod/lrt.cal1.lbar[0:self.nchan]
        lamm = lrt.cal1.lbar[0:self.nchan]/(1+self.z)

        #Do not display the plot.
        plt.ioff()

        plt.autoscale(enable=True)
        plt.clf()
        
        plt.plot(lamm,jym,'cs',markerfacecolor='None')
        plt.errorbar(lam,jy,xerr=None,yerr=ejy,fmt='ro')
        if len(lamu)>0:
            plt.errorbar(lamu,0.9*ejyu,yerr=0.1*ejyu,lolims=True,fmt=None)
        plt.plot(lrt.wavegrid.bcen[0:nw],smod,'k-')
        plt.xscale('log')
        plt.yscale('log')

        plt.autoscale(enable=False)
        plt.plot(lrt.wavegrid.bcen[0:nw],s[0],'b-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[1],'r-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[2],'g-')
        plt.plot(lrt.wavegrid.bcen[0:nw],s[3],'c-')
        plt.plot(lrt.wavegrid.bcen[0:nw],smod,'k-')

        plt.plot(lamm,jym,'cs',markerfacecolor='None')
        plt.errorbar(lam,jy,xerr=None,yerr=ejy,fmt='ro')
        if len(lamu)>0:
        #    print lamu
        #    raw_input()
            plt.errorbar(lamu,0.9*ejyu,yerr=0.1*ejyu,lolims=True,fmt=None)

        #for j in range(len(ejyu)):
        #    plt.arrow(lamu[j], ejyu[j], 0., -0.5*ejyu[j],width=0.1*lamu[j])

        plt.xlabel(r'Rest-Frame $\lambda (\mu m)$')
        plt.ylabel(r'$\propto \nu F_{\nu}$')

        ax = gca()
        xloc = 0.03
        yloc = 1.00
        dyloc = 0.05
        if self.name is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,"%s" % (self.name),
                     transform = ax.transAxes)
        if self.chi2 is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$\chi^2$'+"= %.2f" % (self.chi2),
                     transform = ax.transAxes)
        if self.zspec!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$z_{spec}$'+"= %.2f" % (self.zspec),
                     transform = ax.transAxes)
        if self.zphot!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$z_{phot}$'+"= %.2f" % (self.zphot),
                     transform = ax.transAxes)
        if self.ahat!=None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$\hat{a}$'+" = %.2f" % (self.ahat),
                     transform = ax.transAxes)
        if self.ebv is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$E(B-V)_{AGN}$'+" = %.2f" % (self.ebv),
                     transform = ax.transAxes)

        plt.savefig(fig_name)

        
