import lrt
import numpy as np
import matplotlib.pyplot as plt
from pylab import gca

class StarModel(object):

    def __init__(self, st="MS", bandmag_file="bandmag.dat"):

        #Initialize the stellar fits.
        star_type = {
            "MS":1,
            "GS":2,
            "SGS":3,
        } 
        lrt.starinit(bandmag_file,star_type[st],0)
        self.star_type = star_type
        self.nchan = lrt.data1b.nchan
        self.st = st

        #Initiliaze some parameters needed to pass to the fortran functions.
        self.comp = np.zeros(2)
        self.ns_best = np.array(-1)
        self.chi2 = np.array(-1.0)
        self.jymod = None
        self.op = 0

        self.nspec = 2
        self.z = 0.

        return

    def fit(self):
        #Make sure that the fluxes have been loaded.
        try:
            self.jymod = np.zeros(len(self.jy))
            lrt.stf(self.jy,self.ejy,self.jyuse,self.jymod,self.comp,self.ns_best,self.chi2,self.op)
            frac = self.comp[1]/np.sum(self.comp)
            self.tfit = self.ns_best+frac
        except AttributeError as e:
            print(e)
            print("Could not run lrt.star_fit. Make sure all fluxes have been set before calling.")

        return

    def plot(self, figname=None):

        '''
        Plots the model and the data.
        '''

        #Create the arrays with the templates in the correct scale.
        nw = lrt.wavegrid.nwave
        s = np.ndarray((self.nspec,nw))
        for l in range(self.nspec):
            s[l] = self.comp[l]*lrt.specmod_stars.spec_stars[l+self.ns_best-1][0:nw]/lrt.wavegrid.bcen[0:nw]

        #Best fit model.
        smod = np.copy(s[0])
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
        plt.plot(lrt.wavegrid.bcen[0:nw],smod,'k-')

        plt.plot(lamm,jym,'cs',markerfacecolor='None')
        plt.errorbar(lam,jy,xerr=None,yerr=ejy,fmt='ro')
        if len(lamu)>0:
            plt.errorbar(lamu,0.9*ejyu,yerr=0.1*ejyu,lolims=True,fmt=None)

        plt.xlabel(r'Rest-Frame $\lambda (\mu m)$')
        plt.ylabel(r'$\propto \nu F_{\nu}$')

        ax = gca()
        xloc = 0.03
        yloc = 1.00
        dyloc = 0.05
        if hasattr(self,"name") and self.name is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,"%s" % (self.name),
                     transform = ax.transAxes)
        if self.chi2 is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'$\chi^2$'+"= %.2f" % (self.chi2),
                     transform = ax.transAxes)
        if self.st is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,r'Stellar Type: {}'.format(self.st),
                     transform = ax.transAxes)
        if self.ns_best is not None:
            yloc -= dyloc
            plt.text(xloc, yloc,'Templates: {0:.2f})'.format(self.tfit),
                     transform = ax.transAxes)

        if figname is None:
            plt.show()
        else:        
            plt.savefig(figname)

        return
    
