from platform import mac_ver
from Star_Model import StarModel
import numpy as np

data = np.loadtxt("fake_stellar_catalog.dat",skiprows=1)

for k in range(data.shape[0]):

    jy = data[k,:]
    ejy = 0.1*np.ones(len(jy))
    jyuse = np.ones(len(jy),dtype=np.int32)

    chimin = None
    stype_best = None
    for stype in ["MS","GS","SGS"]:
        star = StarModel(st=stype, bandmag_file="bandmag_stars.dat")
        star.jy = jy
        star.ejy = ejy
        star.jyuse = jyuse
        star.fit()
        if chimin is None or star.chi2<chimin:
            chimin = star.chi2
            stype_best = star.st

    star = StarModel(st=stype_best, bandmag_file="bandmag_stars.dat")
    star.jy = jy
    star.ejy = ejy
    star.jyuse = jyuse
    star.fit()
    #star.plot()
    print("{0} {1:.2f} {2:.3e}".format(star.st, star.tfit, star.chi2))
