from platform import mac_ver
from Star_Model import StarModel
import numpy as np

data = np.loadtxt("fake_stellar_catalog.dat",skiprows=1)

for k in range(data.shape[0]):

    ms_star = StarModel(st="MS", bandmag_file="bandmag_stars.dat")
    ms_star.jy = data[k,:]
    ms_star.ejy = 0.1*np.ones(len(ms_star.jy))
    ms_star.jyuse = np.ones(len(ms_star.jy),dtype=np.int32)

    ms_star.fit()
    ms_star.plot()
