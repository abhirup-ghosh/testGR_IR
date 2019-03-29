#!/usr/bin/env python

"""
Script to check injected dirstibution in uniform in comoving time-volume
"""

import numpy as np
from scipy import integrate
import standard_cosmology as sc
import matplotlib.pyplot as plt

injtxtfile = '../injections/example_inj_SEOBNRv2_ROM_DoubleSpinthreePointFivePN_uniform_comoving_volume.txt'
injdata = np.genfromtxt(injtxtfile, dtype=None, names=True)
d_max = np.max(injdata['distance'])

comov_time_vol_func = lambda d: sc.volume_dHoc(d*sc.H0/sc.c)/(1.+sc.redshift(d*sc.H0/sc.c))
(normalization, n_err) = integrate.quad(comov_time_vol_func, 0., d_max)

d_arr = np.linspace(0., d_max, 100)
vol_arr = np.vectorize(comov_time_vol_func)(d_arr)/normalization

plt.figure()
plt.hist(injdata['distance'], bins=10, normed=True, histtype='stepfilled', color='c', alpha=0.2, label='Injection distances')
plt.plot(d_arr, vol_arr, 'b-', label='Uniform in comov VT')
plt.xlabel('Distance (Mpc)')
plt.legend(loc='lower right')
plt.show()
