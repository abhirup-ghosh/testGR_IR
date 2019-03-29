#!/usr/bin/python


import numpy as np
import lal
import standard_gwtransf as gw

def calc_chirptime(m1, m2, f_min=20.):
  m = float(m1+m2)
  eta = float(m1*m2/m**2)
  return (5./256.)*m*float(lal.MTSUN_SI) / ((np.pi * m * float(lal.MTSUN_SI) * f_min)**(8./3.) * eta) + 1e3*m*float(lal.MTSUN_SI)

def calc_chirptime_from_mcq(mc, q, f_min=20.):
  m = gw.tot_from_mcq(mc, q)
  eta = gw.eta_from_q(q)
  return (5./256.)*m*float(lal.MTSUN_SI) / ((np.pi * m * float(lal.MTSUN_SI) * f_min)**(8./3.) * eta) + 1e3*m*float(lal.MTSUN_SI)

if __name__=='__main__':
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt


	m = np.linspace(10, 100, 100)
	(m1, m2) = np.meshgrid(m, m[::-1], indexing='xy')
	ct = np.vectorize(calc_chirptime)(m1, m2)

	plt.figure()
	plt.imshow(2**np.ceil(np.log2(ct)), extent=(10, 100, 10, 100), cmap=plt.cm.jet)
	plt.colorbar()
	plt.xlabel('$m_1 [M_\odot]$')
	plt.ylabel('$m_2 [M_\odot]$')
	plt.xscale('log')
	plt.yscale('log')
	plt.grid()
	plt.title('seg lengths')
	plt.savefig('chirptimes.png', dpi=200)
