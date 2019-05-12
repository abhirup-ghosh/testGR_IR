"""
calculate the cutoff frequencies for the inspiral-merger-ringdown consistency test 

Abhirup Ghosh, 2019-04-25 
"""

import lalinference.imrtgr.imrtgrutils as tgr
import sys
sys.path.insert(0,'../src')
import nr_fits as nr
import lal, numpy as np

data = np.genfromtxt('./SXS_campaign.dat', names=True, dtype=None, skip_footer=2)

fit_formula = 'bbh_average_fits_precessing'

for idx in range(len(data)):

  m1, m2 = data['m1'][idx], data['m2'][idx]
  s1x, s1y, s1z = data['s1x'][idx], data['s1y'][idx], data['s1z'][idx]
  s2x, s2y, s2z = data['s2x'][idx], data['s2y'][idx], data['s2z'][idx]

  if s1x == 0.:
    s1x = 1e-6

  if s1y == 0.:
    s1y = 1e-6

  if s1z == 0.:
    s1z = 1e-6

  if s2x == 0.:
    s2x = 1e-6
  
  if s2y == 0.:
    s2y = 1e-6

  if s2z == 0.:
    s2z = 1e-6

  s1 = np.sqrt(s1x**2. + s1y**2. + s1z**2.)
  s2 = np.sqrt(s2x**2. + s2y**2. + s2z**2.)
  phi12 = np.arctan(s2y/s2x) - np.arctan(s1y/s1x)

  # calculate the mass and spin of the final BH 
  Mf, af = tgr.calc_final_mass_spin(m1, m2, s1, s2, s1z, s2z, phi12, fit_formula)

  # calculate the Kerr ISCO freq 
  f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

  print m1, m2, s1, s2, Mf[0], af[0], data['name'][idx], f_isco_Kerr[0]
