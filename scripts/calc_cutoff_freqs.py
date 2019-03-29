"""
calculate the cutoff frequencies for the inspiral-merger-ringdown consistency test 

P. Ajith, 2015-09-17 
"""


import imrtgrutils_final as tgr
import nr_fits as nr 
import lal, numpy as np 


# set binary parameters here -- for the time being we neglect the spin paraemter 
# use the final-mass/final-spin formula for non-spinning binaries 
m1, m2, chi1, chi2, tilt1, tilt2, phi12 = 50., 50., 0., 0., 0., 0., 0.
chi1z, chi2z = chi1*np.cos(tilt1), chi2*np.cos(tilt2)

#fit_formula = 'nospin_Pan2011'
#fit_formula = 'nonprecspin_Healy2014'
fit_formula = 'bbh_average_fits_precessing'

# calculate the mass and spin of the final BH 
Mf, af = tgr.calc_final_mass_spin(m1, m2, chi1, chi2, chi1z, chi2z, tilt1, tilt2, phi12, fit_formula)

# calculate the Kerr ISCO freq 
f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

# calculate the dominant QNM freq 
f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)

print 'Mf = %3.2f af = %3.2f f_isco_Kerr = %2.1f Hz f_qnm = %2.1f Hz ' %(Mf, af, f_isco_Kerr, f_qnm)
