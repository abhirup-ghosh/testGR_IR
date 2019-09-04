"""
calculate the cutoff frequencies for the inspiral-merger-ringdown consistency test 

Abhirup Ghosh, P. Ajith, 2015-09-17 
"""


import lal, lalinference, numpy as np 
import lalinference.imrtgr.imrtgrutils as tgr
import nr_fits as nr

# set binary parameters here -- for the time being we neglect the spin paraemter 
# use the final-mass/final-spin formula for non-spinning binaries 
m1, m2, chi1, chi2, tilt1, tilt2, phi12 = 181.27334813487664, 51.50050954651401, 0.9613367288207788, 0.9198680010606237, 2.5895286659491017, 2.730092017643872, 3.2589835114322456
chi1z, chi2z = chi1*np.cos(tilt1), chi2*np.cos(tilt2)

#fit_formula = 'nospin_Pan2011'
#fit_formula = 'nonprecspin_Healy2014'
fit_formula = 'bbh_average_fits_precessing'

# calculate the mass and spin of the final BH 
Mf, af = tgr.calc_final_mass_spin(m1, m2, chi1, chi2, chi1z, chi2z, phi12, fit_formula)
Mf, af = Mf[0], af[0]

# calculate the Kerr ISCO freq 
f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

print 'Mf = %3.2f af = %3.2f f_isco_Kerr = %2.3f Hz ' %(Mf, af, f_isco_Kerr)
