"""
Find the maximum total mass corresponding to a given mass ratio, spin and f_lower 
such that fCut is greater than f_lower. 

P. Ajith, 2015-02-06 
"""

import lalsimulation as lalsim 
import lal 

# fix the injection parameters 
m1_inj = 10. 	# mass 1 (M_sun)  
m2_inj = 20. 	# mass 2 (M_sun) 
chi1_inj = 0.8 
chi2_inj 
chi_inj = chis_s_inj + delta*chi_a	# effective spin parameter 

f_lower = 561.	# low-frequency cutoff (Hz) 

# compute total mass and symmetric mass ratio 
mtot_inj = m1_inj+m2_inj
eta_inj = m1_inj*m2_inj/mtot_inj**2. 

# compute the fCut corresponding to this injection 
fcut_inj = lalsim.SimIMRPhenomBGetFinalFreq(m1_inj, m2_inj, chi_inj)

# max value of m_tot such that fcut is >= f_lower 
mtot_max = fcut_inj*mtot_inj/f_lower

print 'fcut = %2.1f mtot_max = %2.1f' %(fcut_inj, mtot_max) 
