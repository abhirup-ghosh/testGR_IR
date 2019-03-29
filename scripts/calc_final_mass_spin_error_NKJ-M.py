"""
Estimate the errors in the fitting formulae used to compute the final mass/spin by comparing the final mass/spin estimated using the formula with NR results from the SXS catalog. 

P. Ajith, 2015-11-15 
"""

import matplotlib.pyplot as plt 
import plotsettings 
import numpy as np 
import imrtestgr as tgr 
from mpl_toolkits.mplot3d import Axes3D
import pneqns 
import os 

nr_data = '../data/nr/SXS_catalog_params.txt'
#fit_formula = 'nonprecspin_Healy2014'
fit_formula = 'nonprecspin_Husa2015'

"""
# parameters for selecting only the non-precessing spins
out_dir = '2015-12-04_non_prec_spin'
tag = 'non_prec_spin_%s' %fit_formula
MAX_IN_PLANE_SPINS = 1e-3
MAX_Q = 18.
MAX_ECCENRICITY = 1e-3 
"""

# # parameters for selecting only the precessing spins 
# out_dir = '2015-12-04_prec_spin'
# tag = 'prec_spin_incl_in-plane_spin_evol_L0ADM_%s' %fit_formula
# MIN_IN_PLANE_SPINS = 1e-3
# MAX_ECCENRICITY = 1e-3 


# parameters for selecting all quasi-circular simulations 
out_dir = '2015-12-04_allspins'
tag = 'allspins_incl_in-plane_spin_no_evol_L0ADM_%s' %fit_formula
MAX_IN_PLANE_SPINS = 1.
MAX_Q = 18.
MAX_ECCENRICITY = 1e-3 

# make output directory and copy the run script there 
os.system('mkdir -p %s' %out_dir)
os.system('cp %s %s' %(__file__, out_dir))

# read the data 
SXS_run_tag = np.loadtxt(nr_data, dtype=str, skiprows=1, usecols=(0,), unpack=True)
q, chi1, chi2, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, eccentricity, omega0, num_orbs, af, af_x, af_y, af_z, Mf, m1, m2, Jx,	Jy,	Jz,	J = np.loadtxt(nr_data, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), unpack=True)

# select a subset of simulations 
idx = np.logical_and(np.logical_and(np.logical_and(abs(chi1_x) < MAX_IN_PLANE_SPINS, abs(chi2_x) < MAX_IN_PLANE_SPINS), q < MAX_Q), eccentricity < MAX_ECCENRICITY)
#idx = np.logical_and(np.logical_or(abs(chi1_x) > MIN_IN_PLANE_SPINS, abs(chi2_x) > MIN_IN_PLANE_SPINS), eccentricity < MAX_ECCENRICITY)

q, chi1, chi2, chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z, eccentricity, omega0, num_orbs, af, af_x, af_y, af_z, Mf, m1, m2, Jx,	Jy,	Jz,	J, SXS_run_tag = q[idx], chi1[idx], chi2[idx], chi1_x[idx], chi1_y[idx], chi1_z[idx], chi2_x[idx], chi2_y[idx], chi2_z[idx], eccentricity[idx], omega0[idx], num_orbs[idx], af[idx], af_x[idx], af_y[idx], af_z[idx], Mf[idx], m1[idx], m2[idx], Jx[idx],	Jy[idx],	Jz[idx],	J[idx], SXS_run_tag[idx]

# compute the orbital ang momentum 
Lx = Jx-(chi1_x*m1**2 + chi2_x*m2**2)
Ly = Jy-(chi1_y*m1**2 + chi2_y*m2**2)
Lz = Jz-(chi1_z*m1**2 + chi2_z*m2**2)
L = (Lx**2+Ly**2+Lz**2)**0.5

# unit vector along the initial ADM ang momentum 
Lx /= L
Ly /= L
Lz /= L

# rescale all results to unit M 
M = m1+m2 
m1 /= M 
m2 /= M 
af /= Mf**2. 
Mf /= M 

"""
# compute chi_i . L at Schwarzschild ISCO 
chi1dL = np.zeros_like(q)
chi2dL = np.zeros_like(q)
Sperpmag = np.zeros_like(q) 
v_f = np.zeros_like(q)
#L0x, L0y, L0z = 0., 0., 1.
v_final = 6.**-0.5
for i in range(len(q)): 
	v0 = omega0[i]**(1./3)
	print 'i = %d v0 = %f chi1z = %f chi2z = %f. #'  %(i,  v0, chi1_z[i], chi2_z[i]),  

	# in the case of non-precessing spins chi_i . L = chi_z 
	if abs(chi1_x[i]) < 1e-3 and abs(chi2_x[i]) < 1e-3:
		chi1dL[i], chi2dL[i], Sperpmag[i], v_f[i] = chi1_z[i], chi2_z[i], 0., v_final
	# in the case of precessing simulations estimate chi_i . L at v_isco by evolving PN precession eqns 
	else: 
		#chi1dL[i], chi2dL[i], v_f[i] = pneqns.find_SdotL_at_freq(v0, m1[i]*10, m2[i]*10, chi1_x[i], chi1_y[i], chi1_z[i], chi2_x[i], chi2_y[i], chi2_z[i], L0x, L0y, L0z, v_final)
		#chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v, v_f[i] = pneqns.find_S_and_L_at_freq(v0, m1[i]*10, m2[i]*10, chi1_x[i], chi1_y[i], chi1_z[i], chi2_x[i], chi2_y[i], chi2_z[i], L0x, L0y, L0z, v_final)
		chi1x_v, chi1y_v, chi1z_v, chi2x_v, chi2y_v, chi2z_v, Lnx_v, Lny_v, Lnz_v, v_f[i] = pneqns.find_S_and_L_at_freq(v0, m1[i]*10, m2[i]*10, chi1_x[i], chi1_y[i], chi1_z[i], chi2_x[i], chi2_y[i], chi2_z[i], Lx[i], Ly[i], Lz[i], v_final)
		os.system('mv spinev.png %s/spinev_inj%d.png' %(out_dir, i))
		m1sq, m2sq = m1[i]*m1[i], m2[i]*m2[i]
		chi1dL[i] = chi1x_v*Lnx_v + chi1y_v*Lny_v + chi1z_v*Lnz_v
		chi2dL[i] = chi2x_v*Lnx_v + chi2y_v*Lny_v + chi2z_v*Lnz_v
		Sperpx = m1sq*(chi1x_v - chi1dL[i]*Lnx_v) + m2sq*(chi2x_v - chi2dL[i]*Lnx_v)
		Sperpy = m1sq*(chi1y_v - chi1dL[i]*Lny_v) + m2sq*(chi2y_v - chi2dL[i]*Lny_v)
		Sperpz = m1sq*(chi1z_v - chi1dL[i]*Lnz_v) + m2sq*(chi2z_v - chi2dL[i]*Lnz_v)
		Sperpmag[i] = (Sperpx*Sperpx + Sperpy*Sperpy + Sperpz*Sperpz)**0.5
	print 'v_f = %f: chi1dL = %f chi2dL = %f Sperpmag = %f' %(v_f[i], chi1dL[i], chi2dL[i], Sperpmag[i])
"""

# Computing scalar products

chi1dL0 = chi1_x*Lx + chi1_y*Ly + chi1_z*Lz 
chi2dL0 = chi2_x*Lx + chi2_y*Ly + chi2_z*Lz 

# compute in-plane spin vectors and add them together

s1inplanex = m1*m1*(chi1_x - chi1dL0*Lx)
s1inplaney = m1*m1*(chi1_y - chi1dL0*Ly)
s1inplanez = m1*m1*(chi1_z - chi1dL0*Lz)

s2inplanex = m2*m2*(chi2_x - chi2dL0*Lx)
s2inplaney = m2*m2*(chi2_y - chi2dL0*Ly)
s2inplanez = m2*m2*(chi2_z - chi2dL0*Lz)

#s1inplanex = m1*m1*(chi1_x - chi1dL*Lx)
#s1inplaney = m1*m1*(chi1_y - chi1dL*Ly)
#s1inplanez = m1*m1*(chi1_z - chi1dL*Lz)

#s2inplanex = m2*m2*(chi2_x - chi2dL*Lx)
#s2inplaney = m2*m2*(chi2_y - chi2dL*Ly)
#s2inplanez = m2*m2*(chi2_z - chi2dL*Lz)

sinplanemag = ((s1inplanex + s2inplanex)**2. + (s1inplaney + s2inplaney)**2. + (s1inplanez + s2inplanez)**2.)**0.5

# spin angles 
theta_s1 = np.arccos(chi1_z/chi1)
theta_s2 = np.arccos(chi2_z/chi2)

# compute the final mass and spin 
M = m1+m2 
eta = m1*m2/M**2.
Mf_fit, af_fit = tgr.calc_final_mass_spin(m1, m2, chi1dL0, chi2dL0, fit_formula)
#Mf_fit, af_fit = tgr.calc_final_mass_spin(m1, m2, chi1dL, chi2dL, fit_formula)

# Add in the (naive) contribution from in-plane spins to the fit

af_fit_corr = (af_fit*af_fit + sinplanemag*sinplanemag/Mf_fit**4.)**0.5

#af_fit_corr = (af_fit*af_fit + Sperpmag*Sperpmag/Mf_fit**4.)**0.5 

# Add the contribution from the in-plane portion of the initial ADM angular momentum to the fit

# compute the errors in the fits 
Mf_err = Mf_fit-Mf 
Mf_err_norm = Mf_fit/Mf - 1.
#af_err = abs(af_fit)-af 		# the af in the SXS catlog is the magnitude of the final spin. Hence we take the abs of the final spin given by the fit formula 
af_err = af_fit_corr - af
af_err_norm = af_fit_corr/af - 1.
#print '# RMS error in Mf = %3.2e. RMS error in af = %3.2e' %(np.std(Mf_err), np.std(af_err))
print '# RMS normalized error in Mf = %3.2e. RMS normalized error in af = %3.2e' %(np.std(Mf_err_norm), np.std(af_err_norm))

# print the simulations with large fractional errors 
#err_idx_vec = np.where(abs(af_err) > np.std(af_err))[0]
err_idx_vec = np.where(abs(af_err_norm) > np.std(af_err_norm))[0]
#err_idx_vec = np.where(abs(Mf_err) > np.std(Mf_err))[0]
for err_idx in err_idx_vec: 
	#print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] Mf = %f Mf_fit = %f af = %f af_fit = %f af_fit_corr = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], Mf[err_idx], Mf_fit[err_idx], af[err_idx], af_fit[err_idx], af_fit_corr[err_idx])
	#print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_fit_corr = %f af_err = %f af_err_corr = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], af[err_idx], af_fit[err_idx], af_fit_corr[err_idx], abs(af_fit[err_idx]) - af[err_idx], af_err[err_idx])
        print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_fit_corr = %f af_err_norm = %f af_err_norm_corr = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], af[err_idx], af_fit[err_idx], af_fit_corr[err_idx], abs(af_fit[err_idx])/af[err_idx] - 1., af_err_norm[err_idx])
	#print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] Mf = %f Mf_fit = %f Mf_err = %f af_err = %f af_err_corr = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], Mf[err_idx], Mf_fit[err_idx], Mf_err[err_idx], abs(af_fit[err_idx]) - af[err_idx], af_err[err_idx])

# Also do this for the absolute errors

print("------------------")

err_idx_vec = np.where(abs(Mf_err) > np.std(Mf_err))[0]
for err_idx in err_idx_vec: 
	print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_fit_corr = %f af_err = %f af_err_corr = %f' %(SXS_run_tag[err_idx], q[err_idx], chi1_x[err_idx], chi1_y[err_idx], chi1_z[err_idx], chi2_x[err_idx], chi2_y[err_idx], chi2_z[err_idx], af[err_idx], af_fit[err_idx], af_fit_corr[err_idx], abs(af_fit[err_idx]) - af[err_idx], af_err[err_idx])


# Print any simulations for which the new treatment increases the fractional error

#err_anom_idx_vec = np.where(abs(af_err) > abs(abs(af_fit) - af) + 1e-7)[0]
err_anom_idx_vec = np.where(abs(af_err_norm) > abs(abs(af_fit)/af - 1.) + 1e-7)[0]

if len(err_anom_idx_vec) > 0:

	print 'Larger fractional errors with new method in:'
	for err_anom_idx in err_anom_idx_vec:
		#print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_fit_corr = %f af_err = %f af_err_corr = %f' %(SXS_run_tag[err_anom_idx], q[err_anom_idx], chi1_x[err_anom_idx], chi1_y[err_anom_idx], chi1_z[err_anom_idx], chi2_x[err_anom_idx], chi2_y[err_anom_idx], chi2_z[err_anom_idx], af[err_anom_idx], af_fit[err_anom_idx], af_fit_corr[err_anom_idx], abs(af_fit[err_anom_idx]) - af[err_anom_idx], af_err[err_anom_idx])
		print 'tag = %s q = %2.1f chi1 = [%2.1f, %2.1f, %2.1f] chi2 = [%2.1f, %2.1f, %2.1f] af = %f af_fit = %f af_fit_corr = %f af_err_norm = %f af_err_norm_corr = %f' %(SXS_run_tag[err_anom_idx], q[err_anom_idx], chi1_x[err_anom_idx], chi1_y[err_anom_idx], chi1_z[err_anom_idx], chi2_x[err_anom_idx], chi2_y[err_anom_idx], chi2_z[err_anom_idx], af[err_anom_idx], af_fit[err_anom_idx], af_fit_corr[err_anom_idx], abs(af_fit[err_anom_idx])/af[err_anom_idx] - 1., af_err_norm[err_anom_idx])
	

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.hist(Mf_err*100, 20)
plt.grid()
plt.xlabel('$\Delta M_f \\times 100$')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(Mf_err), fontsize=8)
plt.subplot(122)
plt.hist(af_err*100, 20)
plt.grid()
plt.xlabel('$\Delta \chi_f \\times 100 $')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(af_err), fontsize=8)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_hist_%s.png' %(out_dir, tag), dpi=200)
plt.close()

plt.figure(figsize=(8,3.5))
plt.subplot(121)
plt.hist(Mf_err_norm*100, 20)
plt.grid()
plt.xlabel('$\Delta M_f/M_f \\times 100$')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(Mf_err_norm), fontsize=8)
plt.subplot(122)
plt.hist(af_err_norm*100, 20)
plt.grid()
plt.xlabel('$\Delta \chi_f/\chi_f \\times 100 $')
plt.ylabel('$N$')
plt.title('RMS error = %2.1e' %np.std(af_err_norm), fontsize=8)
plt.tight_layout()
plt.savefig('%s/Mf_af_err_norm_hist_%s.png' %(out_dir, tag), dpi=200)
plt.close()

#plt.figure(figsize=(4,8))
#plt.subplot(311)
#plt.semilogx(eccentricity, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
#plt.semilogx(eccentricity, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
#plt.ylabel('$M_f$')
#plt.subplot(312)
#plt.semilogx(eccentricity, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
#plt.semilogx(eccentricity, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
#plt.ylabel('$a_f/M_f$')
#plt.subplot(313)
#plt.semilogx(eccentricity, Mf_err, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f$')
#plt.semilogx(eccentricity, af_err, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f$')
#plt.xlabel('initial eccentricity')
#plt.ylabel('error')
#plt.legend(frameon=False, loc=3)
#plt.savefig('%s/Mf_af_err_vs_eccentricity_%s.png' %(out_dir, tag), dpi=200)
#plt.close()

plt.figure(figsize=(4,8))
plt.subplot(311)
plt.semilogx(eccentricity, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
plt.semilogx(eccentricity, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
plt.ylabel('$M_f$')
plt.subplot(312)
plt.semilogx(eccentricity, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
plt.semilogx(eccentricity, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
plt.ylabel('$a_f/M_f$')
plt.subplot(313)
plt.semilogx(eccentricity, Mf_err_norm, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f/M_f$')
plt.semilogx(eccentricity, af_err_norm, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f/a_f$')
plt.xlabel('initial eccentricity')
plt.ylabel('normalized error')
plt.legend(frameon=False, loc=3)
plt.savefig('%s/Mf_af_err_norm_vs_eccentricity_%s.png' %(out_dir, tag), dpi=200)
plt.close()

cos_theta_s1 =  np.cos(theta_s1)
cos_theta_s2 =  np.cos(theta_s2)

#plt.figure(figsize=(4,8))
#plt.subplot(311)
#plt.plot(cos_theta_s1, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
#plt.plot(cos_theta_s1, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
#plt.ylabel('$M_f$')
#plt.subplot(312)
#plt.plot(cos_theta_s1, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
#plt.plot(cos_theta_s1, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
#plt.ylabel('$a_f/M_f$')
#plt.subplot(313)
#plt.plot(cos_theta_s1, Mf_err, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f$')
#plt.plot(cos_theta_s1, af_err, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f$')
#plt.xlabel('$\cos \\theta_{s1}$')
#plt.ylabel('error')
#plt.legend(frameon=False, loc=4)
#plt.savefig('%s/Mf_af_err_vs_theta_s1_%s.png' %(out_dir, tag), dpi=200)
#plt.close()

plt.figure(figsize=(4,8))
plt.subplot(311)
plt.plot(cos_theta_s1, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
plt.plot(cos_theta_s1, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
plt.ylabel('$M_f$')
plt.subplot(312)
plt.plot(cos_theta_s1, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
plt.plot(cos_theta_s1, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
plt.ylabel('$a_f/M_f$')
plt.subplot(313)
plt.plot(cos_theta_s1, Mf_err_norm, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f/M_f$')
plt.plot(cos_theta_s1, af_err_norm, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f/a_f$')
plt.xlabel('$\cos \\theta_{s1}$')
plt.ylabel('normalized error')
plt.legend(frameon=False, loc=4)
plt.savefig('%s/Mf_af_err_norm_vs_theta_s1_%s.png' %(out_dir, tag), dpi=200)
plt.close()

#plt.figure(figsize=(4,8))
#plt.subplot(311)
#plt.plot(cos_theta_s2, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
#plt.plot(cos_theta_s2, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
#plt.ylabel('$M_f$')
#plt.subplot(312)
#plt.plot(cos_theta_s2, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
#plt.plot(cos_theta_s2, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
#plt.ylabel('$a_f/M_f$')
#plt.subplot(313)
#plt.plot(cos_theta_s2, Mf_err, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f$')
#plt.plot(cos_theta_s2, af_err, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f$')
#plt.xlabel('$\cos \\theta_{s2}$')
#plt.ylabel('error')
#plt.legend(frameon=False, loc=4)
#plt.savefig('%s/Mf_af_err_vs_theta_s2_%s.png' %(out_dir, tag), dpi=200)
#plt.close()

plt.figure(figsize=(4,8))
plt.subplot(311)
plt.plot(cos_theta_s2, Mf, color='c', marker='.', ms=8, ls='none', label='$M_f$')
plt.plot(cos_theta_s2, Mf_fit, color='r', marker='.', ms=4, ls='none', label='$M_f$ fit')
plt.ylabel('$M_f$')
plt.subplot(312)
plt.plot(cos_theta_s2, af, color='c', marker='.', ms=8, ls='none', label='$a_f$')
plt.plot(cos_theta_s2, af_fit, color='r', marker='.', ms=4, ls='none', label='$a_f$ fit')
plt.ylabel('$a_f/M_f$')
plt.subplot(313)
plt.plot(cos_theta_s2, Mf_err_norm, color='orange', marker='.', ms=8, ls='none', label='$\Delta M_f/M_f$')
plt.plot(cos_theta_s2, af_err_norm, color='k', marker='.', ms=4, ls='none', label='$\Delta a_f/a_f$')
plt.xlabel('$\cos \\theta_{s2}$')
plt.ylabel('normalized error')
plt.legend(frameon=False, loc=4)
plt.savefig('%s/Mf_af_err_norm_vs_theta_s2_%s.png' %(out_dir, tag), dpi=200)
plt.close()

#fig = plt.figure(figsize=(12,4))
#ax = fig.add_subplot(121, projection='3d')
#p = ax.scatter(chi1, eta, np.cos(theta_s1), s=40, c=np.log10(abs(Mf_err)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
#ax.set_xlabel('$\chi_1$')
#ax.set_ylabel('$\eta$')
#ax.set_zlabel('$\cos \\theta_{s1}$')
#ax.set_xlim(0,1)
#ax.set_ylim(0, 0.25)
#ax.set_zlim(-1,1)
#fig.colorbar(p)
#for err_idx in err_idx_vec: 
#	p = ax.scatter(chi1[err_idx]*np.ones(100), np.linspace(0, 0.25, 100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
#	p = ax.scatter(chi1[err_idx]*np.ones(100), eta[err_idx]*np.ones(100), np.linspace(-1,1,100), s=3, c='c', marker='.', lw=0, alpha=1)
#	p = ax.scatter(np.linspace(0,1,100), eta[err_idx]*np.ones(100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
#plt.title('log $|\Delta M_f|$')

#ax = fig.add_subplot(122, projection='3d')
#p = ax.scatter(chi1, eta, np.cos(theta_s1), s=40, c=np.log10(abs(af_err)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
#ax.set_xlabel('$\chi_1$')
#ax.set_ylabel('$\eta$')
#ax.set_zlabel('$\cos \\theta_{s1}$')
#ax.set_xlim(0,1)
#ax.set_ylim(0, 0.25)
#ax.set_zlim(-1,1)
#fig.colorbar(p)
#for err_idx in err_idx_vec: 
#	p = ax.scatter(chi1[err_idx]*np.ones(100), np.linspace(0, 0.25, 100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
#	p = ax.scatter(chi1[err_idx]*np.ones(100), eta[err_idx]*np.ones(100), np.linspace(-1,1,100), s=3, c='c', marker='.', lw=0, alpha=1)
#	p = ax.scatter(np.linspace(0,1,100), eta[err_idx]*np.ones(100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
#plt.title('log $|\Delta \chi_f|$')
#plt.savefig('%s/Mf_af_err_eta_chi_theta_%s.png' %(out_dir, tag), dpi=200)

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(121, projection='3d')
p = ax.scatter(chi1, eta, np.cos(theta_s1), s=40, c=np.log10(abs(Mf_err_norm)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
ax.set_xlabel('$\chi_1$')
ax.set_ylabel('$\eta$')
ax.set_zlabel('$\cos \\theta_{s1}$')
ax.set_xlim(0,1)
ax.set_ylim(0, 0.25)
ax.set_zlim(-1,1)
fig.colorbar(p)
for err_idx in err_idx_vec: 
	p = ax.scatter(chi1[err_idx]*np.ones(100), np.linspace(0, 0.25, 100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
	p = ax.scatter(chi1[err_idx]*np.ones(100), eta[err_idx]*np.ones(100), np.linspace(-1,1,100), s=3, c='c', marker='.', lw=0, alpha=1)
	p = ax.scatter(np.linspace(0,1,100), eta[err_idx]*np.ones(100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
plt.title('log $|\Delta M_f/M_f|$')

ax = fig.add_subplot(122, projection='3d')
p = ax.scatter(chi1, eta, np.cos(theta_s1), s=40, c=np.log10(abs(af_err_norm)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
ax.set_xlabel('$\chi_1$')
ax.set_ylabel('$\eta$')
ax.set_zlabel('$\cos \\theta_{s1}$')
ax.set_xlim(0,1)
ax.set_ylim(0, 0.25)
ax.set_zlim(-1,1)
fig.colorbar(p)
for err_idx in err_idx_vec: 
	p = ax.scatter(chi1[err_idx]*np.ones(100), np.linspace(0, 0.25, 100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
	p = ax.scatter(chi1[err_idx]*np.ones(100), eta[err_idx]*np.ones(100), np.linspace(-1,1,100), s=3, c='c', marker='.', lw=0, alpha=1)
	p = ax.scatter(np.linspace(0,1,100), eta[err_idx]*np.ones(100), cos_theta_s1[err_idx]*np.ones(100), s=3, c='c', marker='.', lw=0, alpha=1)
plt.title('log $|\Delta \chi_f/\chi_f|$')
plt.savefig('%s/Mf_af_err_norm_eta_chi_theta_%s.png' %(out_dir, tag), dpi=200)

#plt.close()
#plt.figure(figsize=(12,4))
#plt.subplot(121)
#plt.scatter(cos_theta_s1, cos_theta_s2, s=chi1*200, c=np.log10(abs(Mf_err)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
#plt.xlabel('$\cos \\theta_{s1}$')
#plt.ylabel('$\cos \\theta_{s2}$')
#plt.xlim(-1,1)
#plt.ylim(-1,1)
#plt.colorbar()
#plt.title('log $|\Delta {M_f}|$')

#plt.subplot(122)
#plt.scatter(cos_theta_s1, cos_theta_s2, s=chi1*200, c=np.log10(abs(af_err)), marker='.', lw=0, cmap='jet', alpha=1)
#plt.xlabel('$\cos \\theta_{s1}$')
#plt.ylabel('$\cos \\theta_{s2}$')
#plt.xlim(-1,1)
#plt.ylim(-1,1)
#plt.colorbar()
#plt.title('log $|\Delta \chi_f|$')

#plt.savefig('%s/Mf_af_err_thetas1_thetas2_%s.png' %(out_dir, tag), dpi=200)

plt.close()
plt.figure(figsize=(12,4))
plt.subplot(121)
plt.scatter(cos_theta_s1, cos_theta_s2, s=chi1*200, c=np.log10(abs(Mf_err_norm)), marker='.', lw=0, cmap='YlOrRd', alpha=1)
plt.xlabel('$\cos \\theta_{s1}$')
plt.ylabel('$\cos \\theta_{s2}$')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.colorbar()
plt.title('log $|\Delta {M_f}/M_f|$')

plt.subplot(122)
plt.scatter(cos_theta_s1, cos_theta_s2, s=chi1*200, c=np.log10(abs(af_err_norm)), marker='.', lw=0, cmap='jet', alpha=1)
plt.xlabel('$\cos \\theta_{s1}$')
plt.ylabel('$\cos \\theta_{s2}$')
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.colorbar()
plt.title('log $|\Delta \chi_f/\chi_f|$')

plt.savefig('%s/Mf_af_err_norm_thetas1_thetas2_%s.png' %(out_dir, tag), dpi=200)



