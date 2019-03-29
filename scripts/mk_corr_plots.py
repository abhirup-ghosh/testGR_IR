"""
Plot the correlation between different parameters from the LALInference runs 

P. Ajith, 2015-11-03
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tgrplotsettings 
import os 
import imrtestgr as tgr

base_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-04/uniform_mtot_q_prior'
run_vec = ['IMR', 'inspiral', 'post-inspiral']
inj_vec = [2,4,8,10,12,13,14,16,18,19,20,22,25,26,28,33,34,36,39,43,46,49]
inj_vec = [2]
fit_formula = 'nonprecspin_Healy2014'
#fit_formula = 'nospin_Pan2011'

out_dir = '/home/ajith/working/cbc/testGR_IR/results/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-04/uniform_mtot_q_prior/joint_posteriors_%s' %fit_formula
os.system('mkdir -p %s' %out_dir)

for inj in inj_vec: 
	print '# processing injection ', inj
	plt.figure(figsize=(32,10))
	for ir, run in enumerate(run_vec): 
	
		print '... run ', run
		post_file = '%s/injection_%d/%s/1126285216-0/H1L1/posterior_samples.dat' %(base_dir, inj, run)

		# get injection params 
		dL_inj = float(os.popen("ligolw_print -c distance %s/injection_%d/%s/injection.xml" %(base_dir, inj, run)).read().splitlines()[0])
		m1_inj = float(os.popen("ligolw_print -c mass1 %s/injection_%d/%s/injection.xml" %(base_dir, inj, run)).read().splitlines()[0])
		m2_inj = float(os.popen("ligolw_print -c mass2 %s/injection_%d/%s/injection.xml" %(base_dir, inj, run)).read().splitlines()[0])

		# read the data 
		data = np.genfromtxt(post_file, dtype=None, names=True)

		# posterior samples in various parameters 
		m1, m2, chi1, chi2, dL, logL = data['m1'], data['m2'], data['a1'], data['a2'], data['dist'], data['logl']
		m = m1+m2 
		eta = m1*m2/m**2.
		Mf, af = tgr.calc_final_mass_spin(m1, m2, chi1, chi2, fit_formula)
		chi = (m1*chi1+m2*chi2)/m
		logL = chi

		# injected parameters 
		M_inj = m1_inj+m2_inj
		chi1_inj = 0.
		chi2_inj = 0.
		chi_inj = 0.
		eta_inj = m1_inj*m2_inj/M_inj**2.
		Mf_inj, af_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)

		# plot the scatter plot of the samples colred by logL 
		plt.subplot(3,8,ir*8+1)
		plt.scatter(m, dL, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$M~[M_\odot]$')
		plt.ylabel('$d_L~$[Mpc]')
		plt.ylim(1, 2e3)
		plt.xlim(0.5*M_inj, 2*M_inj)
		plt.axvline(x=M_inj, color='k', lw=0.5)
		plt.axhline(y=dL_inj, color='k', lw=0.5)

		plt.subplot(3,8,ir*8+2)
		plt.scatter(eta, dL, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$\eta$')
		plt.ylabel('$d_L~$[Mpc]')
		plt.ylim(1, 2e3)
		plt.xlim(0., 0.25)
		plt.axvline(x=eta_inj, color='k', lw=0.5)
		plt.axhline(y=dL_inj, color='k', lw=0.5)

		plt.subplot(3,8,ir*8+3)
		plt.scatter(m, eta, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$M~[M_\odot]$')
		plt.ylabel('$\eta$')
		plt.ylim(0, 0.25)
		plt.xlim(0.5*M_inj, 2*M_inj)
		plt.axvline(x=M_inj, color='k', lw=0.5)
		plt.axhline(y=eta_inj, color='k', lw=0.5)
		plt.title(run)

		plt.subplot(3,8,ir*8+4)
		plt.scatter(Mf, af, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$M_f~[M_\odot]$')
		plt.ylabel('$a_f/M_f$')
		plt.ylim(0,1)
		plt.xlim(0.5*Mf_inj, 2*Mf_inj)
		plt.grid()
		plt.axvline(x=Mf_inj, color='k', lw=0.5)
		plt.axhline(y=af_inj, color='k', lw=0.5)

		plt.subplot(3,8,ir*8+5)
		plt.scatter(Mf/m, af, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$M_f/M$')
		plt.ylabel('$a_f/M_f$')
		plt.ylim(0,1)
		plt.xlim(0.9,1)
		plt.grid()
		plt.axvline(x=Mf_inj/M_inj, color='k', lw=0.5)
		plt.axhline(y=af_inj, color='k', lw=0.5)

		plt.subplot(3,8,ir*8+6)
		plt.scatter(m1, m2, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$m_1~[M_\odot]$')
		plt.ylabel('$m_2~[M_\odot]$')
		plt.axvline(x=m1_inj, color='k', lw=0.5)
		plt.axhline(y=m2_inj, color='k', lw=0.5)

		plt.subplot(3,8,ir*8+7)
		plt.scatter(chi1, chi2, s=2, c=logL, lw=0)
		plt.colorbar()
		plt.xlabel('$\chi_1$')
		plt.ylabel('$\chi_2$')
		plt.xlim(-1,1)
		plt.ylim(-1,1)
		plt.axvline(x=chi1_inj, color='k', lw=0.5)
		plt.axhline(y=chi2_inj, color='k', lw=0.5)
		
		plt.subplot(3,8,ir*8+8)
		plt.scatter(chi, eta, s=2, c=logL, lw=0)
		plt.xlim(-1,1)
		plt.colorbar()
		plt.xlabel('$\chi$')
		plt.ylabel('$\eta$')
		plt.axvline(x=chi_inj, color='k', lw=0.5)
		plt.axhline(y=eta_inj, color='k', lw=0.5)


	#plt.savefig('%s/corr_plots_%s_inj%d.png' %(out_dir, fit_formula, inj), dpi=200)
	plt.savefig('%s/corr_plots_%s_inj%d_colorchi.png' %(out_dir, fit_formula, inj), dpi=200)
	plt.close()

