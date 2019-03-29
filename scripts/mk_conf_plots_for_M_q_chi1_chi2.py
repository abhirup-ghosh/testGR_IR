import lal, numpy as np
import matplotlib.pyplot as plt
import os

spin_fit_formula_list = ['nonprecspin_Healy2014']
injection_no_list = list(np.linspace(1,50,50))

in_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-18/uniform_mtot_q_prior_spins'

M_inj_list = []
q_inj_list = []
eta_inj_list = []
chi1_inj_list = []
chi2_inj_list = []
gr_conf_level_list = []

for injection_no in injection_no_list:
	post_loc = '%s/imrtestgr_nonprecspin_Healy2014/injection_%d'%(in_dir_root, injection_no)
	if os.path.isdir(post_loc) == True:
		m1_inj, m2_inj, chi1_inj, chi2_inj = np.loadtxt('%s/injection_%d/IMR/injection.txt'%(in_dir_root,injection_no), usecols=(1,0,6,7))
		M_inj = m1_inj + m2_inj
        	q_inj = m2_inj/m1_inj
        	eta_inj = (m1_inj + m2_inj)/M_inj**2.
		gr_conf_level = np.loadtxt('%s/data/GR_confidence.txt'%post_loc)

		M_inj_list = np.append(M_inj_list, M_inj)
		q_inj_list = np.append(q_inj_list, q_inj)
		eta_inj_list = np.append(eta_inj_list, eta_inj)
		chi1_inj_list = np.append(chi1_inj_list, chi1_inj)
		chi2_inj_list = np.append(chi2_inj_list, chi2_inj)
		gr_conf_level_list = np.append(gr_conf_level_list, gr_conf_level)

plt.figure(figsize=(16,5))
plt.subplot(121)
plt.scatter(M_inj_list, q_inj_list, c=gr_conf_level_list, s=50)
plt.colorbar()
plt.subplot(122)
plt.scatter(chi1_inj_list, chi2_inj_list, c=gr_conf_level_list, s=50)
plt.colorbar()
plt.savefig('%s/imrtestgr_nonprecspin_Healy2014/conf_region_plots_for_M_q_chi1_chi2.png'%in_dir_root)
