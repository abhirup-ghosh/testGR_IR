import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import read_emcee_samples as res
import pycbc
from pycbc  import  detector

postloc_lalinf_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalinference/IMR/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'
postloc_lalinf_insp = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalinference/inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'
postloc_lalinf_ring = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalinference/post-inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'

postloc_tdimr_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/IMR/emcee_samples.dat'
postloc_tdimr_insp = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/inspiral_132Hz_m0p004185/emcee_samples.dat'
postloc_tdimr_ring = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/emcee_samples.dat'

postloc_fdimr_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_FD_reruns_3629/IMR/emcee_samples.dat'
postloc_fdimr_insp = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_FD_reruns_3629/inspiral_132Hz_m0p004185/emcee_samples.dat'
postloc_fdimr_ring = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_FD_reruns_3629/post-inspiral_132Hz_m0p004185/emcee_samples.dat'

data = np.genfromtxt(postloc_lalinf_ring, dtype=None, names=True)
mc1, q1, dL1, i1, t01, phi01, ra1, sin_dec1, pol1 = data['mc'], data['q'], data['distance'], data['theta_jn'], data['time'], data['phase'], data['ra'], np.sin(data['dec']), data['psi']
Fp1,Fc1 = detector.overhead_antenna_pattern(ra1, np.arcsin(sin_dec1), pol1)

nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc2, q2, dL2, i2, t02, phi02, ra2, sin_dec2, pol2 = res.read_emcee_samples_9dim(postloc_tdimr_ring, nwalkers, num_iter, ndim, n_burnin)
Fp2,Fc2 = detector.overhead_antenna_pattern(ra2, np.arcsin(sin_dec2), pol2)

nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc3, q3, dL3, i3, t03, phi03, ra3, sin_dec3, pol3 = res.read_emcee_samples_9dim(postloc_fdimr_ring, nwalkers, num_iter, ndim, n_burnin)
Fp3,Fc3 = detector.overhead_antenna_pattern(ra3, np.arcsin(sin_dec3), pol3)

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)
#ax7 = fig.add_subplot(337)
#ax8 = fig.add_subplot(338)
ax1.hist(mc1, bins=50, normed=True, histtype='step', label='LALInf')
ax1.hist(mc2, bins=50, normed=True, histtype='step', label='Emcee: TD')
ax1.hist(mc3, bins=50, normed=True, histtype='step', label='Emcee: FD')
ax1.legend(loc='upper left')
ax1.set_xlabel('Mc')
ax2.hist(q1, bins=50, normed=True, histtype='step')
ax2.hist(q2, bins=50, normed=True, histtype='step')
ax2.hist(q3, bins=50, normed=True, histtype='step')
ax2.set_xlabel('q')
ax3.hist(dL1, bins=50, normed=True, histtype='step')
ax3.hist(dL2, bins=50, normed=True, histtype='step')
ax3.hist(dL3, bins=50, normed=True, histtype='step')
ax3.set_xlabel('distance')
ax4.hist(i1, bins=50, normed=True, histtype='step')
ax4.hist(i2, bins=50, normed=True, histtype='step')
ax4.hist(i3, bins=50, normed=True, histtype='step')
ax4.set_xlabel('iota')
ax5.hist(Fp1, bins=50, normed=True, histtype='step')
ax5.hist(Fp2, bins=50, normed=True, histtype='step')
ax5.hist(Fp3, bins=50, normed=True, histtype='step')
ax5.set_xlabel('Fp')
ax6.hist(Fc1, bins=50, normed=True, histtype='step')
ax6.hist(Fc2, bins=50, normed=True, histtype='step')
ax6.hist(Fc3, bins=50, normed=True, histtype='step')
ax6.set_xlabel('Fc')
#ax7.hist(t01-1126285216.0, bins=50, normed=True, histtype='step')
#ax7.hist(t02, bins=50, normed=True, histtype='step')
#ax7.set_xlabel('time')
#ax8.hist(phi01, bins=50, normed=True, histtype='step')
#ax8.hist(phi02, bins=50, normed=True, histtype='step')
#ax8.set_xlabel('phase')
plt.tight_layout()
plt.savefig('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20180718_new_runs/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/1D_comparison_ring.png')
