""" 
Plot the evidence of LALInference runs with different upper cutoff frequencies. We are trying to see if the full IMR signal has a larger evidence as compared to an inspiral-merger only signal. 

P. Ajith, 2015-12-09
"""
import numpy as np 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import tgrplotsettings 
import nr_fits as nr
import imrtestgr as tgr
import lal
import waveforms as wave

# median values of the parmaeter estimates of GW150914
# from this run: https://dogmatix.icts.res.in/~abhirup/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-11-18_fhigh_insp134Hz_flow_ring_134Hz_psdfit_margspline/lalinf_imr/posplots.html
m1, m2, a1z, a2z = 39.99, 30.60, 0.27, -0.58

# calculate the mass and spin of the final BH, and the QNM frequency of the least damped mode 
Mf, af = tgr.calc_final_mass_spin(m1, m2, a1z, a2z, fit_formula='nonprecspin_Healy2014')

# complex QNM frequency of the l=m=2 mode (n=0,1,2 overtones) 
Omega0 = nr.qnmfreqs_berti(af, 2, 2, 0)	
Omega1 = nr.qnmfreqs_berti(af, 2, 2, 1)	
Omega2 = nr.qnmfreqs_berti(af, 2, 2, 2)	

# real part of the QNM freqs (Hz) 
f_qnm = np.real(Omega0)/(2*np.pi*Mf*lal.MTSUN_SI)
f_qnm_n1 = np.real(Omega1)/(2*np.pi*Mf*lal.MTSUN_SI)
f_qnm_n2 = np.real(Omega2)/(2*np.pi*Mf*lal.MTSUN_SI)

post_basedir = '/home/abhirup/public_html/ER8/burst_trigger/SEOBNRv2_ROM_DS_thread_safe_lalinference_dl2prior/no_psdfit_marg/'
insp_upp_freq_vec = [1024, 50, 100, 150, 175, 200, 225, 250, 275, 300, 325]
chain = 0 
data = np.zeros((4, 4))

plt.figure(figsize=(5,4))
for insp_upp_freq in insp_upp_freq_vec: 

	# read the data from all parallel chains 
	for chain in range(4): 
		if insp_upp_freq < 350.: 
			file = '%s/inspiral_%dHz/engine/lalinferencenest-0-H1L1-1126259462.39-%d.dat_B.txt' %(post_basedir, insp_upp_freq, chain)
		else: 
			file = '%s/IMR/engine/lalinferencenest-0-H1L1-1126259462.39-%d.dat_B.txt' %(post_basedir, chain)
		data[chain,:] = np.loadtxt(file, unpack=True) 
	
	# bayes factor, signal evidence, noise evidence, and max_L values averaged over all chains 
	bayes_fac, signal_evid, noise_evid, max_L = np.mean(data, axis=0)
	print insp_upp_freq, bayes_fac
		
	if insp_upp_freq > 1000.: 
		bayes_fac0, signal_evid0, noise_evid0, max_L0 = bayes_fac, signal_evid, noise_evid, max_L

	plt.plot(insp_upp_freq, bayes_fac, 'r.', ms=10)
	print insp_upp_freq, bayes_fac

plt.xlabel('$f_\mathrm{upper-cutoff}$ [Hz]')
plt.ylabel('$\log~B^\mathrm{signal}_\mathrm{noise}$')
plt.ylim(0.65*bayes_fac0, 1.*bayes_fac0)
plt.xlim(100, 350)
plt.axvline(x=f_qnm, color='orange', lw=1)
plt.axvline(x=f_qnm_n1, color='b', lw=0.5)
plt.axvline(x=f_qnm_n2, color='m', lw=0.25)

# plot the F(t) of the best fit waveform 
t, hp, hc, A, F = wave.generate_waveforms(m1, m2, a1z, a2z, 20, approx='SEOBNRv2')

max_amp_idx = np.where(A == np.max(A))
F_maxamp = F[max_amp_idx]
plt.axvline(x=F_maxamp, color='c', lw=1)
plt.text(f_qnm+5, 0.9*bayes_fac0, '%2.1f Hz $(f_\mathrm{qnm}$ of $\ell = m = 2, n = 0)$' %f_qnm, color='orange', fontsize=7, rotation=90)
plt.text(f_qnm_n1-6, 0.85*bayes_fac0, '%2.1f Hz $(f_\mathrm{qnm}$ of $\ell = m = 2, n = 1)$' %f_qnm_n1, color='b', fontsize=7, rotation=90)
plt.text(f_qnm_n2-10, 0.82*bayes_fac0, '%2.1f Hz $(f_\mathrm{qnm}$ of $\ell = m = 2, n = 2)$' %f_qnm_n2, color='magenta', fontsize=7, rotation=90)
plt.text(F_maxamp+5, 0.95*bayes_fac0, '%2.1f Hz $(F_\mathrm{max-ampl})$ ' %F_maxamp, color='c', fontsize=8, rotation=90)

plt.tight_layout()
plt.savefig('bases_factor_accumulation.png', dpi=300)
plt.close()

		
