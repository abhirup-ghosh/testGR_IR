import matplotlib.pyplot as plt 
import plotsettings 
import numpy as np 
import imrtgrutils as tgr



m1_targ, m2_targ = 24.0, 25.0
data_dir = 'IMRTGR_WaveformSystStudy_2017-08-29/SEOBNRv2_vs_IMRPhenomPv2/aLIGOZeroDetHighPower'
fname = 'Match_SEOBNRv2_vs_IMRPhenomPv2_aLIGOZeroDetHighPower_m1%.2f_m2%.2f.dat' %(m1_targ, m2_targ)

data_dir = 'IMRTGR_WaveformSystStudy_2017-08-29/NRPNHybridtxt_vs_IMRPhenomPv2/aLIGOZeroDetHighPower'
fname = 'Match_NRPNHybridtxt_vs_IMRPhenomPv2_aLIGOZeroDetHighPower_m127.91_m232.09.txt'

data_dir = 'IMRTGR_WaveformSystStudy_2017-08-29/NRPNHybridtxt_vs_IMRPhenomPv2_post-inspiral/aLIGOZeroDetHighPower/'
fname = 'Match_NRPNHybridtxt_vs_IMRPhenomPv2_post-inspiral_aLIGOZeroDetHighPower_m127.91_m232.09.txt'

data_dir = 'IMRTGR_WaveformSystStudy_2017-08-29/NRPNHybridtxt_vs_IMRPhenomPv2_inspiral/aLIGOZeroDetHighPower/'
fname = 'Match_NRPNHybridtxt_vs_IMRPhenomPv2_inspiral_aLIGOZeroDetHighPower_m127.91_m232.09.dat'

data_dir = 'IMRTGR_WaveformSystStudy_2017-09-11/NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_post-inspiral/aLIGOZeroDetHighPower'
fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_post-inspiral_aLIGOZeroDetHighPower_m127.91_m232.09.txt'
fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_post-inspiral_aLIGOZeroDetHighPower_m1_27.91_m2_32.09_spintempl.dat'
tag = 'post-inspiral_spin'

#data_dir = 'IMRTGR_WaveformSystStudy_2017-09-11/NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_inspiral/aLIGOZeroDetHighPower'
#fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_inspiral_aLIGOZeroDetHighPower_m127.91_m232.09.txt'
#fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_inspiral_aLIGOZeroDetHighPower_m1_27.91_m2_32.09_spintempl.dat'
#tag = 'inspiral_spin'

tag_vec = ['inspiral', 'post-inspiral']
N = 2000 
Mf_data = np.zeros((2, N))
af_data = np.zeros((2, N))

for i_d, tag in enumerate(tag_vec): 

	data_dir = 'IMRTGR_WaveformSystStudy_2017-09-11/NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_%s/aLIGOZeroDetHighPower' %tag 
	fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_%s_aLIGOZeroDetHighPower_m127.91_m232.09.txt' %tag
	fname = 'Match_NRPNHybridtxt_vs_SEOBNRv2_ROM_DoubleSpin_%s_aLIGOZeroDetHighPower_m1_27.91_m2_32.09_spintempl.txt' %tag 
	m1, m2, s1z, s2z, match = np.loadtxt('%s/%s'%(data_dir, fname), unpack=True, usecols=(0,1,4,7,10))
	Mf, af = tgr.calc_final_mass_spin(m1, m2, 0., 0., s1z, s2z, 0., 0., 0., fit_formula='nonprecspin_Healy2014')
	chi_eff = (m1*s1z+m2*s2z)/(m1+m2)

	M = m1+m2 
	eta = m1*m2/M**2 

	m_targ = m1_targ+m2_targ
	eta_targ = m1_targ*m2_targ/m_targ**2 
	q_targ = 1.15 
	m_targ = 60. 
	eta_targ = q_targ/(1+q_targ)**2 
	m1_targ = m_targ/(1.+ q_targ) 
	m2_targ = q_targ*m_targ/(1.+q_targ) 

	mf_targ, af_targ = tgr.calc_final_mass_spin(m1_targ, m2_targ, 0., 0., 0, 0, 0., 0., 0., fit_formula='nonprecspin_Healy2014')
	s1z_targ, s2z_targ = 0., 0.
	chi_eff_targ = 0. 

	max_idx = np.where(match == np.max(match))

	plt.figure(figsize=(10,2.5))
	plt.subplot(141)
	plt.scatter(M, eta, s=1*match**8, c=np.log10(1-match))
	plt.plot(m_targ, eta_targ, 'k<')
	plt.plot(M[max_idx], eta[max_idx], 'rx')
	plt.xlim(min(M), max(M))
	plt.ylim(min(eta), max(eta))
	plt.xlabel('$M~[M_\odot]$')
	plt.ylabel('$\eta$')

	plt.subplot(142)
	plt.scatter(eta, chi_eff, s=1*match**8, c=np.log10(1-match))
	plt.plot(eta_targ, chi_eff_targ, 'k<')
	plt.plot(eta[max_idx], chi_eff[max_idx], 'rx')
	plt.xlim(min(eta), max(eta))
	plt.ylim(min(chi_eff), max(chi_eff))
	plt.xlabel('$\eta$')
	plt.ylabel('$\chi_\mathrm{eff}$')

	plt.subplot(143)
	plt.scatter(s1z, s2z, s=1*match**8, c=np.log10(1-match))
	plt.plot(s1z_targ, s2z_targ, 'k<')
	plt.plot(s1z[max_idx], s2z[max_idx], 'rx')
	plt.xlim(min(s1z), max(s1z))
	plt.ylim(min(s2z), max(s2z))
	plt.xlabel('$\chi_1$')
	plt.ylabel('$\chi_2$')

	plt.subplot(144)
	plt.scatter(Mf, af, s=1*match**8, c=np.log10(1-match))
	plt.plot(mf_targ, af_targ, 'k<')
	plt.plot(Mf[max_idx], af[max_idx], 'rx')
	plt.xlim(min(Mf), max(Mf))
	plt.ylim(min(af), max(af))
	plt.xlabel('$M_f~[M_\odot]$')
	plt.ylabel('$a_f$')
	plt.colorbar()

	plt.tight_layout()
	plt.savefig('%s/%s.png' %(data_dir, fname), dpi=300)
	plt.close()

