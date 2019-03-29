import matplotlib as mpl
mpl.use('Agg')
import commands, os, sys, numpy as np
import matplotlib.pyplot as plt
import imrtestgr as tgr
import nr_fits as nr
import scipy
from scipy import interpolate

postloc_list = commands.getoutput('ls ../runs/waveform_systematics/IMRTGR_WaveformSystStudy_2017-08-31_0.99-1.01_range/SEOBNRv2_*/a*/*/*_29.00.dat').split()

Mf_targ, af_targ = tgr.calc_final_mass_spin(36., 29., 0., 0., 'nonprecspin_Healy2014')

for postloc in  postloc_list:
	print postloc
	data = np.genfromtxt(postloc, names=True, dtype=None, skip_header=3)
	m1 = data['m1']
	m2 = data['m2']
	match = data['match']
	Mf, af = tgr.calc_final_mass_spin(m1, m2, np.zeros(len(m1)), np.zeros(len(m1)), 'nonprecspin_Healy2014')
	Mfi = np.linspace(Mf.min(), Mf.max(), len(Mf))
	afi = np.linspace(af.min(), af.max(), len(af))
	matchi = scipy.interpolate.griddata((Mf, af), match, (Mfi[None,:], afi[:,None]), method='cubic')

	plt.figure()
	plt.scatter(Mf, af, c=np.log10(1.-match), lw=0, alpha=1)
	plt.colorbar()
	CS = plt.contour(Mfi, afi, matchi, levels=(0.68, 0.95, 0.99, 0.997))
	plt.clabel(CS, inline=1, fontsize=10)
	plt.axvline(x=Mf_targ)
	plt.axhline(y=af_targ)
	plt.xlabel('$M_f$')
	plt.ylabel('$a_f$')
	plt.title('$M_f$=%.2f; $a_f$=%.2f; colorbar: log$_{10}$(match)'%(Mf_targ, af_targ))
	plt.grid()
	plt.xlim([Mf.min(), Mf.max()])
	plt.ylim([af.min(), af.max()])
	plt.savefig(postloc.replace('.dat', '_Mfaf.png'))

