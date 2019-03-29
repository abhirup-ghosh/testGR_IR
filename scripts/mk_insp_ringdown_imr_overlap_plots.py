"""
Example Command:
python mk_insp_ringdown_imr_overlap_plots.py -d /home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/IMR_consistency_robustness_tests/imrtestgr/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-26_fhigh_insp50Hz_flow_ring_50Hz_/data -d /home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/IMR_consistency_robustness_tests/imrtestgr/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-26_fhigh_insp95Hz_flow_ring_95Hz_/data -d /home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/IMR_consistency_robustness_tests/imrtestgr/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-26_fhigh_insp100Hz_flow_ring_100Hz_/data -d /home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/IMR_consistency_robustness_tests/imrtestgr/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-26_fhigh_insp150Hz_flow_ring_150Hz_/data -l '(50 Hz)' -l '(95 Hz)' -l '(100 Hz)' -l '(150 Hz)' --fit-formula 'nonprecspin_Healy2014' --m1-inj 50 --m2-inj 50 --chi1-inj 0 --chi2-inj 0 --tag frequency -o .
"""

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import bayesian as ba
import imrtestgr as tgr
import nr_fits as nr
import os
import scipy.ndimage.filters as filter
from optparse import OptionParser
import plotsettings

parser = OptionParser()
parser.add_option("-d", "--postloc-list", dest="post_loc_list", action="append", help="paths to the all the relevant data folders")
parser.add_option("-o", "--outdir", dest="out_dir", help="output directory")
parser.add_option("-l", "--label-list", dest="label_list", action="append", help="label list")
parser.add_option("-t", "--tag", dest="tag", help="tag: approximant, fitformula, frequency")
parser.add_option("--m1-inj", dest="m1_inj", help="injection value of mass 1", default=1.)
parser.add_option("--m2-inj", dest="m2_inj", help="injection value of mass 2", default=1.)
parser.add_option("--chi1-inj", dest="chi1_inj", help="injection value of spin1z", default=0.)
parser.add_option("--chi2-inj", dest="chi2_inj", help="injection value of spin2z", default=0.)
parser.add_option("-f", "--fit-formula", dest="fit_formula", help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_Pan2011', 'nonprecspin_Healy2014'")
parser.add_option("--norm-by", dest="norm_by", help=" [options: norm_by_imr, norm_by_mean]")
(options, args) = parser.parse_args()
post_loc_list = options.post_loc_list
out_dir = options.out_dir
label_list = options.label_list
tag = options.tag
m1_inj = float(options.m1_inj)
m2_inj = float(options.m2_inj)
chi1_inj = float(options.chi1_inj)
chi2_inj = float(options.chi2_inj)
fit_formula = options.fit_formula
norm_by = options.norm_by

Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)

labels = ['line1', 'line2', 'line3', 'line4', 'line5']

#specifying linestyles
ls = ['solid']#['solid','dashed','dotted','dashdot', ':']
color = ['r', 'orange', 'y', 'g', 'c' , 'b', 'k', 'G']
alpha = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]#np.linspace(0.2, 1., 8)

fig = plt.figure(figsize=(10,5))
fig.subplots_adjust(wspace=0.5)
plt.subplot(122)
# plotting inspiral contours
for (i, post_loc) in enumerate(post_loc_list):
        Mf_bins,af_bins = np.loadtxt('%s/Mfchif.dat.gz'%post_loc)
        P_Mfaf_i = np.loadtxt('%s/P_Mfchif_i.dat.gz'%post_loc)
        s1_Mfaf_i = ba.nsigma_value(tgr.gf(P_Mfaf_i), 0.68)
        CSi = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_i), levels=(s1_Mfaf_i,), linewidths=(1,), linestyles=ls[0], colors='orange', alpha=alpha[i])
        CSi.collections[0].set_label('Insp' + label_list[i])
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, frameon=False, fontsize='large')
        plt.hold(True)

# plotting ringdown contours
for (i, post_loc) in enumerate(post_loc_list):
        Mf_bins,af_bins = np.loadtxt('%s/Mfchif.dat.gz'%post_loc)
        P_Mfaf_r = np.loadtxt('%s/P_Mfchif_r.dat.gz'%post_loc)
        s1_Mfaf_r = ba.nsigma_value(tgr.gf(P_Mfaf_r), 0.68)
        CSr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_r), levels=(s1_Mfaf_r, ), linewidths=(1, ), linestyles=ls[0], colors='r', alpha=alpha[i])
        CSr.collections[0].set_label('Ring' + label_list[i])
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, frameon=False, fontsize='large')
        plt.hold(True)

# plotting imr contours
for (i, post_loc) in enumerate(post_loc_list):
        Mf_bins,af_bins = np.loadtxt('%s/Mfchif.dat.gz'%post_loc)
        P_Mfaf_imr = np.loadtxt('%s/P_Mfchif_imr.dat.gz'%post_loc)
        s1_Mfaf_imr = ba.nsigma_value(tgr.gf(P_Mfaf_imr), 0.68)
        CSimr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_imr), levels=(s1_Mfaf_imr, ), linewidths=(1, ), linestyles=ls[0], colors='k', alpha = alpha[i])
        CSimr.collections[0].set_label('IMR' + label_list[i])
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, frameon=False, fontsize='large')

plt.xlabel('$M_f~[M_{\odot}]$', fontsize='xx-large', labelpad=10)
plt.ylabel('$a _f$', fontsize='xx-large', labelpad=10)
#plt.axvline(x=Mf_inj, ls='--', color='k')
#plt.axhline(y=chif_inj, ls='--', color='k')
plt.xlim([40,100])
plt.ylim([0.2,1])
#plt.tight_layout()
#plt.savefig('%s/imr_overlap_%s.pdf'%(out_dir, tag), dpi=300, bbox_inches='tight')

# plotting dMfbyMf dchifbychif plots

#plt.figure(figsize=(5,5))
plt.subplot(121)
for (i, post_loc) in enumerate(post_loc_list):
	if norm_by == 'norm_by_imr':
                P_dMfbyMfdchifbychif = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz' %post_loc)
        elif norm_by == 'norm_by_mean':
                P_dMfbyMfdchifbychif = np.loadtxt('%s/P_dMfbyMfmean_dchifbychifmean.dat.gz' %post_loc)
	dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat.gz' %post_loc)
        dchifbychif_vec = np.loadtxt('%s/dchifbychif_vec.dat.gz' %post_loc)
        s1_dMfbyMfdchifbychif = ba.nsigma_value(tgr.gf(P_dMfbyMfdchifbychif), 0.68)
        CSdMfbyMfdchifbydchif = plt.contour(dMfbyMf_vec, dchifbychif_vec, tgr.gf(P_dMfbyMfdchifbychif), levels=(s1_dMfbyMfdchifbychif, ), linewidths=(1, ), colors=color[i])
        CSdMfbyMfdchifbydchif.collections[0].set_label(label_list[i])
        plt.legend(loc='upper left', frameon=False, fontsize='large')
        plt.hold(True)
plt.xlabel('$\Delta M_f/\\bar{M_f}$', fontsize='xx-large', labelpad=10)
plt.ylabel('$\Delta a_f/\\bar{a_f}$', fontsize='xx-large', labelpad=10)
plt.plot(0, 0, 'k+', ms=12, mew=2)
plt.xlim([-1, 1])
plt.ylim([-1, 1])
plt.tight_layout()
plt.savefig('%s/imr_overlap_%s_%s.pdf'%(out_dir, tag, norm_by), dpi=300, bbox_inches='tight')
#plt.savefig('%s/dMfbyMfdafbydaf_%s.pdf'%(out_dir, tag), dpi=300)
