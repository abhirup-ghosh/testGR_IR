"""
m1, m2, chi1, chi2 = 33.5405393466, 33.5405393466, 0.38579031792, 0.119593857525
"""

import matplotlib as mpl
mpl.use('Agg')
import lal, os, numpy as np
import matplotlib.pyplot as plt
import nr_fits as nr

post_loc_c00 = '../notebooks/Run4c_20171210_HLV_IMRPPv2_cleaneddata_widermcprior_fcut161_seglen16BWPSD_posterior_samples.dat'
post_loc_c02 = '../notebooks/Run43_20180619_HLV_IMRPPv2_cleaneddataC02_widermcprior_fcut161_seglen16BWPSD_padding16_posterior_samples.dat'

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

for post_loc in [post_loc_c00, post_loc_c02]:
    data = np.genfromtxt(post_loc, names=True, dtype=None)
    af, Mf = data['af'], data['mf']

    print 'read data'
    print af	
    
    Q, f_qnm_natural = nr.calc_fqnm_dominant_mode(af)
    f_qnm = f_qnm_natural/(Mf*lal.MTSUN_SI)
    
    ax1.hist(Q,bins=100,histtype='step',normed=True)
    ax2.hist(f_qnm,bins=100,histtype='step',normed=True)
    ax2.axvline(x=287.9,color='r')
    
ax1.xlabel('Q')    
ax2.xlabel('fQNM')
plt.tight_layout()    
plt.savefig('./fqnm_Q_histogram.png')
