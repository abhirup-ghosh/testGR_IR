#import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import imrtestgr as tgr
import bayesian as ba

#inspiral location
insp_sch_isco_loc = '/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-06_v4250_fhigh_insp62Hz_flow_ring_134Hz/data'
insp_kerr_isco_loc = '/home/ajith/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-01_v4252_fhigh_insp134Hz_flow_ring_134Hz/data'

#ringdown location
ring_kerr_isco_loc = '/home/ajith/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-01_v4252_fhigh_insp134Hz_flow_ring_134Hz/data'
ring_80per_qnm_loc = '/home/ajith/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-01_v4252_fhigh_insp134Hz_flow_ring_201Hz/data'
ring_qnm_loc = '/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-06_v4250_fhigh_insp134Hz_flow_ring_251Hz/data'

#imr location
imr_loc = '/home/ajith/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DoubleSpin_nonprecspin_Healy2014_2015-10-01_v4252_fhigh_insp134Hz_flow_ring_134Hz/data'

#locations
insp_loc_list = [insp_sch_isco_loc, insp_kerr_isco_loc]
ring_loc_list = [ring_kerr_isco_loc, ring_80per_qnm_loc]
imr_loc_list = [imr_loc]

#cutoff frequencies
insp_fhigh_list = ['62Hz', '134Hz']
ring_flow_list = ['134 Hz', '201Hz']

#specifying linestyles
ls = ['solid','dashed','dashdot','dotted']


plt.figure(figsize=(5,5))
# plotting inspiral contours
for (i, (insp_fhigh, insp_loc)) in enumerate(zip(insp_fhigh_list, insp_loc_list)):
	Mf_bins,af_bins = np.loadtxt('%s/Mfaf.dat'%insp_loc)
	P_Mfaf_i = np.loadtxt('%s/P_Mfaf_i.dat'%insp_loc)
	s1_Mfaf_i = ba.nsigma_value(tgr.gf(P_Mfaf_i), 0.68)
	CSi = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_i), levels=(s1_Mfaf_i,), linewidths=(1,), linestyles=ls[i], colors='orange')	
	strs_i = [insp_fhigh[i]]
        fmt_i = {}
        for l,s in zip(CSi.levels, strs_i):
                fmt_i[l] = s	
	plt.clabel(CSi,CSi.levels[::2],inline=True,fmt=fmt_i,fontsize=14, use_clabeltext=True)
	plt.hold(True)


# plotting ringdown contours
for (i, (ring_flow, ring_loc)) in enumerate(zip(ring_flow_list, ring_loc_list)):
        Mf_bins,af_bins = np.loadtxt('%s/Mfaf.dat'%ring_loc)
        P_Mfaf_r = np.loadtxt('%s/P_Mfaf_r.dat'%ring_loc)
        s1_Mfaf_r = ba.nsigma_value(tgr.gf(P_Mfaf_r), 0.68)
        CSr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_r), levels=(s1_Mfaf_r, ), linewidths=(1, ), linestyles=ls[i], colors='red')
        strs_r = [ring_flow[i]]
        fmt_r = {}
        for l,s in zip(CSr.levels, strs_r):
                fmt_r[l] = s    
        plt.clabel(CSr,CSr.levels[::2],inline=True,fmt=fmt_r,fontsize=14, use_clabeltext=True)
	plt.hold(True)


# plotting imr contours
for imr_loc in imr_loc_list:
        Mf_bins,af_bins = np.loadtxt('%s/Mfaf.dat'%imr_loc)
        P_Mfaf_imr = np.loadtxt('%s/P_Mfaf_imr.dat'%imr_loc)
        s1_Mfaf_imr = ba.nsigma_value(tgr.gf(P_Mfaf_imr), 0.68)
        CSimr = plt.contour(Mf_bins[:-1], af_bins[:-1], tgr.gf(P_Mfaf_imr), levels=(s1_Mfaf_imr, ), linewidths=(1, ), colors='k')
        plt.hold(True)
plt.xlabel('$M_f(M_{\odot})$')
plt.ylabel('$a_f/M_f$')
plt.xlim([40,100])
plt.ylim([0.2,1])
plt.grid()
plt.savefig('/home/abhirup/public_html/imrtestgr/ER8/G184098/IR_overlap_frequency_variation.png')
