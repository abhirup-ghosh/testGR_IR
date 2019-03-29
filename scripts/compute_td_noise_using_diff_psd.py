import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import td_likelihood_utils as uc
import os

f_start = 10.
srate = 4096.

adligo_zdhp_asd_loc = '/home/abhirup/src/lalsuite/lalsimulation/src/LIGO-T0900288-v3-ZERO_DET_high_P.txt'
adV_design_asd_loc = '/home/abhirup/src/lalsuite/lalsimulation/src/LIGO-P1200087-v18-AdV_DESIGN.txt'

adL_psd = uc.psd_from_txt(adligo_zdhp_asd_loc, delta_f=1.0/16., flow=10., f_high=1024.); adL_asd = np.sqrt(adL_psd)
adV_psd = uc.psd_from_txt(adV_design_asd_loc, delta_f=1.0/16., flow=10., f_high=1024.); adV_asd = np.sqrt(adV_psd)
f = np.arange(0., 1024., 1./16.)

adL_analytc_f, adL_analytc_psd = uc.adLIGO_psd(f_start, srate); adL_analytc_asd = np.sqrt(adL_analytc_psd)
adV_analytc_f, adV_analytc_asd = uc.adVIRGO_asd(f_start, srate)

ce_f, ce_psd = uc.CE_psd(f_start, srate); ce_asd = np.sqrt(ce_psd)
et_f, et_psd = uc.ET_psd(f_start, srate); et_asd = np.sqrt(et_psd)

os.system('cp %s ../PSD'%adligo_zdhp_asd_loc)
os.system('cp %s ../PSD'%adV_design_asd_loc)
np.savetxt('../PSD/CE_Design_ASD_Dwyer_etal_2015.dat', np.c_[ce_f, ce_asd])
np.savetxt('../PSD/ET_Design_ASD_Mishra_etal_2010.dat', np.c_[et_f, et_asd])

plt.figure(figsize=(5,5))
plt.loglog(f, adL_asd[:-1], label='adL-asd', color='b')
plt.loglog(f, adV_asd[:-1], label='adV-asd', color='r')
plt.loglog(adL_analytc_f, adL_analytc_asd, label='adL-analytc-asd', color='g')
plt.loglog(adV_analytc_f, adV_analytc_asd, label='adV-analytc-asd', color='orange')
plt.loglog(ce_f, ce_asd, label='ce-asd', color='c')
plt.loglog(et_f, et_asd, label='et-asd', color='k')

####################################################################
# Estimating PSDs
####################################################################
adL_asd_loc = '/home/abhirup/Documents/Work/testGR_IR/PSD/LIGO-T0900288-v3-ZERO_DET_high_P.txt'
adL_td_noise = uc.TDnoise(delta_t=1./2048., noise_length=16., asd_file=adL_asd_loc, flow=10., f_high=1024.)
adL_psd_estimated = uc.PSD_from_TDnoise(adL_td_noise, delta_t=1./2048.)
adL_asd_estimated = np.sqrt(adL_psd_estimated)

ce_asd_loc = '/home/abhirup/Documents/Work/testGR_IR/PSD/CE_Design_ASD_Dwyer_etal_2015.dat'
ce_td_noise = uc.TDnoise(delta_t=1./2048., noise_length=16., asd_file=ce_asd_loc, flow=10., f_high=1024.)
ce_psd_estimated = uc.PSD_from_TDnoise(ce_td_noise, delta_t=1./2048.)
ce_asd_estimated = np.sqrt(ce_psd_estimated)

et_asd_loc = '/home/abhirup/Documents/Work/testGR_IR/PSD/ET_Design_ASD_Mishra_etal_2010.dat'
et_td_noise = uc.TDnoise(delta_t=1./2048., noise_length=16., asd_file=et_asd_loc, flow=10., f_high=1024.)
et_psd_estimated = uc.PSD_from_TDnoise(et_td_noise, delta_t=1./2048.)
et_asd_estimated = np.sqrt(et_psd_estimated)

plt.loglog(adL_psd_estimated.sample_frequencies, adL_asd_estimated, alpha=0.2, color='b' )
plt.loglog(ce_psd_estimated.sample_frequencies, ce_asd_estimated, alpha=0.2, color='c' )
plt.loglog(et_psd_estimated.sample_frequencies, et_asd_estimated, alpha=0.2, color='k' )
plt.legend(loc='best')
plt.xlim([10, 1024])
plt.ylim([1e-25, 1e-20])
plt.grid(axis='both', which='minor')
plt.savefig('/home/abhirup/Documents/Work/testGR_IR/papers/imrtgr_td_paper/img/PSDs.png')
