#!/usr/bin/env python

import os, numpy as np
from pylal import Fr
import waveforms as wf # This should be in the src
import matplotlib.pyplot as plt
import noiseutils as nu
from copy import copy
import modGR_waveforms as modGR

# A list of locations where things are
datafile_GR = '/home/archis/Work/testGR_IR/waveforms/q1_GR_new.dat.gz'
datafile_modGR = '/home/archis/Work/testGR_IR/waveforms/q1_a2_400_r0_32.dat'

H1_frame_GR_nonoise = '/home/archis/Work/testGR_IR/frames/H-H1_modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise.gwf'
L1_frame_GR_nonoise = '/home/archis/Work/testGR_IR/frames/L-L1_modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise.gwf'

H1_frame_modGR_nonoise = '/home/archis/Work/testGR_IR/frames/H-H1_modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise.gwf'
L1_frame_modGR_nonoise = '/home/archis/Work/testGR_IR/frames/L-L1_modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise.gwf'

recstat_GR_nonoise = '/home/archis/Work/testGR_IR/runs/simulations_modGR/2015-11-27/modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/IMR/1126285216.0-0/H1L1/summary_statistics.dat'
recstat_modGR_nonoise = '/home/archis/Work/testGR_IR/runs/simulations_modGR/2015-12-07/modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/IMR/1126285216.0-0/H1L1/summary_statistics.dat'

# Set the files from one of the above
datafile = datafile_GR

H1_frame = H1_frame_GR_nonoise
L1_frame = L1_frame_GR_nonoise

recstat = recstat_GR_nonoise

## Load the raw waveform ### THIS TAKES TIME ###
#(t_geom, h22r_geom, h22i_geom) = np.loadtxt(data_file, skiprows=1, unpack=True)

# Get the data from the frames
for det in ['H1', 'L1']:
  vars()[det+'_channel'] = '%s:modGR_INJECTED'%(det)
  vars()[det+'_frdata'], vars()[det+'_gps_start'], vars()[det+'_xoffset'], vars()[det+'_xspacing'], vars()[det+'_xunit'], vars()[det+'_yunit'] = Fr.frgetvect(vars()[det+'_frame'], vars()[det+'_channel'])

# Get the 'maP' parameters
data = np.genfromtxt(recstat, dtype=None, names=True, unpack=False)
var_names = [d[0] for d in data]
stat_names = data.dtype.names
for param in ['m1', 'm2', 'a1z', 'a2z', 'distance', 'ra', 'dec', 'costheta_jn', 'psi', 'time']:
  vars()[param+'_maP'] = data[var_names.index(param)][stat_names.index('maP')+1]

# Get the best fit waveforms from the _maP ### NEEDS PyCBC ###
(H1_tarr, H1_bestfitwf) = wf.detector_strain('H1', ra_maP, dec_maP, psi_maP, time_maP, *wf.data_from_TD(wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=costheta_jn_maP, approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))

(L1_tarr, L1_bestfitwf) = wf.detector_strain('L1', ra_maP, dec_maP, psi_maP, time_maP, *wf.data_from_TD(wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=costheta_jn_maP, approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))

## Write the files / make plots

#plt.figure()
#plt.plot(H1_tarr, H1_frdata)
#plt.plot(H1_tarr, H1_bestfitwf)

#plt.figure()
#plt.plot(L1_tarr, L1_frdata)
#plt.plot(L1_tarr, L1_bestfitwf)

#np.savetxt('../residuals/H1_frame_GR_nonoise', np.array([H1_tarr, H1_frdata]), comments='# t hoft')
#np.savetxt('../residuals/L1_frame_GR_nonoise', np.array([L1_tarr, L1_frdata]), comments='# t hoft')

#np.savetxt('../residuals/H1_bestfitwf_GR_nonoise', np.array([H1_tarr, H1_bestfitwf]), comments='# t hoft')
#np.savetxt('../residuals/L1_bestfitwf_GR_nonoise', np.array([L1_tarr, L1_bestfitwf]), comments='# t hoft')

# Frequency domain
(H1_farr_inj, H1_frdata_f) = nu.fd_from_td(H1_tarr, H1_frdata)
(H1_farr_bf, H1_bestfitwf_f) = nu.fd_from_td(H1_tarr, H1_bestfitwf)

#Copy
H1_farr_inj_GR = copy(H1_farr_inj)
H1_frdata_f_GR = copy(H1_frdata_f)

H1_farr_bf_GR = copy(H1_farr_bf)
H1_bestfitwf_f_GR = copy(H1_bestfitwf_f)

# Set the files from one of the above
datafile = datafile_modGR

H1_frame = H1_frame_modGR_nonoise
L1_frame = L1_frame_modGR_nonoise

recstat = recstat_modGR_nonoise

## Load the raw waveform ### THIS TAKES TIME ###
#(t_geom, h22r_geom, h22i_geom) = np.loadtxt(data_file, skiprows=1, unpack=True)

# Get the data from the frames
for det in ['H1', 'L1']:
  vars()[det+'_channel'] = '%s:modGR_INJECTED'%(det)
  vars()[det+'_frdata'], vars()[det+'_gps_start'], vars()[det+'_xoffset'], vars()[det+'_xspacing'], vars()[det+'_xunit'], vars()[det+'_yunit'] = Fr.frgetvect(vars()[det+'_frame'], vars()[det+'_channel'])

# Get the 'maP' parameters
data = np.genfromtxt(recstat, dtype=None, names=True, unpack=False)
var_names = [d[0] for d in data]
stat_names = data.dtype.names
for param in ['m1', 'm2', 'a1z', 'a2z', 'distance', 'ra', 'dec', 'costheta_jn', 'psi', 'time']:
  vars()[param+'_maP'] = data[var_names.index(param)][stat_names.index('maP')+1]

# Get the best fit waveforms from the _maP ### NEEDS PyCBC ###
(H1_tarr, H1_bestfitwf) = wf.detector_strain('H1', ra_maP, dec_maP, psi_maP, time_maP, *wf.data_from_TD(wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=np.arccos(costheta_jn_maP), approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))

(L1_tarr, L1_bestfitwf) = wf.detector_strain('L1', ra_maP, dec_maP, psi_maP, time_maP, *wf.data_from_TD(wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=np.arccos(costheta_jn_maP), approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))

# Frequency domain
(H1_farr_inj, H1_frdata_f) = nu.fd_from_td(H1_tarr, H1_frdata)
(H1_farr_bf, H1_bestfitwf_f) = nu.fd_from_td(H1_tarr, H1_bestfitwf)

#Copy
H1_farr_inj_modGR = copy(H1_farr_inj)
H1_frdata_f_modGR = copy(H1_frdata_f)

H1_farr_bf_modGR = copy(H1_farr_bf)
H1_bestfitwf_f_modGR = copy(H1_bestfitwf_f)


plt.figure()
plt.loglog(H1_farr_inj_GR, np.abs(H1_frdata_f_GR), label='H1 inj GR')
plt.loglog(H1_farr_bf_GR, np.abs(H1_bestfitwf_f_GR), label='H1 rec GR')
plt.loglog(H1_farr_inj_modGR, np.abs(H1_frdata_f_modGR), label='H1 inj modGR')
plt.loglog(H1_farr_bf_modGR, np.abs(H1_bestfitwf_f_modGR), label='H1 rec modGR')
plt.xlabel('Frequency (Hz)')
plt.ylabel('StrainF (Hz$^{-1/2}$)')
plt.legend(loc='lower left')
plt.show()

# Injection parameters

M_inj = 70.
dist_inj = 420.
incl_inj = 2.54
ra_inj = 1.76
dec_inj = -1.23
psi_inj = 1.60
time_inj = 1126285216

GR_waveform = modGR.modGR_waveform('../waveforms/q1_GR_new.dat.gz')
modGR_waveform = modGR.modGR_waveform('../waveforms/q1_a2_400_r0_32.dat')

TD_GR = GR_waveform.TDwaveform(M=M_inj, dist=dist_inj, iota=incl_inj, f_low=20., srate=4096, taper=True, trim_right=True)
TD_modGR = modGR_waveform.TDwaveform(M=M_inj, dist=dist_inj, iota=incl_inj, f_low=20., srate=4096, taper=True, trim_right=True)

(t_arr_GR, hp_arr_GR, hc_arr_GR) = wf.data_from_TD(TD_GR)
(t_arr_modGR, hp_arr_modGR, hc_arr_modGR) = wf.data_from_TD(TD_modGR)

h_arr_GR = hp_arr_GR + 1j*hc_arr_GR
phi_GR = np.unwrap(np.angle(h_arr_GR))
Foft_GR = np.gradient(phi_GR)/np.gradient(t_arr_GR)/(2*np.pi)

h_arr_modGR = hp_arr_modGR + 1j*hc_arr_modGR
phi_modGR = np.unwrap(np.angle(h_arr_modGR))
Foft_modGR = np.gradient(phi_modGR)/np.gradient(t_arr_modGR)/(2*np.pi)

# ### maP ###

recstat = recstat_GR_nonoise

# Get the 'maP' parameters
data = np.genfromtxt(recstat, dtype=None, names=True, unpack=False)
var_names = [d[0] for d in data]
stat_names = data.dtype.names
for param in ['m1', 'm2', 'a1z', 'a2z', 'distance', 'ra', 'dec', 'costheta_jn', 'psi', 'time']:
  vars()[param+'_maP'] = data[var_names.index(param)][stat_names.index('maP')+1]

TD_maP_GR = wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=np.arccos(costheta_jn_maP), approx='SEOBNRv2', f_low=20., srate=4096., taper=False)

recstat = recstat_modGR_nonoise

# Get the 'maP' parameters
data = np.genfromtxt(recstat, dtype=None, names=True, unpack=False)
var_names = [d[0] for d in data]
stat_names = data.dtype.names
for param in ['m1', 'm2', 'a1z', 'a2z', 'distance', 'ra', 'dec', 'costheta_jn', 'psi', 'time']:
  vars()[param+'_maP'] = data[var_names.index(param)][stat_names.index('maP')+1]

TD_maP_modGR = wf.TDwaveform(m1=m1_maP, m2=m2_maP, s1z=a1z_maP, s2z=a2z_maP, dist=distance_maP, incl=np.arccos(costheta_jn_maP), approx='SEOBNRv2', f_low=20., srate=4096., taper=False)

(t_arr_maP_GR, hp_arr_maP_GR, hc_arr_maP_GR) = wf.data_from_TD(TD_maP_GR)
(t_arr_maP_modGR, hp_arr_maP_modGR, hc_arr_maP_modGR) = wf.data_from_TD(TD_maP_modGR)

h_arr_maP_GR = hp_arr_maP_GR + 1j*hc_arr_maP_GR
phi_maP_GR = np.unwrap(np.angle(h_arr_maP_GR))
Foft_maP_GR = np.gradient(phi_maP_GR)/np.gradient(t_arr_maP_GR)/(2*np.pi)

h_arr_maP_modGR = hp_arr_maP_modGR + 1j*hc_arr_maP_modGR
phi_maP_modGR = np.unwrap(np.angle(h_arr_maP_modGR))
Foft_maP_modGR = np.gradient(phi_maP_modGR)/np.gradient(t_arr_maP_modGR)/(2*np.pi)

plt.figure()
plt.plot(t_arr_GR, Foft_GR, label='inj GR')
plt.plot(t_arr_maP_GR, Foft_maP_GR, label='rec GR')
plt.plot(t_arr_modGR, Foft_modGR, label='inj modGR')
plt.plot(t_arr_maP_modGR, Foft_maP_modGR, label='rec modGR')
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.show()



