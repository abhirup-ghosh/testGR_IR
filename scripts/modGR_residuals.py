#!/usr/bin/env python

import os, numpy as np
from pylal import Fr
import waveforms as wf
import matplotlib.pyplot as plt

runtag = 'modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed22347'

nonoise_runtag = 'modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_nonoise'

frame_loc = '../frames'
pos_loc = '../runs/simulations_modGR/2015-12-16'
approx = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN'
seglen = 4
trigtime = 1126285216.0

H1_frame = os.path.join(frame_loc, 'H-H1_%s.gwf'%(nonoise_runtag))
L1_frame = os.path.join(frame_loc, 'L-L1_%s.gwf'%(nonoise_runtag))

pos_file = os.path.join(pos_loc, runtag, '%s_seglen%d_nomargphi'%(approx, seglen), 'IMR', '%.1f-0'%(trigtime), 'H1L1', 'posterior_samples.dat')

for det in ['H1', 'L1']:
  vars()[det+'_channel'] = '%s:modGR_INJECTED'%(det)
  vars()[det+'_frdata'], vars()[det+'_gps_start'], vars()[det+'_xoffset'], vars()[det+'_xspacing'], vars()[det+'_xunit'], vars()[det+'_yunit'] = Fr.frgetvect(vars()[det+'_frame'], vars()[det+'_channel'])

pos_data = np.genfromtxt(pos_file, dtype=None, names=True)

for param in ['m1', 'm2', 'a1z', 'a2z', 'distance', 'ra', 'dec', 'costheta_jn', 'psi', 'time']:
  vars()[param+'_median'] = np.median(pos_data[param])

(H1_tarr, H1_bestfitwf) = wf.detector_strain('H1', ra_median, dec_median, psi_median, time_median, *wf.data_from_TD(wf.TDwaveform(m1=m1_median, m2=m2_median, s1z=a1z_median, s2z=a2z_median, dist=distance_median, incl=costheta_jn_median, approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))

(L1_tarr, L1_bestfitwf) = wf.detector_strain('L1', ra_median, dec_median, psi_median, time_median, *wf.data_from_TD(wf.TDwaveform(m1=m1_median, m2=m2_median, s1z=a1z_median, s2z=a2z_median, dist=distance_median, incl=costheta_jn_median, approx='SEOBNRv2', f_low=20., srate=4096., taper=False)))
