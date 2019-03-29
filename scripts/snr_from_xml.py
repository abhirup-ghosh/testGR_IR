#!/usr/bin/env python

import sys, numpy as np
import noiseutils as nu, waveforms as wf
from pylal import SimInspiralUtils

xml_file = '../injections/popsynth_injections_SEOBNRv4_ROMpseudoFourPN.xml'
out_file = '../injections/snr_popsynth_injections_SEOBNRv4_ROMpseudoFourPN.dat'

injections = SimInspiralUtils.ReadSimInspiralFromFiles([xml_file])

det_list = ['H1', 'L1']
noise_map = {'H1': 'aLIGO_DESIGN', 'L1': 'aLIGO_DESIGN', 'V1': 'AdV_DESIGN'}
f_start_map = {'H1': 20., 'L1': 20., 'V1':20.}

approx = 'SEOBNRv2'
verbose = True
srate = 4096.
f_low = 10.
phi_ref=0.
f_ref=100.
taper=True

ofile = open(out_file, 'w')
ofile.write('#event %s network_snr\n'%(' '.join(['snr_%s'%(det) for det in det_list])))
sys.stdout.write(verbose*('#event %s network_snr\n'%(' '.join(['snr_%s'%(det) for det in det_list]))))

for (ev, inj) in enumerate(injections):
  # Retrieve the injection parameters
  m1 = inj.mass1
  m2 = inj.mass2
  s1x = inj.spin1x
  s1y = inj.spin1y
  s1z = inj.spin1z
  s2x = inj.spin2x
  s2y = inj.spin2y
  s2z = inj.spin2z
  dist = inj.distance
  incl = inj.inclination
  ra = inj.longitude
  dec = inj.latitude
  psi = inj.polarization
  trig_time = inj.geocent_end_time+1e-9*inj.geocent_end_time_ns
  # Generate the wf in TD
  try:
    (t_arr, hp_arr, hc_arr) = wf.data_from_TD(wf.TDwaveform(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dist, incl, approx, srate=srate, f_low=f_low, phi_ref=phi_ref, f_ref=f_ref, taper=taper))
  except RuntimeError:
    sys.stdout.write('%d waveform generation failed\n'%(ev))
  else:
    # Calculate the SNR
    snrsq_opt = 0.
    for det in det_list:
      noise_model = noise_map[det]
      f_start_det = f_start_map[det]
      (t_arr_padded, strain_t) = wf.detector_strain(det, ra, dec, psi, trig_time, t_arr, hp_arr, hc_arr)
      (ff, strain_f) = nu.fd_from_td(t_arr_padded, strain_t)
      vars()[det+'_snr_opt'] = nu.SNR_from_fd(ff, strain_f, noise_model=noise_model, f_low=f_start_det)
      snrsq_opt += (vars()[det+'_snr_opt'])**2.
    snr_opt = np.sqrt(snrsq_opt)
    ofile.write('%d %s %f\n'%(ev, ' '.join([str(vars()[det+'_snr_opt']) for det in det_list]), snr_opt))
    sys.stdout.write(verbose*('%d %s %f\n'%(ev, ' '.join([str(vars()[det+'_snr_opt']) for det in det_list]), snr_opt)))

ofile.close()
