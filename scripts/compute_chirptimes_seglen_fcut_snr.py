import sys, lal, os, commands, numpy as np
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import imrtestgr as tgr
import nr_fits as nr
import matplotlib.pyplot as plt
import optparse as op
import chirptimes as ct
from pylal import SimInspiralUtils
import noiseutils as nu, waveforms as wf


xml_file = '/home/abhirup/Documents/Work/testGR_IR/runs/LIGO_software_injections_H1L1/1128039919_1128072513/SEOBNRv4_ROMpseudoFourPN_128039919_1128072513_m_10_80_SNR_50_100_LALAdLIGO_H1L1_20Hz.xml'

# read injection file
injections = SimInspiralUtils.ReadSimInspiralFromFiles([xml_file])

approximant = 'SEOBNRv4'
verbose = True
srate = 2048.
f_low_wf = 5.
f_start_ligo = 20.#15.
f_start_virgo = 15.
f_high = srate/2.
phi_ref=0.
f_ref=100.
taper=True
ligo_psd = 'aLIGO_EARLY_HIGH'#'ZERO_DET_high_P'
virgo_psd = 'AdV_DESIGN'

det_list = ['H1', 'L1']#['H1', 'L1', 'V1']
noise_map = {'H1': ligo_psd, 'L1': ligo_psd, 'V1': virgo_psd}
f_start_map = {'H1': f_start_ligo, 'L1': f_start_ligo, 'V1': f_start_virgo}

# output file
#outfile = open(xml_file.replace('.xml', '_usefulinfo_flowpsd_15Hz.txt'), 'w')
outfile = open(xml_file.replace('.xml', '_usefulinfo_H1L1_aLIGO_EARLY_HIGH_flowpsd_20Hz.txt'), 'w')

outfile.write('# event_id m1 m2 a1x a2x a1y a2y a1z a2z longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase ctime seglen f_isco %s network_snr_IMR %s network_snr_insp %s network_snr_ring\n'%(' '.join(['snr_%s_IMR'%(det) for det in det_list]), ' '.join(['snr_%s_insp'%(det) for det in det_list]), ' '.join(['snr_%s_ring'%(det) for det in det_list])))
sys.stdout.write(verbose*'# m1 m2 a1x a2x a1y a2y a1z a2z longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase ctime seglen f_isco %s network_snr_IMR %s network_snr_insp %s network_snr_ring\n'%(' '.join(['snr_%s_IMR'%(det) for det in det_list]), ' '.join(['snr_%s_insp'%(det) for det in det_list]), ' '.join(['snr_%s_ring'%(det) for det in det_list])))

for (ev, inj) in enumerate(injections):
    # Retrieve the injection parameters
    m1 = inj.mass1
    m2 = inj.mass2
    a1x = inj.spin1x
    a1y = inj.spin1y
    a1z = inj.spin1z
    a2x = inj.spin2x
    a2y = inj.spin2y
    a2z = inj.spin2z
    distance = inj.distance
    inclination = inj.inclination
    longitude = inj.longitude
    latitude = inj.latitude
    polarization = inj.polarization
    coa_phase = inj.coa_phase
    geocent_end_time = inj.geocent_end_time
    geocent_end_time_ns = inj.geocent_end_time_ns
    trig_time = geocent_end_time+1e-9*geocent_end_time_ns

    # Kerr ISCO frequency computation: Spinning case
    Mf, af = tgr.calc_final_mass_spin(m1, m2, a1z, a2z, 'nonprecspin_Healy2014')
    f_isco = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

    # Chirptime, seglen computation
    ctime = ct.calc_chirptime(m1, m2, f_min=f_start_ligo)
    seglen = 2**np.ceil(np.log2(ctime + 2.)) # lalinference requires a buffer of 2s between the end of the signal and end of the data

    # SNR computation: IMR, inspiral and ringdown
    (t_arr, hp_arr, hc_arr) = wf.data_from_TD(wf.TDwaveform(m1, m2, a1x, a1y, a1z, a2x, a2y, a2z, distance, inclination, approximant, srate=srate, f_low=f_low_wf, phi_ref=phi_ref, f_ref=f_ref, taper=taper))

    snrsq_opt_IMR = 0.
    snrsq_opt_insp = 0.
    snrsq_opt_ring = 0.

    for det in det_list:
      noise_model = noise_map[det]
      f_start_det = f_start_map[det]
      (t_arr_padded, strain_t) = wf.detector_strain(det, longitude, latitude, polarization, trig_time, t_arr, hp_arr, hc_arr)
      (ff, strain_f) = nu.fd_from_td(t_arr_padded, strain_t)
      vars()[det+'_snr_opt_IMR'] = nu.SNR_from_fd(ff, strain_f, noise_model=noise_model, f_low=f_start_det)
      vars()[det+'_snr_opt_insp'] = nu.SNR_from_fd(ff, strain_f, noise_model=noise_model, f_low=f_start_det, f_high=f_isco)
      vars()[det+'_snr_opt_ring'] = nu.SNR_from_fd(ff, strain_f, noise_model=noise_model, f_low=f_isco)
      snrsq_opt_IMR += (vars()[det+'_snr_opt_IMR'])**2.
      snrsq_opt_insp += (vars()[det+'_snr_opt_insp'])**2.
      snrsq_opt_ring += (vars()[det+'_snr_opt_ring'])**2.
  
    snr_opt_IMR = np.sqrt(snrsq_opt_IMR)
    snr_opt_insp = np.sqrt(snrsq_opt_insp)
    snr_opt_ring = np.sqrt(snrsq_opt_ring)

    # Writing to outfile
    outfile.write('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.2f %.2f %.2f %.2f %d %d %.4f %.2f %d %.2f %s %.2f %s %.2f %s %.2f\n'%(ev+1, m1, m2, a1x, a2x, a1y, a2y, a1z, a2z, longitude, latitude, inclination, polarization, distance, geocent_end_time, geocent_end_time_ns, coa_phase, ctime, seglen, f_isco, ' '.join([str(vars()[det+'_snr_opt_IMR']) for det in det_list]), snr_opt_IMR, ' '.join([str(vars()[det+'_snr_opt_insp']) for det in det_list]), snr_opt_insp, ' '.join([str(vars()[det+'_snr_opt_ring']) for det in det_list]), snr_opt_ring))
    sys.stdout.write(verbose*('%d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.2f %.2f %.2f %.2f %d %d %.4f %.2f %d %.2f %s %.2f %s %.2f %s %.2f\n'%(ev + 1, m1, m2, a1x, a2x, a1y, a2y, a1z, a2z, longitude, latitude, inclination, polarization, distance, geocent_end_time, geocent_end_time_ns, coa_phase, ctime, seglen, f_isco, ' '.join([str(vars()[det+'_snr_opt_IMR']) for det in det_list]), snr_opt_IMR, ' '.join([str(vars()[det+'_snr_opt_insp']) for det in det_list]), snr_opt_insp, ' '.join([str(vars()[det+'_snr_opt_ring']) for det in det_list]), snr_opt_ring)))



outfile.close()
