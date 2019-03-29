#!/usr/bin/env python

import os, math, numpy as np, copy
from scipy import random, interpolate
from lal import MSUN_SI, MTSUN_SI, PC_SI, C_SI
import lal, lalsimulation as lalsim
from pycbc.detector import Detector
from pycbc.filter import overlap
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import *
from lalframe import FrWriteREAL8TimeSeries
from optparse import OptionParser
import modGR_waveforms as modGR
import waveforms as wf
import noiseutils as nu

if __name__ == '__main__':
  
  parser = OptionParser()
  parser.add_option('-D', '--data-file', type='string', dest='data_file', help='NR data file to read')
  parser.add_option('-M', '--Mtot', type='float', dest='M', help='total mass', default=70.)
  parser.add_option('--dist', dest='dist', type='float', help='distance', default=420.)
  parser.add_option('--incl', dest='incl', type='float', help='inclination', default=2.54)
  parser.add_option('--ra', dest='ra', type='float', help='right ascension', default=1.76)
  parser.add_option('--dec', dest='dec', type='float', help='declination', default=-1.23)
  parser.add_option('--psi', dest='psi', type='float', help='polarization', default=1.60)
  parser.add_option('--fix-snr', dest='inj_snr', type='float', help='SNR (if provided, overrides distance)', default=None)
  parser.add_option('--phi-ref', dest='phi_ref', type='float', help='phase at reference frequency', default=0.)
  parser.add_option('--f-ref', dest='f_ref', type='float', help='reference frequency', default=100.)
  parser.add_option('--f-low', dest='f_low', type='float', help='lower frequency cut-off (default = 20 Hz)', default=20.)
  parser.add_option('--f-start', dest='f_start', type='float', help='lower frequency for SNR computation (default=f_low)', default=None)
  parser.add_option('--srate', dest='srate', type='float', help='sampling rate (default = 4096 Hz)', default=4096.)
  parser.add_option('--seglen', dest='seglen', type='int', help='length of frame to be injected', default=None)
  parser.add_option('--trig-time', dest='trig_time', type='float', help='trigger time', default=1126285216)
  parser.add_option('-n', '--noise-model', dest='noise_model', type='string', help='noise model from lalinspiral.sbank.psds / lalsimulation', default='aLIGO_EARLY_HIGH')
  parser.add_option('--data-seed', dest='data_seed', type='int', help='random number seed for noise', default=None)
  parser.add_option('--no-noise', action='store_false', dest='add_noise', help='turn off addition of noise', default=True)
  parser.add_option('-c', '--channel-suffix', dest='channel_suffix', type='string', help='suffix for channel name', default='modGR_INJECTED')
  parser.add_option('-i', '--in-folder', dest='in_folder', type='string', help='input folder to read raw waveform', default='../waveforms')
  parser.add_option('-o', '--out-folder', dest='out_folder', help='output folder to dump frame files', default='../frames')
  parser.add_option('--disable-plot', action='store_false', dest='plotting', help='turn off plotting', default=True)
  
  (options, args) = parser.parse_args()
  
  data_file = options.data_file
  M = options.M
  dist = options.dist
  iota = options.incl
  ra = options.ra
  dec = options.dec
  psi = options.psi
  phi_ref = options.phi_ref
  f_ref = options.f_ref
  f_low = options.f_low
  f_start = options.f_start
  srate = options.srate
  seglen = options.seglen
  trig_time = options.trig_time
  add_noise = options.add_noise
  noise_model = options.noise_model
  data_seed = options.data_seed
  channel_suffix = options.channel_suffix
  in_folder = options.in_folder
  inj_snr = options.inj_snr
  out_folder = options.out_folder
  plotting = options.plotting
  
  if data_seed is None:
    data_seed = random.randint(32768)
  
  random.seed(data_seed)
  
  if add_noise:
    noise_str = '%s_seed%d'%(noise_model, data_seed)
  else:
    noise_str = 'nonoise'
  
  if f_start is None:
    f_start = f_low
  
  frame_extn = '%s_M%d_dist%d_incl%s_ra%s_dec%s_psi%s_flow%d_%s'%(os.path.basename(data_file).replace('.gz', '').replace('.dat', ''), round(M), round(dist), str(iota).replace('.', 'p').replace('-', 'm'), str(ra).replace('.', 'p').replace('-', 'm'), str(dec).replace('.', 'p').replace('-', 'm'), str(psi).replace('.', 'p').replace('-', 'm'), round(f_low), noise_str)
  
  os.system('mkdir -p %s'%(out_folder))
  
  # Nathan's EOB: load the data
  modGR_wf_object = modGR.modGR_waveform(data_file)
  
  # Create the waveform in the source frame
  (t_arr, hp_arr, hc_arr) = wf.data_from_TD(modGR_wf_object.TDwaveform(M, dist=dist, iota=iota, f_low=f_low, srate=srate, phi_ref=phi_ref, f_ref=f_ref, taper=True, trim_right=True))
 
  ## Make sure that the coalescence is about 2s before the end of the segment and choose segment length to minimum required power of 2
  wf_len = t_arr[-1] - t_arr[0]
  gps_end = trig_time + 2
  if seglen is None:
    seglen = int(2**math.ceil(np.log2(wf_len + 2)))
  else:
    seglen = int(seglen)
  print 'seglen = %d'%(seglen)
  gps_start = gps_end - seglen
  
  for det in ['H1', 'L1']:
    detector = Detector(det)
    channel = '%s:%s'%(det, channel_suffix)
    
    (t_arr_padded, strain_t) = wf.detector_strain(det, ra, dec, psi, trig_time, t_arr, hp_arr, hc_arr, seglen)
    
    # Create an empty strain vector
    strainT = lal.CreateREAL8TimeSeries(channel, gps_start, 0.0, 1.0/srate, lal.DimensionlessUnit, int(seglen*srate))
    
    strainT.data.data = strain_t
    
    ## If injection SNR is provided, scale distance to get appropriate SNR
    ## ### NOTE ### We are not using this yet!
    #if inj_snr is not None:
      #strain_f_data = np.fft.rfft(strainT.data.data)/strainT.data.length
      #ff = np.fft.rfftfreq(strainT.data.length, strainT.deltaT)
      #df = ff[1] - ff[0]
      
      #psd = aLIGOZeroDetHighPower(len(ff), df, f_low/2.) # FIXME
      #strain_f = FrequencySeries(strain_f_data[np.where(ff>f_low)], delta_f=df)  
      #snr_f = 4.0*overlap(strain_f, strain_f, psd=psd, low_frequency_cutoff=f_low, high_frequency_cutoff=1024, normalized=False)
      
      #strainT.data.data *= inj_snr / snr_f
    
    # Retain signal for plotting if necessary
    vars()[det+'_signal'] = copy.deepcopy(strainT.data.data)
    
    # Calculate optimal SNR
    strainF = wf.FD_from_TD(strainT)
    (ff, strain_f_data) = wf.data_from_FD(strainF)
    df = ff[1] - ff[0]
    
    vars()[det+'_snr_opt'] = nu.SNR_from_fd(ff, strain_f_data, noise_model=noise_model, f_low=f_low)
    
    print '%s\toptimal SNR = %.2f'%(det, vars()[det+'_snr_opt'])
    
    # Add noise
    if add_noise:
      noise_TD = nu.TDnoise(strainT.data.length, strainT.deltaT, noise_model=noise_model, f_low=f_low)
      strainT.data.data += noise_TD
      noise_str = noise_model
    else:
      noise_str = 'nonoise'
    
    # Write the strain to a frame file
    FrWriteREAL8TimeSeries(strainT, 0)
    
    default_filename = '%s-%s-%d-%d.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen)
    frame_filename = '%s-%s-%d-%d_%s.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen, frame_extn)
    
    os.system('mkdir -p %s'%(out_folder))
    os.system('mv %s %s'%(default_filename, os.path.join(out_folder, frame_filename)))
  
  print 'Network optimal SNR = %.2f'%(np.sqrt(H1_snr_opt**2 + L1_snr_opt**2))
  
  # Optional plotting
  if plotting:
    from pylal import Fr
    import matplotlib.pyplot as plt
    
    plt.figure()
    
    if not add_noise:
      colors = {'H1': 'r', 'L1': 'c'}
    else:
      colors = {'H1': 'b', 'L1': 'g'}
    
    for det in ['H1', 'L1']:
      
      channel = '%s:%s'%(det, channel_suffix)
      frame_filename = '%s-%s-%d-%d_%s.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen, frame_extn)
      
      frdata, gps_start, xoffset, xspacing, xunit, yunit = Fr.frgetvect(os.path.join(out_folder, frame_filename), channel)
      
      xx = np.arange(0., len(frdata))*xspacing
      plt.plot(xx, frdata, colors[det], label=det)
    
    if add_noise:
      colors = {'H1': 'r', 'L1': 'c'}
      for det in ['H1', 'L1']:
        plt.plot(xx, vars()[det+'_signal'], colors[det], label=det+'_signal')
    
    plt.ylim(-4e-21, 4e-21)
    plt.legend()
    
    plot_filename = '%s-%d-%d_%s.png'%(channel_suffix, int(gps_start), seglen, frame_extn)
    
    plt.savefig(os.path.join(out_folder, plot_filename), dpi=300)
