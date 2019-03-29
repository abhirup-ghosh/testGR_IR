#!/usr/bin/env python

"""
Script to inject NR waveforms to frames

Contributors: Abhirup Ghosh, Archisman Ghosh, Ashok Choudhary, KaWa Tsang, Laura Van Der Schaaf , Nathan K Johnson-McDaniel, Peter Pang
"""

import time, re, os, sys, glob, math, numpy as np, copy
from math import factorial as fac
from scipy import random, interpolate
from lal import MSUN_SI, MTSUN_SI, PC_SI, C_SI
import lal, lalsimulation as lalsim
from lalframe import FrWriteREAL8TimeSeries
from optparse import OptionParser
from operator import itemgetter
from itertools import groupby
import h5py
import matplotlib as mpl
mpl.use('Agg')

sys.path.insert(0, '../src')
import waveforms as wf
import noiseutils as nu

#define spin-weighted Spherical Harmonics with spin weight s = -2 as a function of  l,m, theta(ange with z axis) and phi(angle in x-y plain with x axis ) 
def spinm2_sphharm(l,m,theta,phi):
  k1 = max(0,m-2)
  k2 = min(l+m,l-2)
  a = sum((((-1)**k)*math.sqrt(fac(l+m)*fac(l-m)*fac(l+2)*fac(l-2))*((np.cos((theta)/2))**(2*l+m-2*k-2))*((np.sin((theta)/2))**(2*k-m+2)))/(fac(k)*fac(k-m+2)*fac(l+m-k)*fac(l-k-2)) for k in range(k1,k2+1))
  ylm = math.sqrt((2*l+1)/(4*np.pi))*a*np.exp(m*(phi)*1j)
  return ylm

def _max_range(indices):
  ranges = []
  for k, g in groupby(enumerate(indices), lambda (i,x):i-x):
    group = map(itemgetter(1), g)
    ranges.append((group[0], group[-1]))
  lengths = [r[1]-r[0]+1 for r in ranges]
  return ranges[lengths.index(max(lengths))]

class NR_waveform(object):
  def __init__(self, data_file, data_format='SXS', N=4, verbose=True):
    t_start = time.time()
    self.data_file = data_file
    self.verbose = verbose
    # Load the data
    if data_format=='SXS':
      sys.stdout.write(verbose*('Loading data from %s\n'%(os.path.realpath(data_file))))
      sys.stdout.write(verbose*('Extrapolation order N = %d\n'%(N)))
      ff = h5py.File(data_file, 'r')
      available_modes = ff.get('Extrapolated_N%d.dir'%(N)).keys()
      highest_mode = available_modes[-1] # ### FIXME ###
      self.lmax = int(re.findall('l\d*', highest_mode)[0][1:])
      self.absmmin = 0
      sys.stdout.write(verbose*('lmax = %d\n'%(self.lmax)))
      for l in range(2, self.lmax+1):
        for m in range(-l, l+1):
          (t_geom, hr_geom, hi_geom) = (ff.get('Extrapolated_N%d.dir/Y_l%d_m%d.dat'%(N, l, m)).value).T
          self.__dict__['t%d%s_geom'%(l, repr(m).replace('-', 'm'))] = t_geom
          self.__dict__['h%d%sr_geom'%(l, repr(m).replace('-', 'm'))] = hr_geom
          self.__dict__['h%d%si_geom'%(l, repr(m).replace('-', 'm'))] = hi_geom
    elif data_format=='IHES_modGR':
      sys.stdout.write(verbose*('Loading data from %s\n'%(os.path.realpath(data_file))))
      (t_geom, hr_geom, hi_geom) = np.loadtxt(data_file, skiprows=1, unpack=True)
      self.lmax = 2
      self.absmmin = 2
      self.t22_geom = t_geom
      self.h22r_geom = hr_geom
      self.h22i_geom = hi_geom
      self.h2m2r_geom = hr_geom
      self.h2m2i_geom = -hi_geom
    sys.stdout.write(verbose*('Time taken = %.3fs.\n\n'%(time.time() - t_start)))
    # Compute number of available cycles
    h22_geom = self.h22r_geom - 1j*self.h22i_geom
    phi22_geom = np.unwrap(np.angle(h22_geom))
    Foft_geom = np.gradient(phi22_geom)/np.gradient(self.t22_geom)/(2.*np.pi)
    idx_geom_raw, = np.where(Foft_geom>0.)
    idx_geom = _max_range(idx_geom_raw)
    idx_peak = np.abs(h22_geom[idx_geom[0]:idx_geom[-1]]).argmax()
    self.ncycles_peak = ((phi22_geom[idx_geom[0]:idx_geom[-1]])[idx_peak]- phi22_geom[idx_geom[0]])/(2.*np.pi)
    self.ncycles = (phi22_geom[idx_geom[-1]] - phi22_geom[idx_geom[0]])/(2.*np.pi)
    self.f_min_geom = Foft_geom[idx_geom[0]]
    sys.stdout.write(verbose*('Number of cycles up to h22 peak amplitude = %.2f\n'%(self.ncycles_peak)))
    sys.stdout.write(verbose*('Total number of phase cycles in h22 = %.2f\n'%(self.ncycles)))
    sys.stdout.write(verbose*('M * f_min = %f M_sun Hz\n\n'%(self.f_min_geom/MTSUN_SI)))
  def TDwaveform(self, M=20., dist=1000., incl=0., f_low=20., srate=4096., phi_ref=0., f_ref=None, taper=False, trim_right=False, trim_tolerance=1e-6, lmax=None, absmmin=None, verbose=None):
    if verbose is None:
      verbose = self.verbose
    if lmax is None:
      lmax = self.lmax
    elif lmax>self.lmax:
      lmax = self.lmax
    if absmmin is None:
      absmmin = self.absmmin
    elif absmmin<self.absmmin:
      absmmin = self.absmmin
    if f_ref is None:
      f_ref = f_low
    t_geom = self.t22_geom
    h22r_geom = self.h22r_geom
    h22i_geom = self.h22i_geom
    # Convert from geometric to SI units
    M_geom = M * MTSUN_SI
    t_SI = t_geom * M_geom
    h22r_SI = self.h22r_geom * M_geom / (dist * 1e6 * PC_SI / C_SI)
    h22i_SI = self.h22i_geom * M_geom / (dist * 1e6 * PC_SI / C_SI)
    # Compute the complex h22
    h22_SI = h22r_SI - 1j*h22i_SI
    # Compute the phase and take derivative to get the frequency
    phi_SI = np.unwrap(np.angle(h22_SI))
    Foft_SI = np.gradient(phi_SI)/np.gradient(t_SI)/(2.*np.pi)
    # Calculate phase at f_ref
    phi_from_F = interpolate.interp1d(Foft_SI, phi_SI)
    phi_at_f_ref = phi_from_F(f_ref)
    phi_at_f_start = phi_SI[0]
    # Correct for phase
    phi0 = np.mod(phi_ref + phi_at_f_start - phi_at_f_ref, 2.*np.pi)
    # Sum over all modes
    hh_geom = np.array([0.+0.j]*len(self.t22_geom))
    for l in range(2, lmax+1):
      for m in range(-l, l+1):
        if abs(m) >= absmmin:
          sys.stdout.write(verbose*('Including (%d,%d) mode.\n'%(l, m)))
          hh_geom +=  ((self.__dict__['h%d%sr_geom'%(l, repr(m).replace('-', 'm'))] - 1j*self.__dict__['h%d%si_geom'%(l, repr(m).replace('-', 'm'))]) * spinm2_sphharm(l, m, incl, phi0) + (self.__dict__['h%d%sr_geom'%(l, repr(m).replace('-', 'm'))] + 1j*self.__dict__['h%d%si_geom'%(l, repr(m).replace('-', 'm'))]) * spinm2_sphharm(l, - m, incl, phi0))/2.
    sys.stdout.write(verbose*'\n')
    hh_SI = hh_geom * M_geom / (dist * 1e6 * PC_SI / C_SI)
    # get hp and hc
    hp_SI = np.real(hh_SI)
    hc_SI = np.imag(hh_SI)
    # Find the frequencies in band
    idx_SI_raw, = np.where(Foft_SI>f_low)
    idx_SI = _max_range(idx_SI_raw)
    # trim if asked to
    idx_SI_max = -1
    if trim_right:
      idx_SI_max = np.where(np.abs(hh_SI)/np.max(np.abs(hh_SI))>trim_tolerance)[0][-1]
      #print 'idx_SI_max = %d'%(idx_SI_max)
    # create an array of times
    t_arr = np.arange(t_SI[idx_SI[0]], t_SI[idx_SI_max], 1./srate)
    # create an interpolation
    hp_interp = interpolate.interp1d(t_SI, hp_SI)
    hc_interp = interpolate.interp1d(t_SI, hc_SI)
    # arrays to return
    hp_arr = hp_interp(t_arr)
    hc_arr = hc_interp(t_arr)
    # create the REAL8TimeSeries
    hp = lal.CreateREAL8TimeSeries('hp', 0, f_low, 1./srate, '', len(hp_arr))
    hc = lal.CreateREAL8TimeSeries('hc', 0, f_low, 1./srate, '', len(hp_arr))
    hp.data.data = hp_arr
    hc.data.data = hc_arr
    # taper waveforms Eq. (3.35) of gr-qc/0001023
    if taper:
      lalsim.SimInspiralREAL8WaveTaper(hp.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
      lalsim.SimInspiralREAL8WaveTaper(hc.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    return hp, hc


if __name__ == '__main__':
  
  parser = OptionParser()
  parser.add_option('-D', '--data-file', type='string', dest='data_file', help='NR data file to read')
  parser.add_option('--data-format', dest='data_format', type='string', help='Available options: SXS (default), IHES_modGR, LALSim', default='SXS')
  parser.add_option('-N', '--extrapolation-order', dest='N', type='int', help='extrapolation order N (default = 4)', default=4)
  parser.add_option('--l-max', dest='lmax', type='int', help='maximum l for injecting higher modes (default = max available; pass 2 for only (2,2)-mode)', default=None)
  parser.add_option('--abs-m-min', dest='absmmin', type='int', help='minimum abs(m) for injecting higher modes (default = zero, that is all modes; pass 2 for only (2,2)-mode)', default=0)
  parser.add_option('-M', '--Mtot', type='float', dest='M', help='total mass in solar masses (REQUIRED ARGUMENT)', default=70.)
  parser.add_option('-q', '--massratio', type='float', dest='q', help='mass ratio (for LALSimulation injections only)', default=1.)
  parser.add_option('--s1x', type='float', dest='s1x', help='x-component of spin of primary (for LALSimulation injections only)', default=0.)
  parser.add_option('--s1y', type='float', dest='s1y', help='y-component of spin of primary (for LALSimulation injections only)', default=0.)
  parser.add_option('--s1z', type='float', dest='s1z', help='e-component of spin of primary (for LALSimulation injections only)', default=0.)
  parser.add_option('--s2x', type='float', dest='s2x', help='x-component of spin of secondary (for LALSimulation injections only)', default=0.)
  parser.add_option('--s2y', type='float', dest='s2y', help='y-component of spin of secondary (for LALSimulation injections only)', default=0.)
  parser.add_option('--s2z', type='float', dest='s2z', help='e-component of spin of secondary (for LALSimulation injections only)', default=0.)
  parser.add_option('--approx', type='string', dest='approx', help='approximant (for LALSimulation injections only; default = "SEOBNRv3")', default='SEOBNRv3')
  parser.add_option('--dist', dest='dist', type='float', help='distance in Mpc', default=420.)
  parser.add_option('--incl', dest='incl', type='float', help='inclination', default=None)
  parser.add_option('--ra', dest='ra', type='float', help='right ascension', default=None)
  parser.add_option('--dec', dest='dec', type='float', help='declination', default=None)
  parser.add_option('--psi', dest='psi', type='float', help='polarization', default=None)
  parser.add_option('--phi-ref', dest='phi_ref', type='float', help='phase at reference frequency', default=0.)
  parser.add_option('--f-ref', dest='f_ref', type='float', help='reference frequency', default=100.)
  parser.add_option('--f-low', dest='f_low', type='float', help='lower frequency cut-off in Hz (default = 20 Hz)', default=20.)
  parser.add_option('--f-start', dest='f_start', type='string', help='lower frequency for SNR computation (default=f_low)', default=None)
  parser.add_option('--f-noise', dest='f_noise', type='string', help='lower frequency cut-off for noise injection in Hz (default=f_start)', default=None)
  parser.add_option('--fix-snr', dest='inj_snr', type='float', help='SNR (if provided, overrides distance)', default=None)
  parser.add_option('--srate', dest='srate', type='float', help='sampling rate in Hz (default = 4096 Hz)', default=4096.)
  parser.add_option('--seglen', dest='seglen', type='int', help='minimum length of frame to be injected (increased up to the required power of 2)', default=None)
  parser.add_option('--trig-time', dest='trig_time', type='float', help='trigger time', default=None)
  parser.add_option('--ifos', dest='ifos', type='string', help='list of detectors [COMMA DELIMITED NO SPACES] (default = H1,L1)', default='H1,L1')
  parser.add_option('-n', '--noise-models', dest='noise_models', type='string', help='noise models from lalinspiral.sbank.psds / lalsimulation: *either* a list in the same order as the detectors [COMMA DELIMITED NO SPACES] *or* a single noise model for all three detectors *or* O1, O2, design (default=design)', default='design')
  parser.add_option('--data-seed', dest='data_seed', type='int', help='random number seed for noise', default=None)
  parser.add_option('--no-noise', action='store_false', dest='add_noise', help='turn off addition of noise', default=True)
  parser.add_option('-c', '--channel-suffix', dest='channel_suffix', type='string', help='suffix for channel name', default='NR_INJECTED')
  parser.add_option('-p', '--frame-prefix', dest='frame_prefix', type='string', help='prefix for frame name (data file name will be used as default)', default=None)
  parser.add_option('-o', '--out-folder', dest='out_folder', help='output folder to dump frame files', default='../frames')
  parser.add_option('--disable-plot', action='store_false', dest='plotting', help='turn off plotting', default=True)
  parser.add_option('--quiet', action='store_false', dest='verbose', help='turn off verbose', default=True)
  
  (options, args) = parser.parse_args()
  
  data_file = options.data_file
  data_format = options.data_format
  N = options.N
  lmax = options.lmax
  absmmin = options.absmmin
  M = options.M
  q = options.q
  s1x = options.s1x
  s1y = options.s1y
  s1z = options.s1z
  s2x = options.s2x
  s2y = options.s2y
  s2z = options.s2z
  approx = options.approx
  dist = options.dist
  incl = options.incl
  ra = options.ra
  dec = options.dec
  psi = options.psi
  phi_ref = options.phi_ref
  f_ref = options.f_ref
  f_low = options.f_low
  f_start = options.f_start
  f_noise = options.f_noise
  srate = options.srate
  seglen = options.seglen
  trig_time = options.trig_time
  add_noise = options.add_noise
  ifos = options.ifos
  noise_models = options.noise_models
  data_seed = options.data_seed
  channel_suffix = options.channel_suffix
  frame_prefix = options.frame_prefix
  inj_snr = options.inj_snr
  out_folder = options.out_folder
  plotting = options.plotting
  verbose = options.verbose
  
  # Print the command line for retaining in the log
  sys.stdout.write(verbose*(' '.join(sys.argv)+'\n\n'))
  
  # O1 run time: 2015-09-12 to 2016-01-19
  if trig_time is None:
    trig_time = random.uniform(1126051217.0, 1137283217.0)
    sys.stdout.write(verbose*('Randomly choosing trig-time = %f\n'%(trig_time)))
  
  if incl is None:
    incl = np.arccos(random.uniform(-1., 1.))
    sys.stdout.write(verbose*('Randomly choosing incl = %f\n'%(incl)))
  
  if ra is None:
    ra = random.uniform(0., 2.*np.pi)
    sys.stdout.write(verbose*('Randomly choosing ra = %f\n'%(ra)))

  if dec is None:
    dec = np.arcsin(random.uniform(-1., 1.))
    sys.stdout.write(verbose*('Randomly choosing dec = %f\n'%(dec)))

  if psi is None:
    psi = random.uniform(0., 2.*np.pi)
    sys.stdout.write(verbose*('Randomly choosing psi = %f\n'%(psi)))
  
  # NOTE: data_seed is only for the noise generation. So initialize the seed within noiseutils.
  
  if data_seed is None:
    data_seed = nu.random.randint(32768)
    sys.stdout.write(verbose*('Randomly choosing data-seed = %d\n'%(data_seed)))
  
  random.seed(data_seed)
  
  det_list = ifos.split(',')
  det_str = ifos.replace(',', '')
  
  if noise_models=='O1':
    noise_map = dict(zip(det_list, [{'H1': 'aLIGO_EARLY_HIGH', 'L1': 'aLIGO_EARLY_HIGH'}[det] for det in det_list]))
  elif noise_models=='O2':
    noise_map = dict(zip(det_list, [{'H1': 'aLIGO_MID_HIGH', 'L1': 'aLIGO_MID_HIGH', 'V1': 'AdV_EARLY_HIGH'}[det] for det in det_list]))
  elif noise_models=='design':
    noise_map = dict(zip(det_list, [{'H1': 'aLIGO_DESIGN', 'L1': 'aLIGO_DESIGN', 'V1': 'AdV_DESIGN', 'I1': 'aLIGO_DESIGN', 'J1': 'aLIGO_DESIGN', 'K1': 'aLIGO_DESIGN'}[det] for det in det_list]))
  elif len(noise_models.split(','))==len(det_list):
    noise_map = dict(zip(det_list, noise_models.split(',')))
  elif len(noise_models.split(','))==1:
    noise_map = dict(zip(det_list, [noise_models]*len(det_list)))
  else:
    sys.stderr.write('Cannot parse noise_models. Exiting.\n')
    sys.exit()
  
  if f_start is None:
    f_start = str(f_low) # Required for the next step
    f_start_map = dict(zip(det_list, [f_low]*len(det_list)))
  elif len(f_start.split(','))==1:
    f_start_map = dict(zip(det_list, [float(f_start)]*len(det_list)))
  elif len(f_start.split(','))==len(det_list):
    f_start_map = dict(zip(det_list, np.vectorize(float)(f_start.split(','))))
  else:
    sys.stderr.write('Cannot parse f_start. Exiting.\n')
    sys.exit()
  
  if f_noise is None:
    f_noise_map = f_start_map
  elif len(f_noise.split(','))==1:
    f_noise_map = dict(zip(det_list, [float(f_noise)]*len(det_list)))
  elif len(f_noise.split(','))==len(det_list):
    f_noise_map = dict(zip(det_list, np.vectorize(float)(f_noise.split(','))))
  else:
    sys.stderr.write('Cannot parse f_noise. Exiting.\n')
    sys.exit()
  
  if frame_prefix is None:
    if data_format=='LALSim':
      frame_prefix = ('%s_q%s_s1x%s_s1y%s_s1z%s_s2x%s_s2y%s_s2z%s'%(approx, str(round(q,2)), str(round(s1x,2)), str(round(s1y,2)), str(round(s1z,2)), str(round(s2x,2)), str(round(s2y,2)), str(round(s2z,2)))).replace('.', 'p')
    else:
      frame_prefix = os.path.basename(data_file).replace('.gz', '').replace('.dat', '').replace('.h5', '').replace('.hdf5', '')
  
  os.system('mkdir -p %s'%(out_folder))
  
  if data_format=='LALSim':
    if q > 1.:
      q = 1./q
    m1 = M/(1.+q)
    m2 = M*q/(1.+q)
    (t_arr, hp_arr, hc_arr) = wf.data_from_TD(wf.TDwaveform(m1=m1, m2=m2, s1x=s1x, s1y=s1y, s1z=s1z, s2x=s2x, s2y=s2y, s2z=s2z, dist=dist, incl=incl, approx=approx, srate=srate, f_low=f_low, phi_ref=phi_ref, f_ref=f_ref, taper=True))
  else:
    # Load the data
    NR_wf_obj = NR_waveform(data_file, data_format=data_format, N=N, verbose=verbose)
    # Create the waveform in the source frame
    (t_arr, hp_arr, hc_arr) = wf.data_from_TD(NR_wf_obj.TDwaveform(M, dist=dist, incl=incl, f_low=f_low, srate=srate, phi_ref=phi_ref, f_ref=f_ref, taper=True, trim_right=True, trim_tolerance=1e-2, lmax=lmax, absmmin=absmmin))
  
  # Make sure that the coalescence is at least 2s before the end of the segment and choose segment length to minimum required power of 2
  wf_len = t_arr[-1] - t_arr[0]
  gps_end = int(math.ceil(trig_time+2))
  required_seglen = wf_len + gps_end - trig_time + 6384*1e3/C_SI # The last term is the radius of the earth in seconds
  if seglen is None:
    seglen = int(2**math.ceil(np.log2(required_seglen)))
  else:
    seglen = int(2**math.ceil(np.log2(max(seglen, required_seglen))))
  sys.stdout.write(verbose*('Injecting with seglen = %d\n\n'%(seglen)))
  gps_start = gps_end - seglen
  
  # Calculate the SNR
  snrsq_opt = 0.
  sys.stdout.write(verbose*('Preliminary SNR calculation:\n'))
  for det in det_list:
    noise_model = noise_map[det]
    f_start_det = f_start_map[det]
    (t_arr_padded, strain_t) = wf.detector_strain(det, ra, dec, psi, trig_time, t_arr, hp_arr, hc_arr, seglen, gps_end)
    (ff, strain_f) = nu.fd_from_td(t_arr_padded, strain_t)
    vars()[det+'_snr_opt'] = nu.SNR_from_fd(ff, strain_f, noise_model=noise_model, f_low=f_start_det)
    sys.stdout.write(verbose*('%s\toptimal SNR = %.2f\t with %s\tf_start=%.2f\n'%(det, vars()[det+'_snr_opt'], noise_model, f_start_det)))
    snrsq_opt += (vars()[det+'_snr_opt'])**2.
  snr_opt = np.sqrt(snrsq_opt)
  sys.stdout.write(verbose*('Network optimal SNR = %.2f\n\n'%(snr_opt)))
  
  # If injection SNR is provided, scale distance to get appropriate SNR
  if inj_snr is not None:
    inj_snr = float(inj_snr)
    hp_arr *= inj_snr / snr_opt
    hc_arr *= inj_snr / snr_opt
    dist *= snr_opt / inj_snr
    sys.stdout.write(verbose*('Scaling to SNR = %.2f\t:\tdist = %f\n\n'%(inj_snr, dist)))
  
  frame_extn = '%s_M%d_dist%d_incl%s_ra%s_dec%s_psi%s_flow%d'%(frame_prefix, round(M), round(dist), str(round(incl,2)).replace('.', 'p').replace('-', 'm'), str(round(ra,2)).replace('.', 'p').replace('-', 'm'), str(round(dec,2)).replace('.', 'p').replace('-', 'm'), str(round(psi,2)).replace('.', 'p').replace('-', 'm'), round(f_low))
    
  # Inject
  snrsq_opt = 0.
  for det in det_list:
    noise_model = noise_map[det]
    f_start_det = f_start_map[det]
    channel = '%s:%s'%(det, channel_suffix)
    
    (t_arr_padded, strain_t) = wf.detector_strain(det, ra, dec, psi, trig_time, t_arr, hp_arr, hc_arr, seglen, gps_end)
    
    # Create an empty strain vector
    strainT = lal.CreateREAL8TimeSeries(channel, gps_start, 0.0, 1.0/srate, lal.DimensionlessUnit, int(seglen*srate))
    
    strainT.data.data = strain_t
    
    # Retain signal for plotting if necessary
    vars()[det+'_signal'] = copy.deepcopy(strainT.data.data)
    
    # Calculate optimal SNR
    strainF = wf.FD_from_TD(strainT)
    (ff, strain_f_data) = wf.data_from_FD(strainF)
    
    vars()[det+'_snr_opt'] = nu.SNR_from_fd(ff, strain_f_data, noise_model=noise_model, f_low=f_start_det)
    
    sys.stdout.write(verbose*('%s\toptimal SNR = %.2f\t with %s\tf_start=%.2f\n'%(det, vars()[det+'_snr_opt'], noise_model, f_start_det)))
    
    snrsq_opt += (vars()[det+'_snr_opt'])**2.
    
    # Add noise
    if add_noise:
      f_noise_det = f_noise_map[det]
      noise_TD = nu.TDnoise(strainT.data.length, strainT.deltaT, noise_model=noise_model, f_low=f_noise_det)
      strainT.data.data += noise_TD
      noise_str = '%s_seed%d_%s'%(noise_model, data_seed, det_str)
      sys.stdout.write(verbose*('For detector %s added %s noise with lower cut-off f_noise = %.2f\n'%(det, noise_model, f_noise_det)))
    else:
      noise_str = 'nonoise'
    
    # Write the strain to a frame file
    FrWriteREAL8TimeSeries(strainT, 0)
    
    default_filename = '%s-%s-%d-%d.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen)
    if glob.glob(default_filename)==[]:
      sys.stderr.write('Cannot find %s. Exiting.\n'%(default_filename))
      sys.exit()
      #default_filename = '%s-%s-%d-%d.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen+1)
    
    frame_filename = '%s-%s-%d-%d_%s_%s.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen, frame_extn, noise_str)
    
    os.system('mkdir -p %s'%(out_folder))
    target_filename = os.path.join(out_folder, frame_filename)
    os.system('mv %s %s'%(default_filename, target_filename))
    sys.stdout.write(verbose*('Moved %s to %s\n'%(default_filename, os.path.realpath(target_filename))))
  
  snr_opt = np.sqrt(snrsq_opt)
  sys.stdout.write(verbose*('Injected with network optimal SNR = %.2f\n'%(snr_opt)))
  
  # Write injection parameters and gps time file
  
  if add_noise:
    noise_str = 'seed%d_%s'%(data_seed, det_str)
  else:
    noise_str = 'nonoise'
  
  injparams_filename = '%s-%d-%d_%s_%s_injparams.dat'%(channel_suffix, int(gps_start), seglen, frame_extn, noise_str)
  
  gpstime_filename = '%s-%d-%d_%s_%s_gpstime.txt'%(channel_suffix, int(gps_start), seglen, frame_extn, noise_str)
  
  np.savetxt(os.path.join(out_folder, injparams_filename), np.array([[trig_time, M, dist, incl, ra, dec, psi, f_low, f_ref, phi_ref]]), header='trig_time M dist incl ra dec psi f_low f_ref phi_ref', fmt='%f')
  
  np.savetxt(os.path.join(out_folder, gpstime_filename), np.array([trig_time]), fmt='%f')
  
  # Optional plotting
  if plotting:
    from pylal import Fr
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    
    noise_colors = {'H1':'salmon', 'L1':'mediumaquamarine', 'V1':'paleturquoise', 'I1':'tan', 'J1':'lightgrey', 'K1':'lightgrey'}
    
    signal_colors = {'H1':'darkred', 'L1':'darkolivegreen', 'V1':'navy', 'I1':'saddlebrown', 'J1':'black', 'K1':'black'}
    
    # Time-domain plot
    plt.figure()
    
    for det in det_list:
      noise_model = noise_map[det]
      channel = '%s:%s'%(det, channel_suffix)
      if add_noise:
        noise_str = '%s_seed%d_%s'%(noise_model, data_seed, det_str)
      else:
        noise_str = 'nonoise'
      frame_filename = '%s-%s-%d-%d_%s_%s.gwf'%(det[0], channel.replace(':', '_'), int(gps_start), seglen, frame_extn, noise_str)
      
      frdata, gps_start, xoffset, xspacing, xunit, yunit = Fr.frgetvect(os.path.join(out_folder, frame_filename), channel)
      
      xx = np.arange(0., len(frdata))*xspacing
      plt.plot(xx, frdata, colors.cnames[noise_colors[det]], label=det+': data')
    
    if add_noise:
      for det in det_list:
        plt.plot(xx, vars()[det+'_signal'], colors.cnames[signal_colors[det]], label=det+': signal')
      noise_str = 'seed%d_%s'%(data_seed, det_str)
    else:
      noise_str = 'nonoise_%s'%(det_str)

    plt.ylim(-4e-21, 4e-21)
    plt.legend(loc='upper left',title='%.4f'%(trig_time))
    
    plt.xlabel(r'Time, $t$ [s] after GPS time of %ds with srate = %dHz'%(gps_start, srate))
    plt.ylabel(r'$h(t)$')
    
    plot_filename = '%s-%d-%d_%s_%s.png'%(channel_suffix, int(gps_start), seglen, frame_extn, noise_str)
    
    plt.savefig(os.path.join(out_folder, plot_filename), dpi=300)
    sys.stdout.write(verbose*('Saved: %s\n\n'%(plot_filename)))
    
    # Frequency-domain plot
    plt.figure()
    
    for det in det_list:
      noise_model = noise_map[det]
      f_start_det = f_start_map[det]
      f_stop_det = srate/2.
      f_arr = np.linspace(f_start_det, f_stop_det, 200)
      
      if add_noise:
        plt.loglog(f_arr, np.sqrt(nu.LIGOPSD(f_arr, noise_model=noise_model)), colors.cnames[noise_colors[det]], label=det+': noise ASD')
      
      (f_arr, strain_f) = nu.fd_from_td(xx, vars()[det+'_signal'])
      plt.loglog(f_arr, 2.*np.sqrt(f_arr)*np.abs(strain_f), colors.cnames[signal_colors[det]], label=det+': signal')
      
    plt.axvline(f_low, color='k')

    plt.xlim(5, 4000)
    plt.ylim(1e-24, 1e-21)
    plt.legend(loc='lower left',title='%.4f'%(trig_time))
    
    plt.xlabel(r'Frequency, $f$ [Hz]')
    if add_noise:
      noise_str = 'seed%d_%s'%(data_seed, det_str)
      plt.ylabel(r'$\sqrt{S(f)}$  and  $2|h(f)|\sqrt{f}$    [Hz$^{-1/2}$]')
    else:
      noise_str = 'nonoise_%s'%(det_str)
      plt.ylabel(r'$2|h(f)|\sqrt{f}$ [Hz$^{-1/2}$]')
    
    plot_filename = '%s-%d-%d_%s_%s_FD.png'%(channel_suffix, int(gps_start), seglen, frame_extn, noise_str)
    
    plt.savefig(os.path.join(out_folder, plot_filename), dpi=300)
    sys.stdout.write(verbose*('Saved: %s\n\n'%(plot_filename)))
