#!/usr/bin/python
"""
Run lalinference_nest for example injection

(C) P. Ajith, Archisman Ghosh, Abhirup Ghosh; Modified: 2015-01-12
"""

import os
import optparse as op
import lal
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import lalsimulation as lalsim



# Executables
lalinferencenest = 'lalinference_nest'
lalappsnest2pos = 'lalapps_nest2pos'
cbcBayesPostProc = 'cbcBayesPostProc.py'

if __name__=="__main__":
  parser = op.OptionParser()
  parser.add_option('-M', type="float", dest="M", help="total mass")
  parser.add_option('-q', type="float", dest="q", help="mass ratio")
  parser.add_option('-n', '--n_isco', type="float", dest="n_isco", help="fraction of f_isco")
  (options, args) = parser.parse_args()
  (M, q, n_isco) = (options.M, options.q, options.n_isco)
  
  # Injection arguments
  inj_xml_file = '../injections/example_inj_IMRPhenomB_snr25_%s_%s.xml'%(str(M), str(q)) 
  #inj_xml_file = '../injections/example_inj_SEOBNRv2_ROM_snr_%s_%s.xml'%(str(M), str(q))
  event_id = 0 # event id
  dataseed = 42 # seed for generating noise realization
  injection_args = '--inj %s --event %d --dataseed %d'%(inj_xml_file, event_id, dataseed)

  # Recovery waveform arguments
  rec_waveform = 'IMRPhenomBpseudoFourPN' #'EOBNRv2threePointFivePN' # approximant
  #rec_waveform = 'SEOBNRv2_ROM_SingleSpinthreePointFivePN' #'TaylorF2threePointFivePN' # approximant
  modeldomain = 'frequency' # time or frequency
  randomseed = 73 # seed for MCMC
  approx_args = '--approx %s --modeldomain="%s" --amporder 0 --randomseed %d'%(rec_waveform, modeldomain, randomseed)

  # Time args
  trigtime = 1041033614 # gps time of injection used as trigger
  psd_start_time = trigtime - 16. # gps start time of psd estimation data
  time_interv = 0.1 # width of time prior centred around trigger (s)
  time_args = '--trigtime %f --dt %f --psdstart %f'%(trigtime, time_interv, psd_start_time)

  # PSD estimation arguments
  psdlength = 1024 # length of psd estimation (s)
  seglen = 32 # length of segments in psd estimation and analysis (s)
  srate = 2048 # sampling rate (Hz)
  psd_args = '--psdlength %d --seglen %d --srate %d'%(psdlength, seglen, srate)

  # Interferometers, time, PSD estimation, nested sampling arguments
  f_low = n_isco*(1./6.)**(3./2)/(LAL_PI*M*LAL_MTSUN_SI) # lower freqency cut-off for overlap integral
  if q == 1.:
	f_qnm = 0.0879
  elif q == 2.:
	f_qnm = 0.0829
  elif q == 4.:
	f_qnm = 0.0744
  elif q == 10.:
	f_qnm = 0.0661
  #f_low = n_isco*f_qnm/(M*LAL_MTSUN_SI)
  f_high = 0.5*srate-1.
  #ifo_args = '--ifo H1 --H1-cache LALSimAdLIGO --H1-channel LALSimAdLIGO --H1-timeslide 0 --H1-flow %f --H1-fhigh %f'%(f_low, f_high)
  ifo_args = '--ifo H1 --H1-cache LALSimAdLIGO --H1-channel LALSimAdLIGO --H1-timeslide 0 --H1-flow %f --H1-fhigh %f --0noise'%(f_low, f_high)

  # fix the injection parameters 
  m1_inj = M/(1.+q)    # mass 1 (M_sun)  
  m2_inj = (M*q)/(1.+q)    # mass 2 (M_sun) 
  
  # Priors on parameters
  #d_min = 10.; d_max = 1500. # distance (Mpc)
  comp_mass_min = 2.; comp_mass_max = 2.0*m2_inj # component mass (M_sun)
  tot_mass_min = 10.; tot_mass_max = 2.0*M # total mass (M_sun)
  #prior_args = '--distance-min %f --distance-max %f --comp-min %f --comp-max %f --mtotal-min %f --mtotal-max %f'%(d_min, d_max, comp_mass_min, comp_mass_max, tot_mass_min, tot_mass_max)
  prior_args = '--mtotal-min %f --mtotal-max %f --comp-min %f --comp-max %f'%(tot_mass_min, tot_mass_max, comp_mass_min, comp_mass_max)

  
  # Fixed parameters
  fix_params = '--disable-spin --margphi --fix-phase 0 --fix-inclination 0 --fix-rightascension 0 --fix-declination 0 --fix-polarisation 0 --pinparams [distance]'
  #fix_params = '--disable-spin --margphi --pinparams [time,ra,dec,psi,theta_jn,dist]'

  # Nested sampling arguments
  nlive = 16*1024
  ns_args = '--nlive %d --resume --progress'%(nlive)

  # Output arguments
  out_dir = '/scratch/runs/abhirup/2015-06-05/ringdown_snr25_flow_fisco/%s_%s/nisco_%s'%(str(M), str(q), str(n_isco))
  dat_file = '%s/engine/lalinferencenest_output.dat'%(out_dir)
  pos_file = '%s/posterior_samples/posterior_samples.dat'%(out_dir)
  proc_dir = '%s/cbc_bayes'%(out_dir)
  output_args = '--outfile %s'%(dat_file)
  os.system('mkdir -p %s %s %s'%(os.path.dirname(dat_file), os.path.dirname(pos_file), proc_dir))

  # Run lalinference_nest
  inf_command = '%s %s %s %s %s %s %s %s %s %s'%(lalinferencenest, injection_args, approx_args, ifo_args, time_args, psd_args, prior_args, fix_params, ns_args, output_args)
  n2p_command = '%s %s --headers %s_params.txt --Nlive %d --pos %s'%(lalappsnest2pos, dat_file, dat_file, nlive, pos_file)
  #n2p_command = '%s %s --headers %s_params.txt --Nlive %d --pos %s --npos=10000'%(lalappsnest2pos, dat_file, dat_file, nlive, pos_file)
  bpp_command = '%s %s --bsn %s_B.txt --skyres 1.0 --inj %s --eventnum %d -o %s'%(cbcBayesPostProc, pos_file, pos_file, inj_xml_file, event_id, proc_dir)

  os.system('cp %s %s'%(__file__,out_dir))
  print(inf_command)
  os.system(inf_command)
  print(n2p_command)
  os.system(n2p_command)
  print(bpp_command)
  os.system(bpp_command)
