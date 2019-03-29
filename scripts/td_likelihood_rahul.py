import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os, lal, commands, glob, sys
import emcee
import optparse as op
import corner
import nr_fits as nr
from lal import MTSUN_SI as LAL_MTSUN_SI
import bayesian as ba
import scipy
import time
import scipy
from scipy.interpolate import interp1d
from pylal import antenna
import td_likelihood_utils as uc
import standard_gwtransf as gw
import imrtgrutils_final as tgr
from pylal import git_version
import noiseutils as nu

"""
This code assumes the LALInference convention of m1>m2 and q=m1/m2<1.
"""
def lnlike_lalsim(theta, t, d):
    mc, q, a1z, a2z = theta
    m1, m2 = gw.comp_from_mcq(mc, q)
    try:
        t_h, h = uc.signalTD(approx, m1, m2, a1z, a2z, dL, flow, srate, gpsTime, ra, dec, iota, psi, 'radians', detector)
        if asd_file is not None:
        	h = uc.whiten(h, psd_interp_obj, delta_t)
	h_interp = scipy.interpolate.interp1d(t_h, h, fill_value=0., bounds_error=False)
        h = h_interp(t)
    except ValueError:
        return -np.inf
    return -0.5*np.dot(d-h, d-h)/(sigma*sigma)

def lnprior_lalsim(theta, mc_min, mc_max, q_min, q_max, a_min, a_max):
    mc, q, a1z, a2z = theta
    m1, m2 = gw.comp_from_mcq(mc, q)
    if mc_min < mc < mc_max and q_min < q < q_max and a_min < a1z < a_max and a_min < a2z < a_max:
        return 0.
    return -np.inf

def lnprob_lalsim(theta, t, data):
    lp_lalsim = lnprior_lalsim(theta, mc_min, mc_max, q_min, q_max, a_min, a_max)
    if not np.isfinite(lp_lalsim):
        return -np.inf
    return lp_lalsim + lnlike_lalsim(theta, t, data)

def lnlike_ring(theta, t, d):
    A, tau, f_qnm, phi0 = theta
    h = uc.qnm(t, t0_ring, A, tau, f_qnm, phi0)
    return -0.5*np.dot(d-h, d-h)/(sigma*sigma)

def lnprior_ring(theta, A_min, A_max, tau_min, tau_max, f_qnm_min, f_qnm_max, phi0_min, phi0_max):
    A, tau, f_qnm, phi0 = theta
    if A_min < A < A_max and tau_min < tau < tau_max and f_qnm_min < f_qnm <f_qnm_max and phi0_min < phi0 < phi0_max:
        return 0.0
    return -np.inf

def lnprob_ring(theta, t, data):
    lp_ring = lnprior_ring(theta, A_min, A_max, tau_min, tau_max, f_qnm_min, f_qnm_max, phi0_min, phi0_max)
    if not np.isfinite(lp_ring):
        return -np.inf
    return lp_ring + lnlike_ring(theta, t, data)

####################################################################################

if __name__ == '__main__':
  start_time = time.time()

  parser = op.OptionParser()
  parser.add_option("--asd-file", dest="asd_file", help="asd file, required only if working with gaussian white noise", default=None)
  parser.add_option("--white-noise-sigma", dest="sigma", type=float, help="standard deviation for white noise injections", default=1.)
  parser.add_option("--num-threads", dest="num_threads", type=int, help="number of threads to be used", default=1)
  parser.add_option("--m1-inj", dest="m1_inj", type=float, help="injected value of component mass m1")
  parser.add_option("--m2-inj", dest="m2_inj", type=float, help="injected value of component mass m2")
  parser.add_option("--a1z-inj", dest="a1z_inj", type=float, help="injected value of z-component of spin of mass m1")
  parser.add_option("--a2z-inj", dest="a2z_inj", type=float, help="injected value of z-component of spin of mass m2")
  parser.add_option("--dL-inj", dest="dL", type=float, help="injected value of luminosity distance")
  parser.add_option("--ra-inj", dest="ra", type=float, help="injected value of RA", default=0.)
  parser.add_option("--dec-inj", dest="dec", type=float, help="injected value of dec", default=0.)
  parser.add_option("--iota-inj", dest="iota", type=float, help="injected value of iota", default=0.)
  parser.add_option("--psi-inj", dest="psi", type=float, help="injected value of polarisation", default=0.)
  parser.add_option("--approximant", dest="approx", help="approximant used for data generation and recovery")
  parser.add_option("-f", "--fit-formula", dest="fit_formula", help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_Pan2011', 'nonprecspin_Healy2014', 'bbh_average_fits_precessing'", default="nonprecspin_Healy2014")
  parser.add_option("--sample-rate", dest="srate", type=float, help="sampling rate for waveform generation", default=2048.)
  parser.add_option("--flow", dest="flow", type=float, help="starting frequency for waveform generation", default=10.)
  parser.add_option("--detector", dest="detector", help="detector")
  parser.add_option("--gpsTime", dest="gpsTime", help="gpsTime of trigger")
  parser.add_option("--outdir", dest="outdir", help="output directory")
  parser.add_option("--ringdown_PE", dest="ringdown_pe", help="PE channel to follow. Options: ['lalsimulation', 'QNM'")
  parser.add_option("--t0-ring", type=float, dest="t0_ring", help="time after the peak waveform used for QNM analysis")
  parser.add_option("--chirpmass-min", type=float, dest="mc_min", help="minimum of chirp mass prior")
  parser.add_option("--chirpmass-max", type=float, dest="mc_max", help="maximum of chirp mass prior")
  parser.add_option("--q-min", type=float, dest="q_min", help="minimum of q prior")
  parser.add_option("--q-max", type=float, dest="q_max", help="maximum of q prior")
  parser.add_option("--az_spin-min", type=float, dest="a_min", help="minimum of z-component spin prior")
  parser.add_option("--az_spin-max", type=float, dest="a_max", help="maximum of z-component spin prior")
  parser.add_option("--A-min", type=float, dest="A_min", help="minimum of A prior")
  parser.add_option("--A-max", type=float, dest="A_max", help="maximum of A prior")
  parser.add_option("--tau-min", type=float, dest="tau_min", help="minimum of tau prior")
  parser.add_option("--tau-max", type=float, dest="tau_max", help="maximum of tau prior")
  parser.add_option("--f_qnm-min", type=float, dest="f_qnm_min", help="minimum of f_qnm prior")
  parser.add_option("--f_qnm-max", type=float, dest="f_qnm_max", help="maximum of f_qnm prior")
  parser.add_option("--phi0-min", type=float, dest="phi0_min", help="minimum of phi0 prior")
  parser.add_option("--phi0-max", type=float, dest="phi0_max", help="maximum of phi0 prior")
  parser.add_option("--insp-mcmc", type=int, dest="n_steps_mcmc_i", help="no. of steps for inspiral MCMC")
  parser.add_option("--insp-burnin", type=int, dest="n_steps_burnin_i", help="no. of steps for inspiral burn-in")
  parser.add_option("--ring-mcmc", type=int, dest="n_steps_mcmc_r", help="no. of steps for ringdown MCMC")
  parser.add_option("--ring-burnin", type=int, dest="n_steps_burnin_r", help="no. of steps for ringdown burn-in")
  (options, args) = parser.parse_args()
  asd_file = options.asd_file
  sigma = options.sigma
  num_threads = options.num_threads
  m1_inj = options.m1_inj
  m2_inj = options.m2_inj
  a1z_inj = options.a1z_inj
  a2z_inj = options.a2z_inj
  dL = options.dL
  ra = options.ra
  dec = options.dec
  iota = options.iota
  psi = options.psi
  approx = options.approx
  fit_formula = options.fit_formula
  srate = options.srate
  flow = options.flow
  detector = options.detector
  gpsTime = options.gpsTime
  outdir = options.outdir
  ringdown_pe = options.ringdown_pe
  mc_min = options.mc_min
  mc_max = options.mc_max
  q_min = options.q_min
  q_max = options.q_max
  a_min = options.a_min
  a_max = options.a_max
  if ringdown_pe == 'QNM':
    A_min = options.A_min
    A_max = options.A_max
    tau_min = options.tau_min
    tau_max = options.tau_max
    f_qnm_min = options.f_qnm_min
    f_qnm_max = options.f_qnm_max
    phi0_min = options.phi0_min
    phi0_max = options.phi0_max
  t0_ring = options.t0_ring
  n_steps_mcmc_i = options.n_steps_mcmc_i
  n_steps_burnin_i = options.n_steps_burnin_i
  n_steps_mcmc_r = options.n_steps_mcmc_r
  n_steps_burnin_r = options.n_steps_burnin_r

  #################################################################################################
  # output directory
  #################################################################################################
  # create output directory and copy current script to data folder
  os.system('mkdir -p %s/data %s/img'%(outdir, outdir))
  os.system('cp %s %s/data' %(__file__, outdir))
  os.system('cp ./run_td_likelihood_pipeline.py %s/data' %(outdir))

  # creating file to save the run command
  run_command = open('%s/command.txt'%(outdir),'w')
  for arg in sys.argv:
    run_command.write('%s\t' %arg)
  run_command.write("\n")
  run_command.write("\n")
  run_command.write("%s"%git_version.verbose_msg)
  run_command.close()

  #################################################################################################
  # Generating Signal
  #################################################################################################
  t, h = uc.signalTD(approx, m1_inj, m2_inj, a1z_inj, a2z_inj, dL, flow, srate, gpsTime, ra, dec, iota, psi, 'radians', detector)
  t_trig = t[np.where(abs(h) == max(abs(h)))]
  delta_t = 1./srate
  delta_f = 1. / len(h) / delta_t

  #################################################################################################
  # SNR computation and Defining ASD if not white noise
  #################################################################################################
  if asd_file is not None:
    psd = uc.psd_from_txt(asd_file=asd_file, delta_f=delta_f, flow=flow, f_high=srate/2.)
    psd_interp_obj = interp1d(psd.sample_frequencies, psd, bounds_error=False)

  #################################################################################################
  # Generating Noise and computing SNR
  #################################################################################################
  if asd_file is not None:
    # adding noise
    nt = uc.TDnoise(delta_t=delta_t, asd_file=asd_file, flow=flow, f_high=srate/2.)

    snr_opt = uc.SNR_from_td(h, delta_t, psd, flow, srate)
    estimated_psd = uc.PSD_from_TDnoise(nt, delta_t)

    d = h + np.asarray(nt[:len(h)])

    # whitening data
    d = uc.whiten(d, psd_interp_obj, delta_t)
    h = uc.whiten(h, psd_interp_obj, delta_t)

  else:
    d = h + np.random.randn(len(h))*sigma
    f, hf = nu.fd_from_td(t, h)
    df = np.mean(np.diff(f))
    band_idx = (f >= flow) * (f <= srate/2.)
    snr_opt = 2*np.sqrt(np.sum(np.nan_to_num(df*abs(hf[band_idx])**2/(sigma*sigma*np.ones(len(hf[band_idx]))))))

  print '... signal injected with SNR %.2f'%(snr_opt)

  # Splitting data for inspiral and ringdown analysis
  t_i = t[np.where(t < t_trig)]; h_i = h[np.where(t < t_trig)]; d_i = d[np.where(t < t_trig)]
  t_r = t[np.where(t > t_trig + t0_ring)]; h_r = h[np.where(t > t_trig + t0_ring)]; d_r = d[np.where(t > t_trig + t0_ring)]

  print '... data read'

  #################################################################################################
  # Recovery Parameters
  #################################################################################################
  # initial conditions for INSPIRAL emcee sampler
  mc_init, q_init = gw.mcq_from_comp(m1_inj, m2_inj)
  a1z_init, a2z_init = a1z_inj, a2z_inj

  # initial conditions for RINGDOWN emcee sampler
  Mf_inj, af_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, 0., 0., a1z_inj, a2z_inj, 0., 0., 0., fit_formula)
  Omega_init = nr.qnmfreqs_berti(af_inj, 2, 2, 0)/(Mf_inj*LAL_MTSUN_SI)
  f_qnm_init = np.real(Omega_init)/(2*np.pi)
  tau_init = -1./np.imag(Omega_init)
  A_init = max(abs(h))
  phi0_init = 1.

  print '... initial conditions: mc_init=%.2f, q_init=%.2f, a1z_init=%.2f, a2z_init=%.2f, A_init=%e, tau_init=%.5f, f_qnm_init=%.2f, phi0_init=%.2f'%(mc_init, q_init, a1z_init, a2z_init, A_init, tau_init, f_qnm_init, phi0_init)

  #################################################################################################
  # MCMC inspiral lalsim
  #################################################################################################
  ndim, nwalkers = 4, 100
  pos = [[mc_init, q_init, a1z_init, a2z_init] + np.random.rand()*np.array([0.1, 0.01, 0.01, 0.01]) for i in range(nwalkers)];

  sampler_i = emcee.EnsembleSampler(nwalkers, ndim, lnprob_lalsim, args=(t_i, d_i), threads=num_threads);
  sampler_i.run_mcmc(pos, n_steps_mcmc_i);
  mc_i_chain, q_i_chain, a1z_i_chain, a2z_i_chain = sampler_i.chain[:, :, 0].T, sampler_i.chain[:, :, 1].T, sampler_i.chain[:, :, 2].T, sampler_i.chain[:, :, 3].T
  samples_i = sampler_i.chain[:, :, :].reshape((-1, ndim));
  lnprob_i = sampler_i.flatlnprobability

  print '... inspiral mcmc done'

  #################################################################################################
  # MCMC merger-ringdown
  #################################################################################################

  if ringdown_pe == 'lalsimulation':
    ndim, nwalkers = 4, 100
    pos = [[mc_init, q_init, a1z_init, a2z_init] + np.random.rand()*np.array([0.1, 0.01, 0.01, 0.01]) for i in range(nwalkers)];

    sampler_r = emcee.EnsembleSampler(nwalkers, ndim, lnprob_lalsim, args=(t_r, d_r), threads=num_threads);
    sampler_r.run_mcmc(pos, n_steps_mcmc_r);
    mc_r_chain, q_r_chain, a1z_r_chain, a2z_r_chain = sampler_r.chain[:, :, 0].T, sampler_r.chain[:, :, 1].T, sampler_r.chain[:, :, 2].T, sampler_r.chain[:, :, 3].T

  elif ringdown_pe == 'QNM':
    ndim, nwalkers = 4, 100
    pos = [[A_init, tau_init, f_qnm_init, phi0_init] + np.random.rand()*np.array([0.1, 1e-3, 1., 0.1]) for i in range(nwalkers)]; ###FIXME###

    sampler_r = emcee.EnsembleSampler(nwalkers, ndim, lnprob_ring, args=(t_r, d_r), threads=num_threads);
    sampler_r.run_mcmc(pos, n_steps_mcmc_r);
    A_r_chain, tau_r_chain, f_qnm_r_chain, phi0_r_chain = sampler_r.chain[:, :, 0].T, sampler_r.chain[:, :, 1].T, sampler_r.chain[:, :, 2].T, sampler_r.chain[:, :, 3].T
  lnprob_r = sampler_r.flatlnprobability
  samples_r = sampler_r.chain[:, :, :].reshape((-1, ndim));

  print '... merger-ringdown mcmc done'

  #################################################################################################
  # Eliminating burnin steps
  #################################################################################################
  samples_i = samples_i[n_steps_burnin_i:]; lnprob_i = lnprob_i[n_steps_burnin_i:]
  samples_r = samples_r[n_steps_burnin_r:]; lnprob_r = lnprob_r[n_steps_burnin_r:]

  mc_i, q_i, a1z_i, a2z_i = samples_i[:,0], samples_i[:,1], samples_i[:,2], samples_i[:,3]
  m1_i, m2_i = gw.comp_from_mcq(mc_i, q_i)

  plt.scatter(mc_i, q_i, c=lnprob_i)
  plt.savefig('sanity_check.png')

  if ringdown_pe == 'lalsimulation':
    mc_r, q_r, a1z_r, a2z_r = samples_r[:,0], samples_r[:,1], samples_r[:,2], samples_r[:,3]
    m1_r, m2_r = gw.comp_from_mcq(mc_r, q_r)

  elif ringdown_pe == 'QNM':
    A_r, tau_r, f_qnm_r, phi0_r = samples_r[:,0], samples_r[:,1], samples_r[:,2], samples_r[:,3]

  print '... eliminated burn-in steps'

  #################################################################################################
  # Converting to (Mf, af) samples
  #################################################################################################
  # Inspiral
  Mf_i, af_i = tgr.calc_final_mass_spin(m1_i, m2_i, 0., 0., a1z_i, a2z_i, 0., 0., 0., fit_formula)

  # Ringdown
  if ringdown_pe == 'lalsimulation':
    Mf_r, af_r = tgr.calc_final_mass_spin(m1_r, m2_r, 0., 0., a1z_r, a2z_r, 0., 0., 0., fit_formula)

  elif ringdown_pe == 'QNM':
    Mf_r, af_r = uc.mfaf_from_qnmfreqs_berti(f_qnm_r, tau_r, 2, 2, 0)
    Mf_r[np.where(Mf_r < 0.)] = 0.
    af_r[np.where(af_r < 0.)] = 0.


  # Best-fit waveform parameters
  mc_i_med, q_i_med, m1_i_med, m2_i_med, a1z_i_med, a2z_i_med = np.median(mc_i), np.median(q_i), np.median(m1_i), np.median(m2_i), np.median(a1z_i), np.median(a2z_i)
  t_h_i_rec, h_i_rec =  uc.signalTD(approx, m1_i_med, m2_i_med, a1z_i_med, a2z_i_med, dL, flow, srate, gpsTime, ra, dec, iota, psi, 'radians', detector)
  if asd_file is not None:
        h_i_rec = uc.whiten(h_i_rec, psd_interp_obj, delta_t)
  h_i_rec_interp = scipy.interpolate.interp1d(t_h_i_rec, h_i_rec, fill_value=0., bounds_error=False)
  h_i_rec = h_i_rec_interp(t_i)

  if ringdown_pe == 'lalsimulation':
    mc_r_med, q_r_med, m1_r_med, m2_r_med, a1z_r_med, a2z_r_med = np.median(mc_r), np.median(q_r), np.median(m1_r), np.median(m2_r), np.median(a1z_r), np.median(a2z_r)
    t_h_r_rec, h_r_rec =  uc.signalTD(approx, m1_r_med, m2_r_med, a1z_r_med, a2z_r_med, dL, flow, srate, gpsTime, ra, dec, iota, psi, 'radians', detector)
    if asd_file is not None:
        h_r_rec = uc.whiten(h_r_rec, psd_interp_obj, delta_t)
    h_r_rec_interp = scipy.interpolate.interp1d(t_h_r_rec, h_r_rec, fill_value=0., bounds_error=False)
    h_r_rec = h_r_rec_interp(t_r)
  elif ringdown_pe == 'QNM':
    A_r_med, tau_r_med, f_qnm_r_med, phi0_r_med = np.median(A_r), np.median(tau_r), np.median(f_qnm_r), np.median(phi0_r)
    h_r_rec = uc.qnm(t_r, t0_ring, A_r_med, tau_r_med, f_qnm_r_med, phi0_r_med)

  #################################################################################################
  # Saving data
  #################################################################################################
  if asd_file is not None:
    asd_file = glob.glob('/home/abhirup/Documents/Work/testGR_IR/PSD/*%s*.txt'%(asd_file))[0]
    os.system('cp -r %s %s/data'%(asd_file, outdir))
  np.savetxt(outdir+'/data/data.dat', np.c_[t, h, d])
  np.savetxt(outdir+'/data/snr.txt', [snr_opt])
  np.savetxt(outdir+'/data/inspiral_data.dat', np.c_[t_i, h_i, d_i])
  np.savetxt(outdir+'/data/ringdown_data.dat', np.c_[t_r, h_r, d_r])
  np.savetxt(outdir+'/data/inspiral_samples_preburnin.dat', samples_i)
  np.savetxt(outdir+'/data/ringdown_samples_preburnin.dat', samples_r)
  np.savetxt(outdir+'/data/emcee_inits.dat', np.c_[mc_init, q_init, a1z_init, a2z_init, A_init, tau_init, f_qnm_init, phi0_init])
  np.savetxt(outdir+'/data/inspiral_samples_postburnin.dat', np.c_[mc_i, q_i, m1_i, m2_i, a1z_i, a2z_i, Mf_i, af_i, lnprob_i], header='mc q m1 m2 a1 a2 mf af lnprob')
  if ringdown_pe == 'lalsimulation':
    np.savetxt(outdir+'/data/ringdown_samples_postburnin.dat', np.c_[mc_r, q_r, m1_r, m2_r, a1z_r, a2z_r, Mf_r, af_r, lnprob_r], header='mc q m1 m2 a1 a2 mf af lnprob')
  elif ringdown_pe == 'QNM':
    np.savetxt(outdir+'/data/ringdown_samples_postburnin.dat', np.c_[A_r, tau_r, f_qnm_r, phi0_r, Mf_r, af_r, lnprob_r], header='A tau f_qnm phi0 mf af lnprob')

  print '... data saved'

  #################################################################################################
  # Plotting 
  #################################################################################################
  # Inspiral Chain plot
  plt.figure(figsize=(16,8))
  plt.subplot(421)
  plt.plot(mc_i_chain, color="k", alpha=0.4, lw=0.5)
  plt.plot(mc_init + np.std(mc_i_chain, axis=1), 'r')
  plt.axhline(y=mc_init, color='g')
  plt.ylabel('$M_c$')
  plt.subplot(423)
  plt.plot(q_i_chain, color="k", alpha=0.4, lw=0.5)
  plt.plot(q_init + np.std(q_i_chain, axis=1), 'r')
  plt.axhline(y=q_init, color='g')
  plt.ylabel('$q$')
  plt.subplot(425)
  plt.plot(a1z_i_chain, color="k", alpha=0.4, lw=0.5)
  plt.plot(a1z_init + np.std(a1z_i_chain, axis=1), 'r')
  plt.axhline(y=a1z_inj, color='g')
  plt.ylabel('$a_{1,z}$')
  plt.subplot(427)
  plt.plot(a2z_i_chain, color="k", alpha=0.4, lw=0.5)
  plt.plot(a2z_init + np.std(a2z_i_chain, axis=1), 'r')
  plt.axhline(y=a2z_inj, color='g')
  plt.ylabel('$a_{2,z}$')
  if ringdown_pe == 'lalsimulation':
    plt.subplot(422)
    plt.plot(mc_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(mc_init + np.std(mc_r_chain, axis=1), 'r')
    plt.axhline(y=mc_init, color='g')
    plt.ylabel('$M_c$')
    plt.subplot(424)
    plt.plot(q_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(q_init + np.std(q_r_chain, axis=1), 'r')
    plt.axhline(y=q_init, color='g')
    plt.ylabel('$q$')
    plt.subplot(426)
    plt.plot(a1z_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(a1z_init + np.std(a1z_r_chain, axis=1), 'r')
    plt.axhline(y=a1z_inj, color='g')
    plt.ylabel('$a_{1,z}$')
    plt.subplot(428)
    plt.plot(a2z_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(a2z_init + np.std(a2z_r_chain, axis=1), 'r')
    plt.axhline(y=a2z_inj, color='g')
    plt.ylabel('$a_{2,z}$')

  elif ringdown_pe == 'QNM':
    plt.subplot(422)
    plt.plot(A_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(A_init + np.std(A_r_chain, axis=1), 'r')
    plt.ylabel('$A$')
    plt.subplot(424)
    plt.plot(tau_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(tau_init + np.std(tau_r_chain, axis=1), 'r')
    plt.ylabel('$\\tau$')
    plt.subplot(426)
    plt.plot(f_qnm_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(f_qnm_init + np.std(f_qnm_r_chain, axis=1), 'r')
    plt.ylabel('$f_{qnm}$')
    plt.subplot(428)
    plt.plot(phi0_r_chain, color="k", alpha=0.4, lw=0.5)
    plt.plot(phi0_init + np.std(phi0_r_chain, axis=1), 'r')
    plt.ylabel('$\phi_0$')
  plt.savefig('%s/img/PE_samples_chain.png'%outdir, dpi=300)

  # Inspiral corner plot
  corner.corner(samples_i, labels=["$M_c (M_{\odot})$", "$q$", "$a_{1,z}$", "$a_{2,z}$"], truths=[mc_i_med, q_i_med, a1z_i_med, a2z_i_med])
  plt.savefig('%s/img/PE_inspiral_corner_postburnin.png'%outdir, dpi=300)

  # Merger-ringdown corner plot
  if ringdown_pe == 'lalsimulation':
    corner.corner(samples_r, labels=["$M_c (M_{\odot})$", "$q$", "$a_{1,z}$", "$a_{2,z}$"], truths=[mc_r_med, q_r_med, a1z_r_med, a2z_r_med])
  elif ringdown_pe == 'QNM':
    corner.corner(samples_r, labels=["A", "$\\tau$", "$f_{qnm}$", "$\phi_0$"], truths=[A_r_med, tau_r_med, f_qnm_r_med, phi0_r_med])
  plt.savefig('%s/img/PE_ringdown_corner_postburnin.png'%outdir, dpi=300)

  # Plotting best-fit waveform
  plt.figure(figsize=(16, 5))
  plt.subplot(131)
  plt.plot(t, h, label='waveform')
  plt.plot(t, d, alpha=0.2, label='data')
  plt.plot(t_i, h_i_rec, label='recovered inspiral')
  plt.plot(t_r, h_r_rec, label='recovered post-inspiral')
  plt.legend(loc='best')
  plt.subplot(132)
  plt.plot(t_i, h_i, label='waveform')
  plt.plot(t_i, h_i_rec, ls='--', label='recovered inspiral')
  plt.plot(t_i, d_i, alpha=0.2, label='inspiral data')
  plt.legend(loc='best')
  plt.subplot(133)
  plt.plot(t_r, h_r, label='waveform')
  plt.plot(t_r, h_r_rec, ls='--',label='recovered post-inspiral')
  plt.plot(t_r, d_r, alpha=0.2, label='ringdown data')
  plt.legend(loc='best')
  plt.savefig('%s/img/PE_bestfit_waveform_postburnin.png'%outdir, dpi=300)

  if ringdown_pe == 'lalsimulation':
    plt.figure()
    plt.loglog(estimated_psd.sample_frequencies, estimated_psd)
    plt.loglog(psd.sample_frequencies, psd)
    plt.xlim([flow, srate/2.])
    plt.ylim([1e-47, 1e-45])
    plt.savefig('%s/img/PE_psd_fits.png'%outdir, dpi=300)

  print '... plotting done'

  end_time = time.time()
  print '... time taken %.2f seconds'%(end_time - start_time)
