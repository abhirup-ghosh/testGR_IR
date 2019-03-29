import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os, lal, commands
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

#################################################################################################
# Injection parameters
#################################################################################################
m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj = 36., 29., 500., 0., 0., 0.
gpsTime, ra, dec, psi, detector =  1126285216.0, 0., 0., 0., 'H1'
approx = 'SEOBNRv4_opt' # common for injection and recovery
fit_formula = 'nospin_Pan2011'#'nonprecspin_Healy2014'
srate, flow = 2048, 20. # common for injection and recovery
#sigma = 1.e-23; asd_file = None # only required for whitenoise injections
asd_file = 'AdL'

#################################################################################################
# Running PE
#################################################################################################
ringdown_pe = 'lalsimulation'
#ringdown_pe = 'QNM'
num_threads = 32
n_steps_mcmc_imr = 50
n_steps_burnin_imr = 10
n_steps_mcmc_i = 50
n_steps_burnin_i = 10
n_steps_mcmc_r = 50
n_steps_burnin_r = 10

# lalsimulation prior
chirpmass_min, chirpmass_max = 1., 100.
q_min, q_max = 0.0001, 1.
dL_min, dL_max = 1., 1000.
iota_min, iota_max = -180., 180.
t0_min, t0_max = -10., 10.
phiref_min, phiref_max = -180., 180.

# ringdown priors
t0_ring = 0.
iA_min, A_max =  0., 10.#1.e-23, 1.e-19
tau_min, tau_max =  0.1e-3, 100.e-3 # seconds
f_qnm_min, f_qnm_max = 100., 400. # Hz ##FIXME##
phi0_min, phi0_max =  0., 2.*np.pi

#################################################################################################
# output directory
#################################################################################################
date = '20180128'
tag = 'normal'
outdir = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/6_param_runs_20Hz/%s_%s_m1_%d_m2_%d_dL_%.2f_iota_%.2f_t0_%.2f_phiref_%.2f_flow_%d_srate_%d_%s_t0_ring_0ms_%s/%s_test'%(date, approx, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, flow, srate, asd_file, ringdown_pe, tag)
#outdir = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/6_param_runs_20Hz/%s_%s_m1_%d_m2_%d_dL_%.2f_iota_%.2f_t0_%.2f_phiref_%.2f_flow_%d_srate_%d_white%s_t0_ring_0ms_%s/%s'%(date, approx, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, flow, srate, sigma, ringdown_pe, tag)

if ringdown_pe == 'lalsimulation':
  if asd_file is not None:
	run_cmd = 'python td_likelihood.py --asd-file %s --num-threads %d --m1-inj %f --m2-inj %f --dL-inj %f --iota-inj %f --t0-inj %f --phiref-inj %f --ra-inj %f --dec-inj %f --psi-inj %f --approximant %s --fit-formula %s --sample-rate %d --flow %f --detector %s --gpsTime %f --outdir %s/pe --ringdown_PE %s --chirpmass-min %f --chirpmass-max %f --q-min %f --q-max %f --dL-min %f --dL-max %f --iota-min %f --iota-max %f --t0-min %f --t0-max %f --phiref-min %f --phiref-max %f --imr-mcmc %d --imr-burnin %d --insp-mcmc %d --insp-burnin %d --ring-mcmc %d --ring-burnin %d --t0-ring %f'%(asd_file, num_threads, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, ra, dec, psi, approx, fit_formula, srate, flow, detector, gpsTime, outdir, ringdown_pe, chirpmass_min, chirpmass_max, q_min, q_max, dL_min, dL_max, iota_min, iota_max, t0_min, t0_max, phiref_min, phiref_max, n_steps_mcmc_imr, n_steps_burnin_imr, n_steps_mcmc_i, n_steps_burnin_i, n_steps_mcmc_r, n_steps_burnin_r, t0_ring)
  else: 
	run_cmd = 'python td_likelihood.py --white-noise-sigma %e --num-threads %d --m1-inj %f --m2-inj %f --dL-inj %f --iota-inj %f --t0-inj %f --phiref-inj %f --ra-inj %f --dec-inj %f --psi-inj %f --approximant %s --fit-formula %s --sample-rate %d --flow %f --detector %s --gpsTime %d --outdir %s/pe --ringdown_PE %s --chirpmass-min %f --chirpmass-max %f --q-min %f --q-max %f --dL-min %f --dL-max %f --iota-min %f --iota-max %f --t0-min %f --t0-max %f --phiref-min %f --phiref-max %f --imr-mcmc %d --imr-burnin %d --insp-mcmc %d --insp-burnin %d --ring-mcmc %d --ring-burnin %d --t0-ring %f'%(sigma, num_threads, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, ra, dec, psi, approx, fit_formula, srate, flow, detector, gpsTime, outdir, ringdown_pe, chirpmass_min, chirpmass_max, q_min, q_max, dL_min, dL_max, iota_min, iota_max, t0_min, t0_max, phiref_min, phiref_max, n_steps_mcmc_imr, n_steps_burnin_imr, n_steps_mcmc_i, n_steps_burnin_i, n_steps_mcmc_r, n_steps_burnin_r, t0_ring)

if ringdown_pe == 'QNM':
  if asd_file is not None:
        run_cmd = 'python td_likelihood.py --asd-file %s --num-threads %d --m1-inj %f --m2-inj %f --dL-inj %f --iota-inj %f --t0-inj %f --phiref-inj %f --ra-inj %f --dec-inj %f --psi-inj %f --approximant %s --fit-formula %s --sample-rate %d --flow %f --detector %s --gpsTime %f --outdir %s/pe --ringdown_PE %s --chirpmass-min %f --chirpmass-max %f --q-min %f --q-max %f --dL-min %f --dL-max %f --iota-min %f --iota-max %f --t0-min %f --t0-max %f --phiref-min %f --phiref-max %f --imr-mcmc %d --imr-burnin %d --insp-mcmc %d --insp-burnin %d --ring-mcmc %d --ring-burnin %d --t0-ring %f --A-min %e --A-max %e --tau-min %f --tau-max %f --f_qnm-min %f --f_qnm-max %f --phi0-min %f --phi0-max %f'%(asd_file, num_threads, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, ra, dec, psi, approx, fit_formula, srate, flow, detector, gpsTime, outdir, ringdown_pe, chirpmass_min, chirpmass_max, q_min, q_max, dL_min, dL_max, iota_min, iota_max, t0_min, t0_max, phiref_min, phiref_max, n_steps_mcmc_imr, n_steps_burnin_imr, n_steps_mcmc_i, n_steps_burnin_i, n_steps_mcmc_r, n_steps_burnin_r, t0_ring, A_min, A_max, tau_min, tau_max, f_qnm_min, f_qnm_max, phi0_min, phi0_max)
  else:
        run_cmd = 'python td_likelihood.py --white-noise-sigma %e --num-threads %d --m1-inj %f --m2-inj %f --dL-inj %f --iota-inj %f --t0-inj %f --phiref-inj %f --ra-inj %f --dec-inj %f --psi-inj %f --approximant %s --fit-formula %s --sample-rate %d --flow %f --detector %s --gpsTime %d --outdir %s/pe --ringdown_PE %s --chirpmass-min %f --chirpmass-max %f --q-min %f --q-max %f --dL-min %f --dL-max %f --iota-min %f --iota-max %f --t0-min %f --t0-max %f --phiref-min %f --phiref-max %f --imr-mcmc %d --imr-burnin %d --insp-mcmc %d --insp-burnin %d --ring-mcmc %d --ring-burnin %d --t0-ring %f --A-min %e --A-max %e --tau-min %f --tau-max %f --f_qnm-min %f --f_qnm-max %f --phi0-min %f --phi0-max %f'%(igma, num_threads, m1_inj, m2_inj, dL_inj, iota_inj, t0_inj, phiref_inj, ra, dec, psi, approx, fit_formula, srate, flow, detector, gpsTime, outdir, ringdown_pe, chirpmass_min, chirpmass_max, q_min, q_max, dL_min, dL_max, iota_min, iota_max, t0_min, t0_max, phiref_min, phiref_max, n_steps_mcmc_imr, n_steps_burnin_imr, n_steps_mcmc_i, n_steps_burnin_i, n_steps_mcmc_r, n_steps_burnin_r, t0_ring, A_min, A_max, tau_min, tau_max, f_qnm_min, f_qnm_max, phi0_min, phi0_max)


print(run_cmd)
os.system(run_cmd)
exit()
#################################################################################################
# IMR consistency test
#################################################################################################
post_samples_i = post_samples_imr = outdir+'/pe/data/inspiral_samples_postburnin.dat'
post_samples_r = outdir+'/pe/data/ringdown_samples_postburnin.dat'
insp_fhigh = ring_flow = 0.

run_cmd = 'python imrtgr_imr_consistency_test_final.py --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s/imrtgr --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401' %(post_samples_i, post_samples_r, post_samples_imr, fit_formula, outdir, m1_inj, m2_inj, 0., 0., insp_fhigh, ring_flow, approx)
print(run_cmd)
os.system(run_cmd)
