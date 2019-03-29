import matplotlib as mpl
mpl.use('Agg')
import os, sys, commands, numpy as np
import extract_from_lalinferencenest_hdf5 as ex
import imrtestgr as tgr
import standard_gwtransf as gw
import matplotlib.pyplot as plt


def get_posterior_samples_dot_dat_file(hdf5_file_loc, outfile, approx='aligned'):
  possamp = ex.get_samples(hdf5_file_loc)
  mc, q, a1, a2 = possamp['mc'], possamp['q'], possamp['a1'], possamp['a2']
  m1, m2 = gw.comp_from_mcq(mc, q)
  if approx == 'aligned':
    a1z, a2z = a1, a2
    Mf, af = tgr.calc_final_mass_spin(m1, m2, a1z, a2z, 'nonprecspin_Healy2014')
  np.savetxt(outfile,np.c_[m1, m2, mc, q, a1z, a2z, Mf, af], header='m1 m2 mc q a1z a2z Mf af')


if __name__ == '__main__':
  base_dir = '../runs/td_likelihood/20170925_semi-final_runs/20170925_GW150914/SEOBNRv4_opt_m1_36_m2_29_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_ZERO_DET_high_P_t0_ring_0ms_FD/post-inspiral_H1_132Hz_fixparam/lalinferencenest/SEOBNRv4_ROMpseudoFourPN'
  hdf5_file_loc = base_dir + '/posterior_samples/posterior_H1_1126285216.000000000-0.hdf5'
  outfile = base_dir + '/1126285216.000000000-0/H1/posterior_samples.dat'
  approx = 'aligned'
  get_posterior_samples_dot_dat_file(hdf5_file_loc, outfile, approx='aligned')
