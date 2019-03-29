import numpy as np
import os, sys, commands
import matplotlib.pyplot as plt
import chirptimes as ct
import scipy.ndimage.filters as filter

snr_V1_r = np.array([])
snr_H1_r = np.array([])
snr_L1_r = np.array([])
snr_network_r = np.array([])
detectable_inj_r = np.array([])

approximant = 'SEOBNRv2_ROM_DoubleSpinpseudoFourPN'
in_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-15/uniform_compmass_spins_comoving_volume'

snr_thres = 12

inj_list = np.linspace(1,5000,5000)
inj_data = np.genfromtxt('../injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt', dtype=None, names=True)

for inj_no in inj_list:
  try:
    print '... reading injection %d'%inj_no
    in_dir = '%s/injection_%d'%(in_dir_root, inj_no)

    snr_IMR_loc = commands.getoutput('ls %s/IMR/engine/*_snr.txt'%(in_dir)).split()[0]
    snr_IMR_data = np.genfromtxt(snr_IMR_loc)
    snr_IMR_V1 = snr_IMR_data[0,1]; snr_IMR_L1 = snr_IMR_data[1,1]; snr_IMR_H1 = snr_IMR_data[2,1]; snr_IMR_network = snr_IMR_data[3,1]

    if (snr_IMR_H1 > snr_thres and snr_IMR_L1 > snr_thres):
	print '... passed threshold'
        detectable_inj_r = np.append(detectable_inj_r, inj_no)

  except:
    print '... data not found'

print len(detectable_inj_r)
np.savetxt('../data/paper2_useful_info/detected_events_H1L1_IMRsnrthreshold_12.txt', np.c_[detectable_inj_r], fmt=['%d'])
