import lal
import json
import numpy as np
import sys
sys.path.append('../src')
import standard_gwtransf as gw
import os
import pycbc
import pycbc.noise
import pycbc.psd
from pycbc.frame import query_and_read_frame
from pycbc.types import timeseries
from pycbc.psd import inverse_spectrum_truncation, interpolate
from pycbc.waveform import get_fd_waveform, get_td_waveform
import pycbc.filter.matchedfilter as mfilter
from pycbc import frame as Fr
import scipy.signal
import imrtgrutils_final as tgr
import nr_fits as nr
import glob

data = np.genfromtxt('./SXS_campaign.dat', names=True, dtype=None)
ifo_map = {'HL':'H1 L1','HLV':'H1 L1 V1'}

for idx in range(len(data)):
    print data[idx]

    name, tag, det_config, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, iota, pol, ra, dec, gps_start_time, gps_end_time, geocentric_end_time, hw_inj_start_time, Mf, af, f_isco_Kerr, snr = data[idx]

    path = glob.glob('/home/abhirup/ligo-nr-data/lvcnr-lfs/SXS/%s*.h5'%name)[-1]

    pycbc_generate_hwinj_command = "pycbc_generate_hwinj --approximant IMRPhenomPv2 --order pseudoFourPN --waveform-low-frequency-cutoff 25 --mass1 %f --mass2 %f --spin1x %f --spin1y %f --spin1z %f --spin2x %f --spin2y %f --spin2z %f --inclination %f --polarization %f --ra %f --dec %f --instruments %s --sample-rate H1:16384 L1:16384 V1:16384 --taper TAPER_START --network-snr %f --low-frequency-cutoff 25.0 --high-frequency-cutoff 1024.0 --psd-model H1:aLIGOZeroDetHighPower L1:aLIGOZeroDetHighPower V1:AdVDesignSensitivityP1200087 --geocentric-end-time %d --gps-start-time %d --gps-end-time %d --strain-high-pass 1 --psd-output H1:%s_H1_psd.txt L1:%s_L1_psd.txt V1:%s_V1_psd.txt --tag %s"%(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, iota, pol, ra, dec, ifo_map[det_config], snr, geocentric_end_time, gps_start_time, gps_end_time, tag, tag, tag, tag)

    os.system(pycbc_generate_hwinj_command)
    print '... ran pycbc_generate_hwinj command: ', pycbc_generate_hwinj_command
    
    if not glob.glob('%s_*.xml.gz'%tag):
      hw_inj_start_time = 0
      print '... file does not exist'
    else:
      hw_inj_start_time = int(glob.glob('%s_*.xml.gz'%tag)[0][-17:-7])
      print '... file exists'

    pycbc_insert_frame_hwinj_h1 = 'pycbc_insert_frame_hwinj --fake-strain zeroNoise --low-frequency-cutoff 25 --gps-start-time %d --gps-end-time %d --pad-data 8 --sample-rate 16384 --hwinj-file %s_%d_H1.txt  --hwinj-start-time %d --ifo H1 --output-file %s_H-H1HWINJ.gwf --strain-high-pass 1'%(gps_start_time, gps_end_time, tag, hw_inj_start_time,  hw_inj_start_time, tag)

    pycbc_insert_frame_hwinj_l1 = 'pycbc_insert_frame_hwinj --fake-strain zeroNoise --low-frequency-cutoff 25 --gps-start-time %d --gps-end-time %d --pad-data 8 --sample-rate 16384 --hwinj-file %s_%d_L1.txt  --hwinj-start-time %d --ifo L1 --output-file %s_L-L1HWINJ.gwf --strain-high-pass 1'%(gps_start_time,gps_end_time, tag, hw_inj_start_time, hw_inj_start_time, tag)

    os.system(pycbc_insert_frame_hwinj_h1)
    print '... ran pycbc_insert_frame_hwinj for H1: ', pycbc_insert_frame_hwinj_h1

    os.system(pycbc_insert_frame_hwinj_l1)
    print '... ran pycbc_insert_frame_hwinj for L1: ', pycbc_insert_frame_hwinj_l1


    if det_config == 'HLV':
      pycbc_insert_frame_hwinj_v1 = 'pycbc_insert_frame_hwinj --fake-strain zeroNoise --low-frequency-cutoff 25 --gps-start-time %d --gps-end-time %d --pad-data 8 --sample-rate 16384 --hwinj-file %s_%d_V1.txt  --hwinj-start-time %d --ifo V1 --output-file %s_V-V1HWINJ.gwf --strain-high-pass 1'%(gps_start_time,gps_end_time, tag, hw_inj_start_time, hw_inj_start_time, tag)

      os.system(pycbc_insert_frame_hwinj_v1)
      print '... ran pycbc_insert_frame_hwinj for V1: ', pycbc_insert_frame_hwinj_v1
