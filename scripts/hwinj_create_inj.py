#!/home/haris.k/devel/VirtualPycbc_Oct2017/bin/python

'''
Script to create injection
'''
import argparse
import os
import numpy as np
import logging

parser = argparse.ArgumentParser(description='Script to create injection txt')
parser.add_argument('-geocentric_end_time', '--geocentric_end_time', help='geocentric_end_time',type=int, required=True)
parser.add_argument('-gps_start_time', '--gps_start_time', help='gps_start_time',type=int, required=True)
parser.add_argument('-gps_end_time', '--gps_end_time', help='gps_end_time',type=int, required=True)

args = parser.parse_args()
geocentric_end_time =  args.geocentric_end_time
gps_start_time = args.gps_start_time
gps_end_time = args.gps_end_time



m1 = 34.0687565856
m2 = 28.1859215189
inclination = 0.576431783756
ra = 0.82572753427 # longitude
dec = -0.789415698817    #latitude
a1 = 0.463873715754
a2 = 0.43164035178
a1z = 0.0678055755045
a2z = 0.0161760018521
snr = 17.7961237513
distance = 539.31687166
#evet_gps_time = 1186741688.99
a1xysqr = a1**2 - a1z**2
a2xysqr = a2**2 - a2z**2

a1x = 0. #np.sqrt(a1xysqr/2.)
a2x = 0. #np.sqrt(a2xysqr/2.)
a1y= 0.
a2y= 0.

cline='pycbc_generate_hwinj  --instruments H1 L1 V1 --approximant IMRPhenomPv2 --order pseudoFourPN --waveform-low-frequency-cutoff 15 --mass1 34.0687565856 --mass2 28.1859215189 --spin1x 0. --spin1y 0. --spin1z 0.0678055755045 --spin2x 0. --spin2y 0. --spin2z  0.0161760018521 --inclination 0.576431783756 --polarization 1.09326457842 --ra 0.82572753427 --dec -0.789415698817 --sample-rate H1:16384 L1:16384 V1:16384 --frame-type H1:H1_HOFT_C00 L1:L1_HOFT_C00 V1:V1Online  --channel-name H1:GDS-CALIB_STRAIN L1:GDS-CALIB_STRAIN V1:Hrec_hoft_16384Hz --taper TAPER_START --network-snr 17.7961237513 --low-frequency-cutoff 20.0 --high-frequency-cutoff 1024.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --pad-data 8 --geocentric-end-time %d --gps-start-time %d --gps-end-time %d --strain-high-pass 0.1'%(geocentric_end_time,gps_start_time,gps_end_time)

os.system(cline)
