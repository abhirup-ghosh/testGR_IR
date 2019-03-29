import numpy as np
import lal
import lalsimulation as lalsim
import h5py
from lal import MSUN_SI, PC_SI, PI

filepath = '/home/abhirup/Documents/Work/testGR_IR/NR_frames/lvcnr-lfs/SXS/SXS_BBH_0126_Res5.h5'

f = h5py.File(filepath, 'r')

params = lal.CreateDict()
lalsim.SimInspiralWaveformParamsInsertNumRelData(params, filepath)

# Metadata parameters masses:
# Note there is an intrinsic check the make sure the mass ratio is correct

m1 = f.attrs['mass1']
m2 = f.attrs['mass2']

# Choose extrinsic parameters:

mtotal = 56.00
m1SI = m1 * mtotal / (m1 + m2) * lal.MSUN_SI
m2SI = m2 * mtotal / (m1 + m2) * lal.MSUN_SI

deltaT = 1.0/4096.

f_lower = f.attrs['f_lower_at_1MSUN']/mtotal  # this generates the whole NR waveforms

fRef = f_lower   #beginning of the waveform
fStart = fRef

inclination = PI/3
distance = 100. * PC_SI * 1.0e6
phiRef = 0.0

f.close()

# Spins
# The NR spins need to be transformed into the lal frame. We wrote a function that does that for you:
spins = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(fRef, mtotal,filepath)

s1x = spins[0]
s1y = spins[1]
s1z = spins[2]
s2x = spins[3]
s2y = spins[4]
s2z = spins[5]

approx = lalsim.NR_hdf5

# Generating the waveform via ChooseTDWaveform():

hp, hc = lalsim.SimInspiralChooseTDWaveform(m1SI, m2SI, s1x, s1y, s1z,
           s2x, s2y, s2z, distance, inclination, phiRef, np.pi/2., 0.0, 0.0, 
           deltaT, fStart, fRef, params, approx)

print m1, m2
print s1x, s2x
print s1y, s2y
print s1z, s2z

print 'pycbc_generate_hwinj --numrel-data %s --instruments H1 L1 --approximant NR_hdf5 --order pseudoFourPN --waveform-low-frequency-cutoff 25 --mass1 %f --mass2 %f --spin1x %f --spin2x %f --spin1y %f --spin2y %f --spin1z %f --spin2z %f --inclination 2.83701491 --polarization 1.42892206 --ra -1.26157296 --dec 1.94972503 --sample-rate H1:16384 L1:16384 --frame-type H1:H1_HOFT_C00 L1:L1_HOFT_C00 --channel-name H1:GDS-CALIB_STRAIN L1:GDS-CALIB_STRAIN --taper TAPER_START --network-snr 25 --psd-low-frequency-cutoff 25.0 --psd-high-frequency-cutoff 1000.0 --psd-estimation median --psd-segment-length 16 --psd-segment-stride 8 --pad-data 8 --geocentric-end-time 1126259462.0 --gps-start-time 1126258960 --gps-end-time 1126259472 --strain-high-pass 1'%(filepath, m1*mtotal, m2*mtotal, s1x, s2x, s1y, s2y, s1z, s2z)

