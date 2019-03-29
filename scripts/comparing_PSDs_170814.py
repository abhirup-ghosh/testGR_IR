import matplotlib as mpl
mpl.use('Agg')
import os, numpy as np
import matplotlib.pyplot as plt
import lal
from pylal import Fr
import scipy
from scipy import signal as ss
from scipy import interpolate

laptop = False

if laptop:
  home = '/Users/abhirupghosh'
  data_loc = home+'/data'
else:
  home = '/home/abhirup'
  data_loc = '/home/data'

outdir = home+'/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/PSD_comparisons'
os.system('mkdir -p %s'%outdir)

H1_channel = 'H1:GDS-CALIB_STRAIN_T1700403_v1'
L1_channel = 'L1:GDS-CALIB_STRAIN_T1700403_v1'

H1_frame = data_loc+'/clean/GW170814/T1700403/H-H1_HOFT_C00_T1700403_v1-1186739813-4096.gwf'
L1_frame = data_loc+'/clean/GW170814/T1700403/L-L1_HOFT_C00_T1700403_v1-1186739813-4096.gwf'

H1_BW_psd_T1700403_v1_loc = home+'/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/PSDs/BayesWave_PSD_T1700403_v1_H1.dat'
H1_pycbc_psd_loc = home+'/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/PSDs/pycbc_H1-psd.txt'

H1_BW_freq, H1_BW_psd = np.loadtxt(H1_BW_psd_T1700403_v1_loc, unpack=True)
H1_pycbc_freq, H1_pycbc_psd = np.loadtxt(H1_pycbc_psd_loc, unpack=True)
H1_BW_psd_interp = scipy.interpolate.interp1d(H1_BW_freq, H1_BW_psd, fill_value=0., bounds_error=False)
H1_pycbc_psd_interp = scipy.interpolate.interp1d(H1_pycbc_freq, H1_pycbc_psd, fill_value=0., bounds_error=False)

ht_H1, gps_start_H1, xoffset_H1, dt_H1, xunit_H1, yunit_H1 = Fr.frgetvect(H1_frame, H1_channel)
t_H1 = np.arange(0., len(ht_H1))*dt_H1

plt.figure()
plt.plot(t_H1, ht_H1)
plt.savefig(outdir+'/H1_data.png')
plt.close()

ht_H1_4s_segment_list = np.split(ht_H1, indices_or_sections=1024)
t_H1_4s_segment_list = np.split(t_H1, indices_or_sections=1024)

plt.figure()
for idx in range(len(t_H1_4s_segment_list[:1024])):
  plt.plot(t_H1_4s_segment_list[idx], ht_H1_4s_segment_list[idx])
  plt.hold(True)
plt.plot(t_H1, ht_H1, 'k-', alpha=0.1)
plt.savefig(outdir+'/td_noise.png')
plt.close()

hf_H1_4s_segment_list = []
hf_H1_4s_segment_BW_whitened_list = []
hf_H1_4s_segment_pycbc_whitened_list = []

for (t_H1_4s_segment, ht_H1_4s_segment) in zip(t_H1_4s_segment_list,ht_H1_4s_segment_list):
  hf_H1_4s_segment = np.fft.fft(ht_H1_4s_segment)
  f_H1_4s_segment = np.fft.fftfreq(len(t_H1_4s_segment), d=dt_H1[0])
  hf_H1_4s_segment_list.append(hf_H1_4s_segment)
  
  hf_H1_4s_segment_BW_whitened = hf_H1_4s_segment/np.sqrt(H1_BW_psd_interp(f_H1_4s_segment))
  hf_H1_4s_segment_BW_whitened_list.append(hf_H1_4s_segment_BW_whitened)

  hf_H1_4s_segment_pycbc_whitened = hf_H1_4s_segment/np.sqrt(H1_pycbc_psd_interp(f_H1_4s_segment))
  hf_H1_4s_segment_pycbc_whitened_list.append(hf_H1_4s_segment_pycbc_whitened)

plt.figure()
plt.loglog(f_H1_4s_segment, abs(hf_H1_4s_segment_list[32])**2.)
plt.loglog(f_H1_4s_segment, abs(hf_H1_4s_segment_list[64])**2.)
plt.loglog(f_H1_4s_segment, abs(hf_H1_4s_segment_list[128])**2.)
plt.loglog(f_H1_4s_segment, abs(hf_H1_4s_segment_list[512])**2.)
plt.savefig(outdir+'/fd_noise.png')

hf_H1_4s_segment_list = np.asarray(hf_H1_4s_segment_list)
hf_H1_4s_segment_BW_whitened_list = np.asarray(hf_H1_4s_segment_BW_whitened_list)
hf_H1_4s_segment_pycbc_whitened_list = np.asarray(hf_H1_4s_segment_pycbc_whitened_list)

idx, = np.where(f_H1_4s_segment == 100)
nk_100 = hf_H1_4s_segment_list[:, idx]
nk_100_BW_whitened = hf_H1_4s_segment_BW_whitened_list[:, idx]
nk_100_pycbc_whitened = hf_H1_4s_segment_pycbc_whitened_list[:, idx]

nk_stand_normal = np.linspace(-4., 4., 10000)
p_nk_stand_normal = (1./np.sqrt(2.*np.pi))*np.exp(-nk_stand_normal**2./2.)

plt.figure(figsize=(16,5))
plt.subplot(131)
plt.hist(nk_100.real, bins=100)
plt.subplot(132)
plt.hist(nk_100_BW_whitened.real, bins=100, color='k', histtype='step', label='BW whitened')
plt.hist(nk_100_pycbc_whitened.real, bins=100, color='r', histtype='step', label='pycbc whitened')
plt.legend(loc='best')
plt.subplot(133)
plt.hist(nk_100_BW_whitened.real/65536/4, bins=100, color='k', histtype='step', label='BW whitened', normed=True)
plt.hist(nk_100_pycbc_whitened.real/65536/4, bins=100, color='r', histtype='step', label='pycbc whitened', normed=True)
plt.plot(nk_stand_normal, p_nk_stand_normal, 'k-')
plt.savefig(outdir+'/nk_100Hz.png')
