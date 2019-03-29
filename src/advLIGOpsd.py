"""
To generate Adv-LIGO PSD and noise strain.

(C) Walter Del Pozzo. Modified by Archisman Ghosh, 2014-09-18.
"""

import numpy as np
import lalinspiral.sbank.psds as psds

def AdvLIGOPSD_wdp(f):
   if f<0.001: return 1
   x = f/215.;
   x2 = x*x;
   f10=f/10.0;
   f50=f/50.0;
   f100=f/100.0;
   f200=f/200.0;
   f300=f/300.0;
   f1000=f/1000.0;
   f2000=f/2000.0;
   x1=f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10*f10;
   x2=f50*f50*f50*f50*f50*f50;
   psd = f*((60000.0/x1)+5.0/x2+1.07*pow(f100,-3.25)+
            3.7*pow(f200,-1.25)+0.9*pow(f300,-0.08)+0.85*pow(f1000,0.8)+0.35*f2000*f2000*f2000);
   return  1.35e-50*psd

def AdvLIGOPSD(f):
  return psds.noise_models['AdvLIGOPsd'](f)

def mydiff(f):
  return np.append(np.append([f[1]-f[0]],(f[2:]-f[:-2])/2.),[f[-1]-f[-2]])

def AdvLIGOnoise_old(f):
  df = mydiff(f)
  noise = np.array([np.random.normal(0,0.5*np.sqrt(AdvLIGOPSD(f[ii]))*df[ii]) for ii in range(len(f))])
  return noise

def AdvLIGOnoise(f):
  df = mydiff(f)
  noise = np.array([np.random.normal(0,np.sqrt(AdvLIGOPSD(f[ii]))) for ii in range(len(f))])
  return noise

""" Usage """
if __name__=="__main__":
  import matplotlib.pyplot as plt
  
  #ff = np.exp(np.linspace(np.log(8),np.log(4096),4096))
  ff = np.linspace(8, 4096, 32768/4)
  ff_avg = np.mean(np.reshape(ff, (len(ff)/16, 16)), axis=1)
  
  psd = AdvLIGOPSD(ff)
  
  #noise_real_old = AdvLIGOnoise_old(ff)
  #noise_imag_old = AdvLIGOnoise_old(ff)
  #noise_power_old = 2*(noise_real_old**2 + noise_imag_old**2)/mydiff(ff)**2
  #noise_power_avg_old = np.mean(np.reshape(noise_power_old, (len(ff)/16, 16)), axis=1)
  
  noise_real = AdvLIGOnoise(ff)/np.sqrt(2)
  noise_imag = AdvLIGOnoise(ff)/np.sqrt(2)
  noise_complex = np.vectorize(complex)(noise_real, noise_imag)
  noise_power = np.abs(noise_complex)**2 # vnoise_real**2 + noise_imag**2
  noise_power_avg = np.mean(np.reshape(noise_power, (len(ff)/16, 16)), axis=1)
  
  noise_TD = np.real(np.sqrt(len(ff))*np.fft.ifft(noise_complex))
  
  plt.figure()
  plt.plot(ff_avg, noise_power_avg, 'g.')
  #plt.plot(ff_avg, noise_power_avg_old, 'c.')
  plt.loglog(ff, psd, 'r-')
  plt.xlabel('Frequency (Hz)')
  plt.ylabel('Noise PSD')
  
  #plt.figure()
  #plt.loglog(ff, np.abs(noise_real_old), 'c.')
  #plt.loglog(ff, np.sqrt(psd)*mydiff(ff), 'r-')
  #plt.xlabel('Frequency (Hz)')
  #plt.ylabel('Strain')
  
  plt.figure()
  plt.loglog(ff, np.abs(noise_real), 'g.')
  plt.loglog(ff, np.sqrt(psd), 'r-')
  plt.xlabel('Frequency (Hz)')
  plt.ylabel('Noise Strain')
  
  plt.figure()
  plt.plot(noise_TD, 'b.')
  plt.xlabel('Time (s)')
  plt.ylabel('Noise Strain')
  
  plt.show()
  
  
