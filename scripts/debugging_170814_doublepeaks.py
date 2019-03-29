import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import corner

data_quality_list = ['cleaneddata', 'uncleaneddata']

i = 1
"""
plt.figure(figsize=(10,5))
for data_quality in data_quality_list:

  post_loc = '/home/abhirup/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/20170904_HLV_BWPSD_%s_finalruns/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1186741861.53-0/V1H1L1'%data_quality

  H1_loc = post_loc + '/H1'
  L1_loc = post_loc + '/L1'
  V1_loc = post_loc + '/V1'

  data_hlv = np.genfromtxt(post_loc+'/posterior_samples.dat', names=True, dtype=None)
  data_h = np.genfromtxt(H1_loc+'/posterior_samples.dat', names=True, dtype=None)
  data_l = np.genfromtxt(L1_loc+'/posterior_samples.dat', names=True, dtype=None)
  data_v = np.genfromtxt(V1_loc+'/posterior_samples.dat', names=True, dtype=None)

  plt.subplot(1,2,i)
  plt.scatter(data_h['m1'], data_h['m2'], lw=0, color='r', s=1, alpha=0.2, label='H1')
  plt.scatter(data_l['m1'], data_l['m2'], lw=0, color='b', s=1, alpha=0.2, label='L1')
  #plt.scatter(data_v['m1'], data_v['m2'], lw=0, color='orange', s=1,  alpha=0.2, label='V1')
  plt.legend(loc='best')
  plt.title('%s'%data_quality)
  plt.xlabel('m1')
  plt.ylabel('m2')

  i += 1

plt.savefig('debugging_170814_doublepeaks_hlv_scatterplots.png')
"""
post_loc = '/home/abhirup/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/20170904_HLV_BWPSD_cleaneddata_finalruns/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1186741861.53-0/V1H1L1'
data_hlv = np.genfromtxt(post_loc+'/posterior_samples.dat', names=True, dtype=None)

print np.shape(data_hlv)

corner.corner(data_hlv)
plt.savefig('debugging_170814_doublepeaks_hlv_cornerplot.png')
