import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

data_mr = np.genfromtxt('/home/abhirup/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/20170904_HLV_BWPSD_cleaneddata_finalruns/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1186741861.53-0/V1H1L1/posterior_samples.dat', names=True, dtype=None)

m1, m2, mc, mtot = data_mr['m1'], data_mr['m2'], data_mr['mc'], data_mr['mtotal']

idxhigh, = np.where(mc > 16.)
idxlow, = np.where(mc <= 16.)

np.savetxt('/home/abhirup/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/20170904_HLV_BWPSD_cleaneddata_finalruns/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1186741861.53-0/V1H1L1/posterior_samples_mc_above_16.dat',data_mr[idxhigh])
np.savetxt('/home/abhirup/Documents/Work/O2/2017/August/14/1186741861p5268/G297595/lalinference/20170904_HLV_BWPSD_cleaneddata_finalruns/post-inspiral/lalinferencenest/IMRPhenomPv2pseudoFourPN/1186741861.53-0/V1H1L1/posterior_samples_mc_below_16.dat',data_mr[idxlow])

plt.figure(figsize=(8,8))
plt.subplot(221)
plt.hist(m1[idxhigh],histtype='step', alpha=0.2, color='r', bins=np.linspace(10,100,20),ls='dashed')
plt.hist(m1[idxlow],histtype='step', alpha=0.2, color='g', bins=np.linspace(10,100,20),ls='dashed')
plt.xlabel('m1')
plt.subplot(222)
plt.hist(m2[idxhigh],histtype='step', alpha=0.2, color='r',bins=np.linspace(0,50,20),ls='dashed')
plt.hist(m2[idxlow],histtype='step', alpha=0.2, color='g',bins=np.linspace(0,50,20),ls='dashed')
plt.xlabel('m2')
plt.subplot(223)
plt.hist(mc[idxhigh],histtype='step', alpha=0.2, color='r',bins=np.linspace(0,50,20),ls='dashed')
plt.hist(mc[idxlow],histtype='step', alpha=0.2, color='g',bins=np.linspace(0,50,20),ls='dashed')
plt.xlabel('mc')
plt.subplot(224)
plt.hist(mtot[idxhigh],histtype='step', alpha=0.2, color='r',ls='dashed', bins=np.linspace(10,100,20))
plt.hist(mtot[idxlow],histtype='step', alpha=0.2, color='g', ls='dashed',bins=np.linspace(10,100,20))
plt.xlabel('mtotal')
plt.tight_layout()
plt.savefig('170814_separated_post_mc_16.png')
