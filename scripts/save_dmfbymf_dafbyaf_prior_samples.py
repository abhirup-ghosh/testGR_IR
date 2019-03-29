import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import bayesian as ba
import lalinference.imrtgr.imrtgrutils as tgr

N_bins = 200

Mf_i = np.random.uniform(1., 500., 10000000)
Mf_r = np.random.uniform(1., 500., 10000000)

af_i = np.random.uniform(-1., 1., 10000000)
af_r = np.random.uniform(-1., 1., 10000000)

dMf = Mf_i - Mf_r
daf = af_i - af_r

Mfbar = (Mf_i + Mf_r)/2.
afbar = (af_i + af_r)/2.

dMfbyMf_vec = np.linspace(-2., 2., N_bins)
dafbyaf_vec = np.linspace(-2., 2., N_bins)

dMfbyMf_intp = (dMfbyMf_vec[:-1] + dMfbyMf_vec[1:])/2.
dafbyaf_intp = (dafbyaf_vec[:-1] + dafbyaf_vec[1:])/2.

dMfbyMf = dMf / Mfbar
dafbyaf = daf / afbar

P_dMfbyMf_dafbyaf_pr, dMfbyMf_vec, dafbyaf_vec = np.histogram2d(dMfbyMf, dafbyaf, bins=(dMfbyMf_vec, dafbyaf_vec), normed=True)
P_dMfbyMf_dafbyaf_pr = P_dMfbyMf_dafbyaf_pr.T

s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf_pr, 0.68)
s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf_pr, 0.95)

plt.figure(figsize=(15,4))
plt.subplot(131)
plt.pcolormesh(dMfbyMf_intp, dafbyaf_intp, P_dMfbyMf_dafbyaf_pr, cmap='YlOrBr')
plt.contour(dMfbyMf_intp, dafbyaf_intp, P_dMfbyMf_dafbyaf_pr, levels=(s2_v1v2, s1_v1v2), cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim([-2, 2])
plt.ylim([-2, 2])
plt.grid()
plt.colorbar()
plt.title('Original')
plt.savefig('prior.png')
