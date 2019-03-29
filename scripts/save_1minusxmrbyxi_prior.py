""" 
Compute the prior distribution in the fractional mass and spin parameters (dMfbyMf, dafbyaf) 
corresponding to a uniform prior in final and spin of the final black hole. 

Abhirup Ghosh, 2016-06-24 

$Id:$
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tgrplotsettings, time
import imrtestgr as tgr
import scipy
import scipy.signal as ss
from scipy import interpolate as interp
import pickle, gzip
from optparse import OptionParser


t0 = time.time()

# read inputs from command line 
parser = OptionParser()
parser.add_option("--Mf-min", dest="Mf_min",
                  help="minimum value of the Mf prior [M_sun]")
parser.add_option("--Mf-max", dest="Mf_max",
                  help="maximum value of the Mf prior [M_sun]")
parser.add_option("--af-min", dest="af_min",
                  help="minimum value of the af prior")
parser.add_option("--af-max", dest="af_max",
                  help="maximum value of the af prior")
parser.add_option("--dMfbyMf-min", dest="dMfbyMf_min",
                  help="the prior distribution will be calculated from -dMfbyMf_min to 1")
parser.add_option("--dafbyaf-min", dest="dafbyaf_min",
                  help="the prior distribution will be calculated from -dafbyaf_min to 1")
parser.add_option("--num-threads", dest="num_threads",
                  help="number of threads to be used")
parser.add_option("--num-samples", dest="N_sampl",
                  help="number of prior samples to be used for computing the histogram")
parser.add_option("--num-bins", dest="N_bins",
                  help="number of bins to be used for computing the histogram")
(options, args) = parser.parse_args()

Mf_min = float(options.Mf_min)
Mf_max = float(options.Mf_max)
af_min = float(options.af_min)
af_max = float(options.af_max)
N_sampl = float(options.N_sampl)
N_bins = int(options.N_bins)
dMfbyMf_min = float(options.dMfbyMf_min)
dafbyaf_min = float(options.dafbyaf_min)
num_threads = int(options.num_threads)

print '... N_sampl = %e N_bins = %d' %(N_sampl, N_bins)

# compute the limits of integration for computing delta_Mf and delta_af
Mf_lim = max(Mf_min, Mf_max)
af_lim = max(af_min, af_max)

# the integral used to compute (Delta Mf, Delta af) has limits from -infinity to +infinity. We 
# are approximating this by setting the limits to (-Mf_lim to Mf_lim) and (-af_lim to af_lim)
# where Mf_lim and af_lim are the max values of Mf and af where the posteriors have nonzero 
# support. The scipy.signal.correlate2d function requires arguments x_bins and y_bins that need 
# to be symmetric around zero

Mf_bins = np.linspace(-Mf_max, Mf_max, N_bins)
af_bins = np.linspace(-af_max, af_max, N_bins)

dMf = np.mean(np.diff(Mf_bins))
daf = np.mean(np.diff(af_bins))

Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
af_intp = (af_bins[:-1] + af_bins[1:])/2.

# sample (Mf, af) from a distribution uniform in Mf between (Mf_min, Mf_max) and af between
# (af_min, af_max)

Mf_sampl = np.random.uniform(Mf_min, Mf_max, N_sampl)
af_sampl = np.random.uniform(af_min, af_max, N_sampl)

# compute the 2D posterior distributions
P_Mfaf_pr, Mf_bins, af_bins = np.histogram2d(Mf_sampl, af_sampl, bins=(Mf_bins, af_bins), normed=True)

# transpose to go from (X,Y) indexing returned by np.histogram2d() to array (i,j) indexing for further
# computations. From now onwards, different rows (i) correspond to different values of Mf and different 
# columns (j) correspond to different values of af 
P_Mfaf_pr = P_Mfaf_pr.T

P_Mfaf_i_pr_interp_object = scipy.interpolate.interp2d(Mf_intp, af_intp, P_Mfaf_pr, fill_value=0., bounds_error=False)
P_Mfaf_r_pr_interp_object = scipy.interpolate.interp2d(Mf_intp, af_intp, P_Mfaf_pr, fill_value=0., bounds_error=False)

# defining limits of delta_Mf/Mf and delta_af/af
MfRbyMfI_vec = np.linspace(-dMfbyMf_min+1., 0., N_bins)
afRbyafI_vec = np.linspace(-dafbyaf_min+1., 0., N_bins)
P_MfRbyMfI_afRbyafI_pr  = np.ones(shape=(N_bins,N_bins))

for i, v2 in enumerate(afRbyafI_vec):
    for j, v1 in enumerate(MfRbyMfI_vec):
        P_MfRbyMfI_afRbyafI_pr[i,j] = tgr.calc_sum(Mf_intp, af_intp, v1, v2, P_Mfaf_r_pr_interp_object, P_Mfaf_i_pr_interp_object) * dMf * daf


P_MfRbyMfI_afRbyafI_pr /= np.sum(P_MfRbyMfI_afRbyafI_pr) * dMf * daf
print '... computed (delta_Mf/Mf, delta_af/af) prior'

# create an interpolation object and save it 
outfile = '../data/Prior_MfRbyMfI_afRbyafI_Mfmin_%.1f_Mfmax_%.1f_afmin%.1f_afmax%.1f_dMfbyMfmin%.1f_dafbyafmin%.1f'%(Mf_min, Mf_max, af_min, af_max, dMfbyMf_min, dafbyaf_min)
P_MfRbyMfI_afRbyafI_pr_interp_obj = interp.interp2d(MfRbyMfI_vec, afRbyafI_vec, P_MfRbyMfI_afRbyafI_pr, fill_value=0., bounds_error=False)
f = gzip.open(outfile+".pklz",'wb')
pickle.dump(P_MfRbyMfI_afRbyafI_pr_interp_obj, f)

print '... saved the interpolation object.'

# read the interpolation object, reconstruct the data from the interpolation object 
f = gzip.open(outfile+".pklz",'rb')
P_MfRbyMfI_afRbyafI_pr_interp_obj = pickle.load(f)

MfRbyMfI_vec = np.flipud(MfRbyMfI_vec)
afRbyafI_vec = np.flipud(afRbyafI_vec)
P_MfRbyMfI_afRbyafI_pr_interp = P_MfRbyMfI_afRbyafI_pr_interp_obj(MfRbyMfI_vec, afRbyafI_vec)
P_MfRbyMfI_afRbyafI_pr_interp = np.fliplr(P_MfRbyMfI_afRbyafI_pr_interp)
P_MfRbyMfI_afRbyafI_pr_interp = np.flipud(P_MfRbyMfI_afRbyafI_pr_interp)
MfRbyMfI_vec = np.flipud(MfRbyMfI_vec)
afRbyafI_vec = np.flipud(afRbyafI_vec)

# difference between the original and interpolated data 
interp_err = abs(P_MfRbyMfI_afRbyafI_pr - P_MfRbyMfI_afRbyafI_pr_interp)

print '... maximum difference between the original and interpolated data is %e' %np.max(interp_err)

plt.figure(figsize=(15,4))
plt.subplot(131)
plt.pcolormesh(1-MfRbyMfI_vec, 1-afRbyafI_vec, P_MfRbyMfI_afRbyafI_pr, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta a_f / a_f$')
plt.xlim(dMfbyMf_min, 1.)
plt.ylim(dafbyaf_min, 1.)
plt.grid()
plt.colorbar()
plt.title('Original')

plt.subplot(132)
plt.pcolormesh(1-MfRbyMfI_vec, 1-afRbyafI_vec, P_MfRbyMfI_afRbyafI_pr_interp, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta a_f / a_f$')
plt.xlim(dMfbyMf_min, 1.)
plt.ylim(dafbyaf_min, 1.)
plt.grid()
plt.colorbar()
plt.title('Interpolation')

plt.subplot(133)
plt.pcolormesh(1-MfRbyMfI_vec, 1-afRbyafI_vec, interp_err, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta a_f / a_f$')
plt.xlim(dMfbyMf_min, 1.)
plt.ylim(dafbyaf_min, 1.)
plt.grid()
plt.clim(np.min(interp_err)-1e-16,np.max(interp_err)+1e-16)
plt.colorbar()
plt.title('Difference')

plt.savefig(outfile+'.png', dpi=300)
print '... saved fig (time taken: %f secs)' %(time.time()-t0)

