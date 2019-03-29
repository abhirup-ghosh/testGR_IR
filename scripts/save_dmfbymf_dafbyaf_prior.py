""" 
Compute the prior distribution in the fractional mass and spin parameters (dMfbyMf, dchifbychif) 
corresponding to a uniform prior in final and spin of the final black hole. 

Abhirup Ghosh, 2016-06-24 

$Id:$
"""

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tgrplotsettings, time
import imrtestgr as tgr
import scipy
import scipy.signal as ss
from scipy import interpolate as interp
import pickle, gzip
from optparse import OptionParser

def P_Mfmean_afmean(Mf, af, Mf_min, Mf_max):

	dx = 2*Mf 
	dy = 2*af 
	xmax = Mf_max 
	xmin = Mf_min 	
	""" compute the distribution of the mean of two random variables """
	if dx > xmax + xmin && dx < 2 xmax && 1 < dy < 2:
		P = (-2 + dy)*(dx - 2 xmax)
	elif dx > 2 xmin && dx <= xmax + xmin && 1 < dy < 2:
		P = -(-2 + dy)*(dx - 2 xmin)
	elif dx > xmax + xmin && dx < 2 xmax && 0 < dy <= 1:
		P = -dy*(dx - 2 xmax)
	elif dx > 2 xmin && dx <= xmax + xmin && 0:
		P dy*(dx - 2 xmin)
	else:
		P = 0. 

t0 = time.time()

# read inputs from command line 
parser = OptionParser()
parser.add_option("--Mf-min", dest="Mf_min",
                  help="minimum value of the Mf prior [M_sun]")
parser.add_option("--Mf-max", dest="Mf_max",
                  help="maximum value of the Mf prior [M_sun]")
parser.add_option("--chif-min", dest="chif_min",
                  help="minimum value of the chif prior")
parser.add_option("--chif-max", dest="chif_max",
                  help="maximum value of the chif prior")
parser.add_option("--dMfbyMf-lim", dest="dMfbyMf_lim",
                  help="the prior distribution will be calculated from -dMfbyMf_lim to dMfbyMf_lim")
parser.add_option("--dchifbychif-lim", dest="dchifbychif_lim",
                  help="the prior distribution will be calculated from -dchifbychif_lim to -dchifbychif_lim")
parser.add_option("--num-threads", dest="num_threads",
                  help="number of threads to be used")
parser.add_option("--num-samples", dest="N_sampl",
                  help="number of prior samples to be used for computing the histogram")
parser.add_option("--num-bins", dest="N_bins",
                  help="number of bins to be used for computing the histogram")
(options, args) = parser.parse_args()

Mf_min = float(options.Mf_min)
Mf_max = float(options.Mf_max)
chif_min = float(options.chif_min)
chif_max = float(options.chif_max)
N_sampl = float(options.N_sampl)
N_bins = int(options.N_bins)
dMfbyMf_lim = float(options.dMfbyMf_lim)
dchifbychif_lim = float(options.dchifbychif_lim)
num_threads = int(options.num_threads)

print '... N_sampl = %e N_bins = %d' %(N_sampl, N_bins)

# compute the limits of integration for computing delta_Mf and delta_chif
Mf_lim = max(Mf_min, Mf_max)
chif_lim = max(chif_min, chif_max)

# the integral used to compute (Delta Mf, Delta af) has limits from -infinity to +infinity. We 
# are approximating this by setting the limits to (-Mf_lim to Mf_lim) and (-chif_lim to chif_lim)
# where Mf_lim and chif_lim are the max values of Mf and chif where the posteriors have nonzero 
# support. The scipy.signal.correlate2d function requires arguments x_bins and y_bins that need 
# to be symmetric around zero

Mf_bins = np.linspace(-Mf_max, Mf_max, N_bins)
chif_bins = np.linspace(-chif_max, chif_max, N_bins)

dMf = np.mean(np.diff(Mf_bins))
dchif = np.mean(np.diff(chif_bins))

Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
chif_intp = (chif_bins[:-1] + chif_bins[1:])/2.

# sample (Mf, chif) from a distribution uniform in Mf between (Mf_min, Mf_max) and chif between
# (chif_min, chif_max)

Mf_sampl_i = np.random.uniform(Mf_min, Mf_max, N_sampl)
chif_sampl_i = np.random.uniform(chif_min, chif_max, N_sampl)
Mf_sampl_r = np.random.uniform(Mf_min, Mf_max, N_sampl)
chif_sampl_r = np.random.uniform(chif_min, chif_max, N_sampl)

epsilon_sampl = 2*(Mf_sampl_i-Mf_sampl_r)/(Mf_sampl_i+Mf_sampl_r)
sigma_sampl = 2*(chif_sampl_i-chif_sampl_r)/(chif_sampl_i+chif_sampl_r)

# compute the 2D posterior distributions
P_Mfchif_pr, Mf_bins, chif_bins = np.histogram2d(Mf_sampl_i, chif_sampl_i, bins=(Mf_bins, chif_bins), normed=True)


# transpose to go from (X,Y) indexing returned by np.histogram2d() to array (i,j) indexing for further
# computations. From now onwards, different rows (i) correspond to different values of Mf and different 
# columns (j) correspond to different values of chif 
P_Mfchif_pr = P_Mfchif_pr.T

# compute the posterior of (delta_Mf/Mf, delta_chif/chif)
P_dMfdchif_pr = dMf*dchif*ss.correlate2d(P_Mfchif_pr, P_Mfchif_pr, boundary='fill', mode='same')
P_dMfdchif_pr = P_dMfdchif_pr.T
print '... computed P(delta_Mf, delta_chif)'

# compute interpolation objects for the Mf,chif posterior and delta_Mf and delta_chif posterior
P_dMfdchif_pr_interp_object = scipy.interpolate.interp2d(Mf_intp, chif_intp, P_dMfdchif_pr, fill_value=0., bounds_error=False)
P_Mfchif_imr_pr_interp_object = scipy.interpolate.interp2d(Mf_intp, chif_intp, P_Mfchif_pr, fill_value=0., bounds_error=False)

P_dMfdchif_pr_interp = P_dMfdchif_pr_interp_object(Mf_intp, chif_intp)

# defining limits of delta_Mf/Mf and delta_chif/chif
dMfbyMf_vec = np.linspace(-dMfbyMf_lim, dMfbyMf_lim, N_bins)
dchifbychif_vec = np.linspace(-dchifbychif_lim, dchifbychif_lim, N_bins)

# compute the P(dMf/Mf, dchif/chif) by evaluating the integral
diff_dMfbyMf = np.mean(np.diff(dMfbyMf_vec))
diff_dchifbychif = np.mean(np.diff(dchifbychif_vec))
P_dMfbyMf_dchifbychif_pr = np.zeros(shape=(N_bins,N_bins))

# compute the distribution of epsilon and sigma by directly histogramming hte samples 
P_epssigma_pr, eps_bins, sigma_bins =  np.histogram2d(epsilon_sampl, sigma_sampl, bins=(dMfbyMf_vec, dchifbychif_vec), normed=True)
P_epssigma_pr = P_epssigma_pr.T

eps_intp = (eps_bins[:-1] + eps_bins[1:])/2.
sig_intp = (sigma_bins[:-1] + sigma_bins[1:])/2.

P_epssigma_pr_interp_obj = interp.interp2d(eps_intp, sig_intp, P_epssigma_pr, fill_value=0., bounds_error=False)
P_epssigma_pr_interp = P_epssigma_pr_interp_obj(eps_intp, sig_intp)

print '... started computing (delta_Mf/Mf, delta_chif/chif) prior'
# compute the posterior on the fractional deviation parameters (delta_Mf/Mf, delta_chif/chif). 
# Approximate the integral in Eq.(6) of the document LIGO-P1500185-v5 by a discrete sum
for i, v2 in enumerate(dchifbychif_vec):
    for j, v1 in enumerate(dMfbyMf_vec):
      P_dMfbyMf_dchifbychif_pr[i,j] = tgr.calc_sum(Mf_intp, chif_intp, v1, v2, P_dMfdchif_pr_interp_object, P_Mfchif_imr_pr_interp_object)*diff_dMfbyMf*diff_dchifbychif

P_dMfbyMf_dchifbychif_pr /= np.sum(P_dMfbyMf_dchifbychif_pr) * diff_dMfbyMf * diff_dchifbychif
print '... computed (delta_Mf/Mf, delta_chif/chif) prior'

# create an interpolation object and save it 
outfile = '../data/Prior_dMfbyMf_dchifbychif_Mfmin_%.1f_Mfmax_%.1f_chifmin%.1f_chifmax%.1f_dMfbyMflim%.1f_dchifbychiflim%.1f'%(Mf_min, Mf_max, chif_min, chif_max, dMfbyMf_lim, dchifbychif_lim)
P_dMfbyMf_dchifbychif_pr_interp_obj = interp.interp2d(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr, fill_value=0., bounds_error=False)
f = gzip.open(outfile+".pklz",'wb')
pickle.dump(P_dMfbyMf_dchifbychif_pr_interp_obj, f)

print '... saved the interpolation object.'

# read the interpolation object, reconstruct the data from the interpolation object 
f = gzip.open(outfile+".pklz",'rb')
P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)
P_dMfbyMf_dchifbychif_pr_interp = P_dMfbyMf_dchifbychif_pr_interp_obj(dMfbyMf_vec, dchifbychif_vec)

# difference between the original and interpolated data 
interp_err = abs(P_dMfbyMf_dchifbychif_pr - P_dMfbyMf_dchifbychif_pr_interp)

print '... maximum difference between the original and interpolated data is %e' %np.max(interp_err)

dMfbyMf_lim = 0.9
dchifbychif_lim = 0.9

plt.figure(figsize=(15,4))

plt.subplot(131)
plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr, cmap='YlOrBr')
plt.contour(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('Original')

plt.subplot(132)
plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_interp, cmap='YlOrBr')
plt.contour(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_interp, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('Interpolation')

plt.subplot(133)
plt.pcolormesh(eps_intp, sig_intp, P_epssigma_pr, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('direct histogram')
plt.savefig(outfile+'.png', dpi=300)
plt.close()


# compute the power of N 
N = 50 
P_dMfbyMf_dchifbychif_pr_powN = P_dMfbyMf_dchifbychif_pr**N
P_dMfbyMf_dchifbychif_pr_powN /= np.sum(P_dMfbyMf_dchifbychif_pr_powN) * diff_dMfbyMf * diff_dchifbychif

P_dMfbyMf_dchifbychif_pr_interp_powN = P_dMfbyMf_dchifbychif_pr_interp**N
P_dMfbyMf_dchifbychif_pr_interp_powN /= np.sum(P_dMfbyMf_dchifbychif_pr_interp_powN) * diff_dMfbyMf * diff_dchifbychif

P_dMfbyMf_dchifbychif_pr_dir_powN = P_epssigma_pr_interp**N
P_dMfbyMf_dchifbychif_pr_dir_powN /= np.sum(P_dMfbyMf_dchifbychif_pr_dir_powN) * diff_dMfbyMf * diff_dchifbychif

# save power of N 
plt.figure(figsize=(15,4))
plt.subplot(131)
plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_powN, cmap='YlOrBr')
plt.contour(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_powN, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('Original')

plt.subplot(132)
plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_interp_powN, cmap='YlOrBr')
plt.contour(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_interp_powN, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('Interpolation')

plt.subplot(133)
plt.pcolormesh(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif_pr_dir_powN, cmap='YlOrBr')
plt.contour(eps_intp, sig_intp, P_dMfbyMf_dchifbychif_pr_dir_powN, cmap='YlOrBr')
plt.xlabel('$\Delta M_f / M_f$')
plt.ylabel('$\Delta \chi _f / \chi _f$')
plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
plt.ylim(-dchifbychif_lim, dchifbychif_lim)
plt.grid()
plt.colorbar()
plt.title('Difference')
plt.savefig(outfile+'_powN.png', dpi=300)
plt.close()


print '... saved fig (time taken: %f secs)' %(time.time()-t0)
plt.show()
