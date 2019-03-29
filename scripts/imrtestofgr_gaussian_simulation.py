""" 
Test the calculation in the IMR test of GR using Gaussian data 

- Generate Gaussian posteriors P_I(Mf, af), P_R(Mf, af) 
- Compute P_Delta(Delta Mf, Delta af)
- Compute P(Delta Mf/Mf, Delta af/af) 

P. Ajith, 2015-08-10
"""

import numpy.random as rand 
import numpy as np 
import matplotlib.pyplot as plt 
import plotsettings 
from scipy import interpolate as intp
from scipy import integrate as intg
import scipy.signal as ss

""" compute the integrant of P(dMf/Mf, daf/af). """
def P_integrant(af, Mf, v1, v2): 
		
	# Create dMf and daf vectors corresponding to the given v1 and v2. These vectors have to be 
	# monotonically increasing in order to evaluate the interpolated prob densities. Hence, for 
	# v1, v2 < 0, flip them, evaluate the prob density (in column or row) and flip it back 
	dMf = v1*Mf 
	daf = v2*af 

	if v1 < 0.:
		dMf = np.flipud(dMf) 
	if v2 < 0.:
		daf = np.flipud(daf) 

	P_delta = P_dMfdaf_intp(dMf, daf)

	if v1 < 0.:
		P_delta = np.fliplr(P_delta) 
	if v2 < 0.:
		P_delta = np.flipud(P_delta) 

	P_imr = P_Mfaf_imr_intp(Mf, af)
	return P_imr*P_delta*abs(Mf_mat*af_mat), P_imr, P_delta  

""" compute P(dMf/Mf, daf/af). """
def calc_sum(Mf, af, v1, v2,): 

	Pintg, P_imr, P_delta = P_integrant(af, Mf, v1, v2)
	return np.sum(Pintg)


######################################################################################
################################# MAIN PROGRAM #######################################
######################################################################################

# parameters of (Mf, af) distribution - inspiral analysis  
MfI0 = 50. 
afI0 = 0.75
sigma_Mf_i = 3 
sigma_af_i = 0.1 
c0I = 0.5

# parameters of (Mf, af) distribution - ringdown analysis  
MfR0 = 50. 
afR0 = 0.75
sigma_Mf_r = 3 
sigma_af_r = 0.15 
c0R = 0.75

# parameters of (Mf, af) distribution - IMR analysis  
Mf0 = 50. 
af0 = 0.75
sigma_Mf = 1 
sigma_af = 0.02 
c0 = 0.1

marksize = 1 
N_sampl = 10000 # number of posterior samples 
N_bins = 201 			# number of bins in histograms 
tag = '2017-03-19_norm_by_mean_test%d' %N_bins

####################################################################################
################################## generate samples ################################
####################################################################################

# generate samples of (Mf, af) - inspiral analysis 
C_i = [[sigma_Mf_i**2, c0I*sigma_Mf_i*sigma_af_i],[c0I*sigma_Mf_i*sigma_af_i, sigma_af_i**2]] 
Mf_i, af_i = np.random.multivariate_normal([MfI0, afI0], C_i, N_sampl).T 

# generate samples of (Mf, af) - ringdown analysis 
C_r = [[sigma_Mf_r**2, c0R*sigma_Mf_r*sigma_af_r],[c0R*sigma_Mf_r*sigma_af_r, sigma_af_r**2]] 
Mf_r, af_r = np.random.multivariate_normal([MfR0, afR0], C_r, N_sampl).T 

# generate samples of (Mf, af) - IMR analysis 
C = [[sigma_Mf**2, c0*sigma_Mf*sigma_af],[c0*sigma_Mf*sigma_af, sigma_af**2]] 
Mf_imr, af_imr = np.random.multivariate_normal([Mf0, af0], C, N_sampl).T 

print '... generated posterior samples' 

# impose a cut on af -- to make the distrbutions non-Gaussian, and to mimic the eta = 0.25 cutoff 
idx = np.logical_and(af_i < 1, af_r < 1) 
Mf_i = Mf_i[idx]
Mf_r = Mf_r[idx]
af_i = af_i[idx]
af_r = af_r[idx]
Mf_imr = Mf_imr[idx]
af_imr = af_imr[idx]

# parameter vectors at which histograms are computed 
Mf_limit = max(abs(np.append([Mf_i, Mf_r], Mf_imr)))
af_limit = max(abs(np.append([af_i, af_r], af_imr)))
Mf_bins = np.linspace(-Mf_limit, Mf_limit, N_bins)
af_bins = np.linspace(-af_limit, af_limit, N_bins)
dMfbyMf_bins = np.linspace(-1., 1., N_bins)
dafbyaf_bins = np.linspace(-1., 1., N_bins) 

# compute the prob density of (Mf, af) from various analyses 
P_Mfaf_i, Mf_edges, af_edges = np.histogram2d(Mf_i, af_i, bins=(Mf_bins, af_bins), normed=True)
P_Mfaf_r, Mf_edges, af_edges = np.histogram2d(Mf_r, af_r, bins=(Mf_bins, af_bins), normed=True)
P_Mfaf_imr, Mf_edges, af_edges = np.histogram2d(Mf_imr, af_imr, bins=(Mf_bins, af_bins), normed=True)
P_Mfaf_i = P_Mfaf_i.T
P_Mfaf_r = P_Mfaf_r.T
P_Mfaf_imr = P_Mfaf_imr.T

# generate sampels of Delta_Mf, Delta_af (where Delta Mf := Mf_i - Mf_r and af := af_i - af_r) and computes its prob density 
delta_Mf = Mf_i-Mf_r
delta_af = af_i-af_r 
P_dMfdaf, dMf_edges, daf_edges = np.histogram2d(delta_Mf, delta_af, bins=(Mf_bins, af_bins), normed=True)
P_dMfdaf = P_dMfdaf.T

# generate sampels of mean_Mf, mean_af (where mean_Mf := (Mf_i + Mf_r)/2 and mean_af := (af_i + af_r)/2 and computes its prob density 
mean_Mf = (Mf_i+Mf_r)/2.
mean_af = (af_i+af_r)/2.
P_meanMfdaf, meanMf_edges, meanaf_edges = np.histogram2d(mean_Mf, mean_af, bins=(Mf_bins, af_bins), normed=True)
P_meanMfdaf = P_meanMfdaf.T

# Generate the samples of delta_Mf/Mf and delta_af/af, and their histograms  
z1 = delta_Mf/Mf_imr
z2 = delta_af/af_imr
P_z1z2, z1_edges, z2_edges = np.histogram2d(z1, z2, bins=(dMfbyMf_bins, dafbyaf_bins), normed=True)
P_z1z2 = P_z1z2.T

print '... comptued histograms' 
####################################################################################
############### compute the prob densities semi analytically #######################
####################################################################################

# compute interpolation objects for the Mf,af posterior 
Mf_intp = (Mf_edges[:-1] + Mf_edges[1:])/2.
af_intp = (af_edges[:-1] + af_edges[1:])/2.
dMf_intp = (dMf_edges[:-1] + dMf_edges[1:])/2.
daf_intp = (daf_edges[:-1] + daf_edges[1:])/2.

P_Mfaf_i_intp = intp.interp2d(Mf_intp, af_intp, P_Mfaf_i, fill_value=0., bounds_error=False)
P_Mfaf_r_intp = intp.interp2d(Mf_intp, af_intp, P_Mfaf_r, fill_value=0., bounds_error=False)
P_Mfaf_imr_intp = intp.interp2d(Mf_intp, af_intp, P_Mfaf_imr, fill_value=0., bounds_error=False)
#P_dMfdaf_intp = intp.interp2d(dMf_intp, daf_intp, P_dMfdaf, fill_value=0., bounds_error=False)

print '... constructed interpolation objects' 

# now compute the P_Delta(Delta Mf, Delta af) from the inspiral and rigndown posteriors 
d_Mf = np.mean(np.diff(Mf_bins))
d_af = np.mean(np.diff(af_bins))

# compute using the 2d correlation function 
P_dMfdaf_calc = d_Mf*d_af*ss.correlate2d(P_Mfaf_i, P_Mfaf_r, boundary='fill', mode='same')

P_dMfdaf_calc = P_dMfdaf_calc.T
P_dMfdaf_intp = intp.interp2d(dMf_intp, daf_intp, P_dMfdaf_calc, fill_value=0., bounds_error=False)

# compute using the wd convolution (evaluate the second function g at g(-x)) 
#P_dMfdaf_calc = d_Mf*d_af*ss.convolve2d(P_Mfaf_i, np.flipud(np.fliplr(P_Mfaf_r)), boundary='fill', mode='same')

print '... computed P(delta_Mf, delta_af)' 

# compute the P(dMf/Mf, daf/af) by evaluating the integral 
dx = np.mean(np.diff(dMfbyMf_bins))
dy = np.mean(np.diff(dafbyaf_bins))
P_dMfbyMf_dafbyaf = np.zeros(shape=(N_bins,N_bins))
Mf_mat, af_mat = np.meshgrid(Mf_intp, af_intp)

for i, v2 in enumerate(dafbyaf_bins): 
	for j, v1 in enumerate(dMfbyMf_bins): 
		P_dMfbyMf_dafbyaf[i,j] = calc_sum(Mf_intp, af_intp, v1, v2)*dx*dy

print '... computed P(delta_Mf/Mf, delta_af/af)' 
####################################################################################
################################ plot the results  #################################
####################################################################################

plt.close()
plt.figure(figsize=(22,12))

# samples 
plt.subplot(361)
plt.plot(Mf_i, af_i, 'k.', mew=0, ms=marksize)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_i),max(Mf_i))
plt.ylim(min(af_i),max(af_i))
plt.xlim(30,80)
plt.ylim(0,1)
plt.title('Inspiral')

plt.subplot(362)
plt.plot(Mf_r, af_r, 'k.', mew=0, ms=marksize)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_r),max(Mf_r))
plt.ylim(min(af_r),max(af_r))
plt.xlim(30,80)
plt.ylim(0,1)
plt.title('Ringdown')

plt.subplot(363)
plt.plot(Mf_imr, af_imr, 'k.', mew=0, ms=marksize)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_imr),max(Mf_imr))
plt.ylim(min(af_imr),max(af_imr))
plt.ylim(0,1)
plt.xlim(30,80)
plt.title('IMR')

plt.subplot(364)
plt.plot(delta_Mf, delta_af, 'k.', mew=0, ms=marksize)
plt.xlabel('$\Delta M_f$')
plt.ylabel('$\Delta a_f$')
plt.xlim(min(delta_Mf), max(delta_Mf))
plt.ylim(min(delta_af), max(delta_af))
plt.title('$P_\Delta(\Delta M_f, \Delta a_f)$')

plt.subplot(365)
plt.plot(z1, z2, 'k.', mew=0, ms=marksize)
plt.xlim(min(z1), max(z1))
plt.ylim(min(z2), max(z2))
plt.ylim(-1,1)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')
plt.title('$P(\Delta M_f/M_f, \Delta a_f/a_f)$')

plt.subplot(366)
plt.plot(mean_Mf, mean_af, 'k.', mew=0, ms=marksize)
plt.ylim(-1,1)
plt.xlabel('$\\bar M_f$')
plt.ylabel('$\\ a_f$')
plt.title('$P(\\bar M_f, \\bar a_f)$')

plt.subplot(367)
plt.pcolormesh(Mf_edges, af_edges, P_Mfaf_i, cmap='YlOrBr')
plt.xlim(min(Mf_i),max(Mf_i))
plt.ylim(min(af_i),max(af_i))
plt.xlim(30,80)
plt.ylim(0,1)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
#plt.colorbar()

plt.subplot(368)
plt.pcolormesh(Mf_edges, af_edges, P_Mfaf_r, cmap='YlOrBr')
plt.xlim(min(Mf_r),max(Mf_r))
plt.ylim(min(af_r),max(af_r))
plt.xlim(30,80)
plt.ylim(0,1)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
#plt.colorbar()

plt.subplot(369)
plt.pcolormesh(Mf_edges, af_edges, P_Mfaf_imr, cmap='YlOrBr')
plt.xlim(min(Mf_imr),max(Mf_imr))
plt.ylim(min(af_imr),max(af_imr))
plt.xlim(30,80)
plt.ylim(0,1)
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
#plt.colorbar()

plt.subplot(3,6,10)
plt.pcolormesh(dMf_edges, daf_edges, P_dMfdaf, cmap='YlOrBr')
plt.xlabel('$\Delta M_f$')
plt.ylabel('$\Delta a_f$')
plt.xlim(min(delta_Mf), max(delta_Mf))
plt.ylim(min(delta_af), max(delta_af))
plt.clim(0,0.25)
#plt.colorbar()

plt.subplot(3,6,11)
plt.pcolormesh(z1_edges, z2_edges, P_z1z2, cmap='YlOrBr')
plt.xlim(min(z1), max(z1))
plt.ylim(min(z2), max(z2))
plt.ylim(-1,1)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')
#plt.colorbar()

plt.subplot(3,6,12)

plt.subplot(3,6,13)
plt.pcolormesh(Mf_intp, af_intp, P_Mfaf_i_intp(Mf_intp, af_intp), cmap='YlOrBr')
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_i),max(Mf_i))
plt.ylim(min(af_i),max(af_i))
plt.xlim(30,80)
plt.ylim(0,1)
#plt.colorbar()

plt.subplot(3,6,14)
plt.pcolormesh(Mf_intp, af_intp, P_Mfaf_r_intp(Mf_intp, af_intp), cmap='YlOrBr')
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_r),max(Mf_r))
plt.ylim(min(af_r),max(af_r))
plt.xlim(30,80)
plt.ylim(0,1)
#plt.colorbar()

plt.subplot(3,6,15)
plt.pcolormesh(Mf_intp, af_intp, P_Mfaf_imr_intp(Mf_intp, af_intp), cmap='YlOrBr')
plt.xlabel('$M_f$')
plt.ylabel('$a_f$')
plt.xlim(min(Mf_imr),max(Mf_imr))
plt.ylim(min(af_imr),max(af_imr))
plt.xlim(30,80)
plt.ylim(0,1)
#plt.colorbar()

plt.subplot(3,6,16)
plt.pcolormesh(dMf_edges, daf_edges, P_dMfdaf_calc, cmap='YlOrBr')
plt.xlabel('$\Delta M_f$')
plt.ylabel('$\Delta a_f$')
plt.xlim(min(delta_Mf), max(delta_Mf))
plt.ylim(min(delta_af), max(delta_af))
plt.clim(0,0.25)
#plt.colorbar()

plt.subplot(3,6,17)
plt.pcolormesh(dMfbyMf_bins, dafbyaf_bins, P_dMfbyMf_dafbyaf, cmap='YlOrBr')
#plt.colorbar()
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')
plt.xlim(min(z1), max(z1))
plt.ylim(min(z2), max(z2))
plt.ylim(-1,1)

plt.subplot(3,6,18)

print '... made plots' 

plt.savefig('2d_division_gaussian_%s_v3.png' %tag, dpi=300)
print '... saved plots' 

