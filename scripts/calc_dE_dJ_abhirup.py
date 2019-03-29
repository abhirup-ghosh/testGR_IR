import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import loadline as ll
import scipy
import scipy.ndimage.filters as filter
import scipy.signal as ss
from optparse import OptionParser
import plotsettings
import bayesian as ba
import time
import os
import calc_dE_dJ_debug_plots as dp

def calc_final_mass_spin(mtot, eta, fit_formula):
        """"
        Calculate the mass and spin of the final black hole using a fitting formula obtained 
        from NR simulations 
        """
        if fit_formula == 'nospin_pan2011':
                mf = mtot*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.))
                af = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)
        else:
                print '# unknown spin fit formula'
                exit()
        return mf, af

def gf(P):
	return filter.gaussian_filter(P, sigma=2.0)


### MAIN PROGRAM #### 
start_time = time.time()

# read inputs from command line 
parser = OptionParser()
parser.add_option("-i", "--insp-post", dest="insp_post",
                  help="file containing the posterior samples from the lalinference inspiral run")
parser.add_option("-r", "--ring-post", dest="ring_post",
                  help="file containing the posterior samples from the lalinference ringdown run")
parser.add_option("-m", "--imr-post", dest="imr_post",
                  help="file containing the posterior samples from the full lalinference IMR run")
parser.add_option("-f", "--fit-formula", dest="fit_formula",
                  help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_pan2011'")
parser.add_option("-o", "--out-dir", dest="out_dir",
                  help="output directory")
parser.add_option("-d", "--debug-plots", dest="debug_plots", help="")
parser.add_option("-M", "--mtot_inj", dest="M", help="total injected mass")
parser.add_option("-Q", "--mratio_inj", dest="Q", help="injected mass ratio")
parser.add_option("--comp_mass_prior_min", dest="comp_mass_prior_min", help="component mass prior min")
parser.add_option("--comp_mass_prior_max", dest="comp_mass_prior_max", help="component mass prior max")




(options, args) = parser.parse_args()
insp_post = options.insp_post
ring_post = options.ring_post
imr_post = options.imr_post
out_dir = options.out_dir
fit_formula = options.fit_formula
debug_plots = options.debug_plots
M = options.M
Q = options.Q
comp_mass_prior_min = options.comp_mass_prior_min 
comp_mass_prior_max = options.comp_mass_prior_max

debug_plots = False


# create output directory and copy the script 
os.system('mkdir -p %s' %out_dir)
os.system('cp %s %s' %(__file__, out_dir))


M = float(M)
Q = float(Q)
ETA = Q/(1.+Q)**2.
M_F, A_F = calc_final_mass_spin(M, ETA, fit_formula)

###############################################################################################
# Read the posteriors from the inspiral, ringdown and imr lalinference runs (after post-processing) 
###############################################################################################

legend = ll.load_line(insp_post, dtype='string')
m1_i, m2_i, logl_i = np.loadtxt(insp_post, skiprows =1, usecols = (list(legend).index('m1'), list(legend).index('m2'), list(legend).index('logl')), unpack=True)
mtot_i = m1_i+m2_i
eta_i = m1_i*m2_i/mtot_i**2.
mf_i, af_i = calc_final_mass_spin(mtot_i, eta_i, fit_formula)

legend = ll.load_line(ring_post, dtype='string')
m1_r, m2_r, logl_r = np.loadtxt(ring_post, skiprows =1, usecols = (list(legend).index('m1'), list(legend).index('m2'), list(legend).index('logl')), unpack=True)
mtot_r = m1_r+m2_r
eta_r = m1_r*m2_r/mtot_r**2.
mf_r, af_r = calc_final_mass_spin(mtot_r, eta_r, fit_formula)

legend = ll.load_line(imr_post, dtype='string')
m1_imr, m2_imr, logl_imr = np.loadtxt(imr_post, skiprows =1, usecols = (list(legend).index('m1'), list(legend).index('m2'), list(legend).index('logl')), unpack=True)
mtot_imr = m1_imr+m2_imr
eta_imr = m1_imr*m2_imr/mtot_imr**2.
mf_imr, af_imr = calc_final_mass_spin(mtot_imr, eta_imr, fit_formula)

print min(mf_i), max(mf_i), min(af_i), max(af_i)
print min(mf_r), max(mf_r), min(af_r), max(af_r)
print min(mf_imr), max(mf_imr), min(af_imr), max(af_imr)

print 'read posteriors...'

print 'correcting posteriors...'
###############################################################################################
# METHOD 1
###############################################################################################
mf_min_or = min(min(mf_i), min(mf_r), min(mf_imr))
mf_max_or = max(max(mf_i), max(mf_r), max(mf_imr))
mf_lim_or = max(abs(mf_min_or), abs(mf_max_or))

af_min_or = min(min(af_i), min(af_r), min(af_imr))
af_max_or = max(max(af_i), max(af_r), max(af_imr))
af_lim_or = max(abs(af_min_or), abs(af_max_or))

mf_bins_or = 500
af_bins_or = 500


mf_an_or = np.linspace(-mf_lim_or, mf_lim_or, mf_bins_or)
af_an_or = np.linspace(-af_lim_or, af_lim_or, af_bins_or)

dmf_or = np.mean(np.diff(mf_an_or))
daf_or = np.mean(np.diff(af_an_or))

# generate random samples of the prior uniform in component masses 
m1_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)
m2_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)

# compute the corrresponding samples of prior in mf, af 
mtot_pr =  m1_pr + m2_pr
eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)
mf_pr, af_pr = calc_final_mass_spin(mtot_pr, eta_pr, fit_formula)

# compute the 2D prior distribution in Mf and af 
P_mfaf_pr, mf_an_or, af_an_or = np.histogram2d(mf_pr, af_pr, bins=(mf_an_or, af_an_or), normed=True)

# compute the 2D posterior distributions for the inspiral, ringodwn and IMR analyses 
P_mfaf_i, mf_an_or, af_an_or = np.histogram2d(mf_i, af_i, bins=(mf_an_or, af_an_or), normed=True)
P_mfaf_r, mf_an_or, af_an_or = np.histogram2d(mf_r, af_r, bins=(mf_an_or, af_an_or), normed=True)
P_mfaf_imr, mf_an_or, af_an_or = np.histogram2d(mf_imr, af_imr, bins=(mf_an_or, af_an_or), normed=True)

# compute the corrected 2D posteriors in Mf and af by dividing by the prior distribution 
P_mfaf_i = P_mfaf_i/P_mfaf_pr
P_mfaf_r = P_mfaf_r/P_mfaf_pr
P_mfaf_imr = P_mfaf_imr/P_mfaf_pr

P_mfaf_i[np.isnan(P_mfaf_i)] = 0.
P_mfaf_r[np.isnan(P_mfaf_r)] = 0.
P_mfaf_imr[np.isnan(P_mfaf_imr)] = 0.

P_dmfdaf_or = dmf_or*daf_or*ss.correlate2d(P_mfaf_i, P_mfaf_r, boundary='fill', mode='same')
#########################################################################################
# plotting                                                                              #
#########################################################################################

# details for contours
s1_i = ba.nsigma_value(gf(P_mfaf_i))
s2_i = ba.nsigma_value(gf(P_mfaf_i),2)

s1_r = ba.nsigma_value(gf(P_mfaf_r))
s2_r = ba.nsigma_value(gf(P_mfaf_r),2)

s1_imr = ba.nsigma_value(gf(P_mfaf_imr))
s2_imr = ba.nsigma_value(gf(P_mfaf_imr),2)

plt.figure(figsize=(16,16))
plt.subplot(431)
plt.pcolormesh(mf_an_or, af_an_or, gf(P_mfaf_i).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_or[:-1], af_an_or[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.xlim([min(mf_i), max(mf_i)])
plt.ylim([min(af_i), max(af_i)])
plt.title('Inspiral')
plt.grid()

plt.subplot(432)

plt.pcolormesh(mf_an_or, af_an_or, gf(P_mfaf_r).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_or[:-1], af_an_or[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_r), max(mf_r)])
plt.ylim([min(af_r), max(af_r)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('Ringdown')
plt.grid()

plt.subplot(433)
plt.pcolormesh(mf_an_or, af_an_or, gf(P_mfaf_imr).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_or[:-1], af_an_or[:-1], gf(P_mfaf_imr).T, levels=(s1_imr,s2_imr), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_imr), max(mf_imr)])
plt.ylim([min(af_imr), max(af_imr)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('IMR')
plt.grid()

print 'Method 1 done...'


###############################################################################################
# METHOD 2
###############################################################################################
mf_bins_i = 50
af_bins_i = 50

mf_bins_r = 100
af_bins_r = 100

mf_bins_imr = 20
af_bins_imr = 20

mf_an_i = np.linspace(min(mf_i), max(mf_i), mf_bins_i)
af_an_i = np.linspace(min(af_i), max(af_i), af_bins_i)

mf_an_r = np.linspace(min(mf_r), max(mf_r), mf_bins_r)
af_an_r = np.linspace(min(af_r), max(af_r), af_bins_r)

mf_an_imr = np.linspace(min(mf_imr), max(mf_imr), mf_bins_imr)
af_an_imr = np.linspace(min(af_imr), max(af_imr), af_bins_imr)

# generate random samples of the prior uniform in component masses 
m1_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)
m2_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)

# compute the corrresponding samples of prior in mf, af 
mtot_pr =  m1_pr + m2_pr
eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)
mf_pr, af_pr = calc_final_mass_spin(mtot_pr, eta_pr, fit_formula)

P_mfaf_pr_i, mf_an_i, af_an_i = np.histogram2d(mf_pr, af_pr, bins=(mf_an_i, af_an_i), normed=True)
P_mfaf_i, mf_an_i, af_an_i = np.histogram2d(mf_i, af_i, bins=(mf_an_i, af_an_i), normed=True)
P_mfaf_i = P_mfaf_i/P_mfaf_pr_i

P_mfaf_pr_r, mf_an_r, af_an_r = np.histogram2d(mf_pr, af_pr, bins=(mf_an_r, af_an_r), normed=True)
P_mfaf_r, mf_an_r, af_an_r = np.histogram2d(mf_r, af_r, bins=(mf_an_r, af_an_r), normed=True)
P_mfaf_r = P_mfaf_r/P_mfaf_pr_r

P_mfaf_pr_imr, mf_an_imr, af_an_imr = np.histogram2d(mf_pr, af_pr, bins=(mf_an_imr, af_an_imr), normed=True)
P_mfaf_imr, mf_an_imr, af_an_imr = np.histogram2d(mf_imr, af_imr, bins=(mf_an_imr, af_an_imr), normed=True)
P_mfaf_imr = P_mfaf_imr/P_mfaf_pr_imr

P_mfaf_i[np.isnan(P_mfaf_i)] = 0.
P_mfaf_r[np.isnan(P_mfaf_r)] = 0.
P_mfaf_imr[np.isnan(P_mfaf_imr)] = 0.

#########################################################################################
# plotting                                                                              #
#########################################################################################

# details for contours
s1_i = ba.nsigma_value(gf(P_mfaf_i))
s2_i = ba.nsigma_value(gf(P_mfaf_i),2)

s1_r = ba.nsigma_value(gf(P_mfaf_r))
s2_r = ba.nsigma_value(gf(P_mfaf_r),2)

s1_imr = ba.nsigma_value(gf(P_mfaf_imr))
s2_imr = ba.nsigma_value(gf(P_mfaf_imr),2)

plt.subplot(434)
plt.pcolormesh(mf_an_i, af_an_i, gf(P_mfaf_i).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_i[:-1], af_an_i[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.xlim([min(mf_i), max(mf_i)])
plt.ylim([min(af_i), max(af_i)])
plt.title('Inspiral')
plt.grid()

plt.subplot(435)

plt.pcolormesh(mf_an_r, af_an_r, gf(P_mfaf_r).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_r[:-1], af_an_r[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_r), max(mf_r)])
plt.ylim([min(af_r), max(af_r)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('Ringdown')
plt.grid()

plt.subplot(436)
plt.pcolormesh(mf_an_imr, af_an_imr, gf(P_mfaf_imr).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_imr[:-1], af_an_imr[:-1], gf(P_mfaf_imr).T, levels=(s1_imr,s2_imr), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_imr), max(mf_imr)])
plt.ylim([min(af_imr), max(af_imr)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('IMR')
plt.grid()

print 'Method 2 done...'

###############################################################################################

P_mfaf_i_interp_object_RBS = scipy.interpolate.RectBivariateSpline(mf_an_i[:-1], af_an_i[:-1], P_mfaf_i)
P_mfaf_r_interp_object_RBS = scipy.interpolate.RectBivariateSpline(mf_an_r[:-1], af_an_r[:-1], P_mfaf_r)
P_mfaf_imr_interp_object_RBS = scipy.interpolate.RectBivariateSpline(mf_an_imr[:-1], af_an_imr[:-1], P_mfaf_imr)

P_mfaf_i_interp_object_I2D = scipy.interpolate.interp2d(mf_an_i[:-1], af_an_i[:-1], P_mfaf_i, fill_value=0., bounds_error=False, kind='cubic')
P_mfaf_r_interp_object_I2D = scipy.interpolate.interp2d(mf_an_r[:-1], af_an_r[:-1], P_mfaf_r, fill_value=0., bounds_error=False, kind='cubic')
P_mfaf_imr_interp_object_I2D = scipy.interpolate.interp2d(mf_an_imr[:-1], af_an_imr[:-1], P_mfaf_imr, fill_value=0., bounds_error=False, kind='cubic')


################################################################################################
# METHOD 3
################################################################################################

mf_min = min(min(mf_i), min(mf_r), min(mf_imr))
mf_max = max(max(mf_i), max(mf_r), max(mf_imr))
mf_lim = max(abs(mf_min), abs(mf_max))

af_min = min(min(af_i), min(af_r), min(af_imr))
af_max = max(max(af_i), max(af_r), max(af_imr))
af_lim = max(abs(af_min), abs(af_max))

mf_bins = 500
af_bins = 500


mf_an_RBS = np.linspace(-mf_lim, mf_lim, mf_bins)
af_an_RBS = np.linspace(-af_lim, af_lim, af_bins)

dmf = np.mean(np.diff(mf_an_RBS))
daf = np.mean(np.diff(af_an_RBS))

P_mfaf_i = P_mfaf_i_interp_object_RBS(mf_an_RBS[:-1], af_an_RBS[:-1])
P_mfaf_r = P_mfaf_r_interp_object_RBS(mf_an_RBS[:-1], af_an_RBS[:-1])
P_mfaf_imr = P_mfaf_imr_interp_object_RBS(mf_an_RBS[:-1], af_an_RBS[:-1])

#########################################################################################
# plotting                                                                              #
#########################################################################################

# details for contours
s1_i = ba.nsigma_value(gf(P_mfaf_i))
s2_i = ba.nsigma_value(gf(P_mfaf_i),2)

s1_r = ba.nsigma_value(gf(P_mfaf_r))
s2_r = ba.nsigma_value(gf(P_mfaf_r),2)

s1_imr = ba.nsigma_value(gf(P_mfaf_imr))
s2_imr = ba.nsigma_value(gf(P_mfaf_imr),2)

plt.subplot(437)
plt.pcolormesh(mf_an_RBS, af_an_RBS, gf(P_mfaf_i).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_RBS[:-1], af_an_RBS[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.xlim([min(mf_i), max(mf_i)])
plt.ylim([min(af_i), max(af_i)])
plt.title('Inspiral')
plt.grid()

plt.subplot(438)

plt.pcolormesh(mf_an_RBS, af_an_RBS, gf(P_mfaf_r).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_RBS[:-1], af_an_RBS[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_r), max(mf_r)])
plt.ylim([min(af_r), max(af_r)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('Ringdown')
plt.grid()

plt.subplot(439)
plt.pcolormesh(mf_an_RBS, af_an_RBS, gf(P_mfaf_imr).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_RBS[:-1], af_an_RBS[:-1], gf(P_mfaf_imr).T, levels=(s1_imr,s2_imr), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_imr), max(mf_imr)])
plt.ylim([min(af_imr), max(af_imr)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('IMR')
plt.grid()

P_dmfdaf_RBS = dmf*daf*ss.correlate2d(P_mfaf_i, P_mfaf_r, boundary='fill', mode='same')

print 'Method 3 done...'

################################################################################################
# METHOD 4
################################################################################################

mf_min = min(min(mf_i), min(mf_r), min(mf_imr))
mf_max = max(max(mf_i), max(mf_r), max(mf_imr))
mf_lim = max(abs(mf_min), abs(mf_max))

af_min = min(min(af_i), min(af_r), min(af_imr))
af_max = max(max(af_i), max(af_r), max(af_imr))
af_lim = max(abs(af_min), abs(af_max))

mf_bins = 500
af_bins = 500


mf_an_I2D = np.linspace(-mf_lim, mf_lim, mf_bins)
af_an_I2D = np.linspace(-af_lim, af_lim, af_bins)

dmf = np.mean(np.diff(mf_an_I2D))
daf = np.mean(np.diff(af_an_I2D))



P_mfaf_i = P_mfaf_i_interp_object_I2D(mf_an_I2D[:-1], af_an_I2D[:-1])
P_mfaf_r = P_mfaf_r_interp_object_I2D(mf_an_I2D[:-1], af_an_I2D[:-1])
P_mfaf_imr = P_mfaf_imr_interp_object_I2D(mf_an_I2D[:-1], af_an_I2D[:-1])

#########################################################################################
# plotting                                                                              #
#########################################################################################

# details for contours
s1_i = ba.nsigma_value(gf(P_mfaf_i))
s2_i = ba.nsigma_value(gf(P_mfaf_i),2)

s1_r = ba.nsigma_value(gf(P_mfaf_r))
s2_r = ba.nsigma_value(gf(P_mfaf_r),2)

s1_imr = ba.nsigma_value(gf(P_mfaf_imr))
s2_imr = ba.nsigma_value(gf(P_mfaf_imr),2)

plt.subplot(4,3,10)
plt.pcolormesh(mf_an_I2D, af_an_I2D, gf(P_mfaf_i).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_I2D[:-1], af_an_I2D[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.xlim([min(mf_i), max(mf_i)])
plt.ylim([min(af_i), max(af_i)])
plt.title('Inspiral')
plt.grid()

plt.subplot(4,3,11)

plt.pcolormesh(mf_an_I2D, af_an_I2D, gf(P_mfaf_r).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_I2D[:-1], af_an_I2D[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_r), max(mf_r)])
plt.ylim([min(af_r), max(af_r)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('Ringdown')
plt.grid()

plt.subplot(4,3,12)
plt.pcolormesh(mf_an_I2D, af_an_I2D, gf(P_mfaf_imr).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_I2D[:-1], af_an_I2D[:-1], gf(P_mfaf_imr).T, levels=(s1_imr,s2_imr), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([min(mf_imr), max(mf_imr)])
plt.ylim([min(af_imr), max(af_imr)])
plt.axvline(x=M_F, ls='--', lw=1, color='k')
plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('IMR')
plt.grid()
plt.tight_layout()
plt.savefig('%s/data.png'%(out_dir))

P_dmfdaf_I2D = dmf*daf*ss.correlate2d(P_mfaf_i, P_mfaf_r, boundary='fill', mode='same')


print 'Method 4 done...'

################################################################################################
# compute the posterior of (delta_mf, delta_af)
################################################################################################
s1_dmfdaf_or = ba.nsigma_value(gf(P_dmfdaf_or))
s2_dmfdaf_or = ba.nsigma_value(gf(P_dmfdaf_or),2)

s1_dmfdaf_RBS = ba.nsigma_value(gf(P_dmfdaf_RBS))
s2_dmfdaf_RBS = ba.nsigma_value(gf(P_dmfdaf_RBS),2)

s1_dmfdaf_I2D = ba.nsigma_value(gf(P_dmfdaf_I2D))
s2_dmfdaf_I2D = ba.nsigma_value(gf(P_dmfdaf_I2D),2)

plt.figure(figsize=(16,6))

plt.subplot(131)
plt.pcolormesh(mf_an_or, af_an_or, gf(P_dmfdaf_or).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_or[:-1], af_an_or[:-1], gf(P_dmfdaf_or).T, levels=(s1_dmfdaf_or,s2_dmfdaf_or), linewidth=0.1, cmap='hot')
plt.axvline(x=0, ls='--', lw=1, color='k')
plt.axhline(y=0, ls='--', lw=1, color='k')
plt.xlabel('$\Delta M_f~[M_\odot]$')
plt.ylabel('$\Delta a_f$')
plt.xlim([-M/2.,M/2.])
plt.ylim([-1.,1.])
plt.grid()

plt.subplot(132)
plt.pcolormesh(mf_an_RBS, af_an_RBS, gf(P_dmfdaf_RBS).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_RBS[:-1], af_an_RBS[:-1], gf(P_dmfdaf_RBS).T, levels=(s1_dmfdaf_RBS,s2_dmfdaf_RBS), linewidth=0.1, cmap='hot')
plt.axvline(x=0, ls='--', lw=1, color='k')
plt.axhline(y=0, ls='--', lw=1, color='k')
plt.xlabel('$\Delta M_f~[M_\odot]$')
plt.ylabel('$\Delta a_f$')
plt.xlim([-M/2.,M/2.])
plt.ylim([-1.,1.])
plt.grid()

plt.subplot(133)
plt.pcolormesh(mf_an_I2D, af_an_I2D, gf(P_dmfdaf_I2D).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an_I2D[:-1], af_an_I2D[:-1], gf(P_dmfdaf_I2D).T, levels=(s1_dmfdaf_I2D,s2_dmfdaf_I2D), linewidth=0.1, cmap='hot')
plt.axvline(x=0, ls='--', lw=1, color='k')
plt.axhline(y=0, ls='--', lw=1, color='k')
plt.xlabel('$\Delta M_f~[M_\odot]$')
plt.ylabel('$\Delta a_f$')
plt.xlim([-M/2.,M/2.])
plt.ylim([-1.,1.])
plt.grid()

plt.savefig('%s/delta.png'%(out_dir))

print 'calculated delta_mf, delta_af posteriors...'
