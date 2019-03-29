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
"""
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

debug_plots = True


# create output directory and copy the script 
os.system('mkdir -p %s' %out_dir)
os.system('cp %s %s' %(__file__, out_dir))

M = float(M)
Q = float(Q)
ETA = Q/(1.+Q)**2.
M_F, A_F = calc_final_mass_spin(M, ETA, fit_formula)
"""
###############################################################################################
# Read the posteriors from the inspiral, ringdown and imr lalinference runs (after post-processing) 
###############################################################################################
"""
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

"""
M = 100.

mean_i = [100.,0.5]; cov_i = [[100.,-0.35],[-0.35,0.0025]]
mean_r = [100.,0.5]; cov_r = [[225.,0.75],[0.75,0.01]]
mean_imr = [100.,0.5]; cov_imr = [[25.,0.],[0.,0.0004]]

mf_i, af_i = np.random.multivariate_normal(mean_i,cov_i,100000).T
mf_r, af_r = np.random.multivariate_normal(mean_r,cov_r,100000).T
mf_imr, af_imr = np.random.multivariate_normal(mean_imr,cov_imr,100000).T


print '... read posteriors'
###############################################################################################

###############################################################################################
# compute the limits of integration for computing delta_Mf and delta_af 
###############################################################################################
mf_min = min(min(mf_i), min(mf_r), min(mf_imr))
mf_max = max(max(mf_i), max(mf_r), max(mf_imr))
mf_lim = max(abs(mf_min), abs(mf_max))

af_min = min(min(af_i), min(af_r), min(af_imr))
af_max = max(max(af_i), max(af_r), max(af_imr))
af_lim = max(abs(af_min), abs(af_max))

mf_bins = 200
af_bins = 200


mf_an = np.linspace(-mf_lim, mf_lim, mf_bins)
af_an = np.linspace(-af_lim, af_lim, af_bins)

dmf = np.mean(np.diff(mf_an))
daf = np.mean(np.diff(af_an))
###############################################################################################

"""
###############################################################################################
# Undo the effect of the prior from the lalinference posterior. Lalinference assumes a        #
# uniform prior in component masses. We need to assume a uniform prior in Mf, af              #
###############################################################################################

# generate random samples of the prior uniform in component masses 
m1_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)
m2_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, 1e7)

# compute the corrresponding samples of prior in mf, af 
mtot_pr =  m1_pr + m2_pr
eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)
mf_pr, af_pr = calc_final_mass_spin(mtot_pr, eta_pr, fit_formula)

# compute the 2D prior distribution in Mf and af 
P_mfaf_pr, mf_an, af_an = np.histogram2d(mf_pr, af_pr, bins=(mf_an, af_an), normed=True)
for i in range(len(P_mfaf_pr)):
  for j in range(len(P_mfaf_pr)):
    if P_mfaf_pr[i,j]==0.:
	P_mfaf_pr[i,j]=1e-16
"""
# compute the 2D posterior distributions for the inspiral, ringodwn and IMR analyses 
P_mfaf_i, mf_an, af_an = np.histogram2d(mf_i, af_i, bins=(mf_an, af_an), normed=True)
P_mfaf_r, mf_an, af_an = np.histogram2d(mf_r, af_r, bins=(mf_an, af_an), normed=True)
P_mfaf_imr, mf_an, af_an = np.histogram2d(mf_imr, af_imr, bins=(mf_an, af_an), normed=True)

"""
if debug_plots == True: 
  plt.figure(figsize=(24,12))
  plt.subplot(241)
  plt.plot(m1_pr, m2_pr, 'k.', ms=0.5, mew=0.01) 
  plt.xlabel('$m_1~[M_\odot]$')
  plt.ylabel('$m_2~[M_\odot]$')
  plt.grid()
  plt.title('Prior samples')
  plt.subplot(242)
  plt.plot(mf_pr, af_pr, 'k.', ms=0.5, mew=0.01) 
  plt.xlabel('$M_f~[M_\odot]$')
  plt.ylabel('$a_f/M_f$')
  plt.xlim([M-0.2*M,M+0.2*M])
  plt.grid()
  plt.title('Prior samples')
  plt.subplot(243)
  plt.pcolormesh(mf_an, af_an, P_mfaf_pr.T, cmap='YlOrBr')
  plt.colorbar()
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.axvline(x=M_F, ls='--', lw=1, color='k')
  plt.axhline(y=A_F, ls='--', lw=1, color='k')
  plt.xlim([M-0.2*M,M+0.2*M])
  plt.ylim([0, max(af_an)])
  plt.clim(np.min(P_mfaf_pr), np.max(P_mfaf_pr))
  plt.title('Prior dist')
  plt.grid()

  plt.subplot(245)
  plt.plot(m1_r, m2_r, 'k.', ms=1, mew=0.01)
  plt.xlabel('$m_1 [M_{\odot}]$')
  plt.xlabel('$m_2 [M_{\odot}]$')
  plt.title('Ringdown posterior samples')
  plt.grid()

  plt.subplot(246)
  plt.plot(mf_r, af_r, 'k.', ms=1, mew=0.01)
  plt.axvline(x=M_F, ls='--', lw=1, color='k')
  plt.axhline(y=A_F, ls='--', lw=1, color='k')
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.grid()
  plt.title('Ringdown posterior samples')

  plt.subplot(247)
  plt.pcolormesh(mf_an[:-1], af_an[:-1], P_mfaf_r.T, cmap='YlOrBr')
  plt.colorbar()
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.axvline(x=M_F, ls='--', lw=1, color='k')
  plt.axhline(y=A_F, ls='--', lw=1, color='k')
  plt.xlim([M-0.2*M,M+0.2*M])
  plt.ylim([0, max(af_an)])
  plt.clim(np.min(P_mfaf_pr), np.max(P_mfaf_pr))
  plt.title('Ring down posterior (before correction)')
  plt.grid()

# compute the corrected 2D posteriors in Mf and af by dividing by the prior distribution 
P_mfaf_i = P_mfaf_i/P_mfaf_pr
P_mfaf_r = P_mfaf_r/P_mfaf_pr
P_mfaf_imr = P_mfaf_imr/P_mfaf_pr

if debug_plots == True:
  plt.subplot(248)
  plt.pcolormesh(mf_an[:-1], af_an[:-1], P_mfaf_r.T, cmap='YlOrBr')
  plt.colorbar()
  plt.xlabel('$M_f [M_{\odot}]$')
  plt.ylabel('$a_f/M_f$')
  plt.axvline(x=M_F, ls='--', lw=1, color='k')
  plt.axhline(y=A_F, ls='--', lw=1, color='k')
  plt.xlim([M-0.2*M,M+0.2*M])
  plt.ylim([0, max(af_an)])
  plt.clim(np.min(P_mfaf_pr), np.max(P_mfaf_pr))
  plt.savefig('%s/DebugPlots_Prior_M%2.1f_q%2.1f.png' %(out_dir, M, Q), dpi=200)
"""
print '... computed (prior) corrected posteriors'

###############################################################################################




################################################################################################
# compute the posterior of (delta_mf, delta_af)
################################################################################################
P_dmfdaf = dmf*daf*ss.correlate2d(P_mfaf_i, P_mfaf_r, boundary='fill', mode='same')
###############################################################################################

s1_dmfdaf = ba.nsigma_value(gf(P_dmfdaf))
s2_dmfdaf = ba.nsigma_value(gf(P_dmfdaf),2)


################################################################################################
# compute the posterior of (delta_mf/Mf, delta_af/af)
################################################################################################
# compute interpolation objects for the mf,af posterior and delta_mf and delta_af posterior 
P_dmfdaf_interp_object = scipy.interpolate.interp2d(mf_an[:-1], af_an[:-1], P_dmfdaf, fill_value=0., bounds_error=False)
P_mfaf_imr_interp_object = scipy.interpolate.interp2d(mf_an[:-1], af_an[:-1], P_mfaf_imr, fill_value=0., bounds_error=False)

# defining limits of delta_mf/Mf and delta_af/af. limits are currently set arbitrarily FIXME 
v1_an = np.linspace(-1.0, 1.0, 100)
v2_an = np.linspace(-1.0, 1.0, 100)

# create a matrix of zeros for allocating the posterior P(delta_mf/Mf, delta_af/af) 
P_v1v2_intrp = np.zeros([len(v1_an), len(v2_an)])

# requirements for the definition of dblquad integranttion
def integrant(af, mf):
  return P_mfaf_imr_interp_object(mf,af)*P_dmfdaf_interp_object(mf*v1, af*v2)*mf*af

mf_max = max(abs(min(mf_imr)), abs(max(mf_imr)))
af_max = max(abs(min(af_imr)), abs(max(af_imr)))

mf_an_from_imr = np.linspace(-mf_max,mf_max,1000)
af_an_from_imr = np.linspace(-af_max,af_max,1000)
dmf_an_from_imr = mf_an_from_imr[1]-mf_an_from_imr[0]
daf_an_from_imr = af_an_from_imr[1]-af_an_from_imr[0]


for i, v1 in enumerate(v1_an):
  for j, v2 in enumerate(v2_an):
    P_v1v2_intrp[i,j] = np.sum(integrant(af_an_from_imr, mf_an_from_imr))*dmf_an_from_imr*daf_an_from_imr

# save results 
#np.savetxt(out_dir+'/v1v2.dat', (v1_an,v2_an))
#np.savetxt(out_dir+'/P_v1v2.dat', P_v1v2_intrp)
#########################################################################################





#########################################################################################
# plotting										#
#########################################################################################

# details for contours
s1_i = ba.nsigma_value(gf(P_mfaf_i))
s2_i = ba.nsigma_value(gf(P_mfaf_i),2)

s1_r = ba.nsigma_value(gf(P_mfaf_r))
s2_r = ba.nsigma_value(gf(P_mfaf_r),2)

s1_imr = ba.nsigma_value(gf(P_mfaf_imr))
s2_imr = ba.nsigma_value(gf(P_mfaf_imr),2)

s1_dmfdaf = ba.nsigma_value(gf(P_dmfdaf))
s2_dmfdaf = ba.nsigma_value(gf(P_dmfdaf),2)

s1_v1v2 = ba.nsigma_value(gf(P_v1v2_intrp))
s2_v1v2 = ba.nsigma_value(gf(P_v1v2_intrp),2)

plt.figure(figsize=(16,8))
plt.subplot(231)
plt.pcolormesh(mf_an, af_an, gf(P_mfaf_i).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an[:-1], af_an[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
#plt.axvline(x=M_F, ls='--', lw=1, color='k')
#plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.xlim([M-0.2*M,M+0.2*M])
plt.ylim([0, 1.])
#plt.ylim([0, max(af_an)])
plt.title('Inspiral')
plt.grid()

plt.subplot(232)
plt.pcolormesh(mf_an, af_an, gf(P_mfaf_r).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an[:-1], af_an[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([M-0.2*M,M+0.2*M])
plt.ylim([0, 1.])
#plt.ylim([0, max(af_an)])
#plt.axvline(x=M_F, ls='--', lw=1, color='k')
#plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('Ringdown')
plt.grid()

plt.subplot(233)
plt.pcolormesh(mf_an, af_an, gf(P_mfaf_imr).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an[:-1], af_an[:-1], gf(P_mfaf_imr).T, levels=(s1_imr,s2_imr), linewidth=0.1, cmap='hot')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([M-0.2*M,M+0.2*M])
plt.ylim([0, 1.])
plt.ylim([0, max(af_an)])
#plt.axvline(x=M_F, ls='--', lw=1, color='k')
#plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.title('IMR')
plt.grid()

plt.subplot(236)
plt.contour(mf_an[:-1], af_an[:-1], gf(P_mfaf_i).T, levels=(s1_i,s2_i), linewidth=0.1, colors='c')
plt.contour(mf_an[:-1], af_an[:-1], gf(P_mfaf_r).T, levels=(s1_r,s2_r), linewidth=0.1, colors='orange')
plt.xlabel('$M_f~[M_\odot]$')
plt.ylabel('$a_f/M_f$')
plt.xlim([M-0.2*M,M+0.2*M])
#plt.ylim([0, max(af_an)])
plt.ylim([0, 1.])
#plt.axvline(x=M_F, ls='--', lw=1, color='k')
#plt.axhline(y=A_F, ls='--', lw=1, color='k')
plt.grid()

plt.subplot(234)
plt.pcolormesh(mf_an, af_an, gf(P_dmfdaf).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(mf_an[:-1], af_an[:-1], gf(P_dmfdaf).T, levels=(s1_dmfdaf,s2_dmfdaf), linewidth=0.1, cmap='hot')
plt.axvline(x=0, ls='--', lw=1, color='k')
plt.axhline(y=0, ls='--', lw=1, color='k')
plt.xlabel('$\Delta M_f~[M_\odot]$')
plt.ylabel('$\Delta a_f$')
plt.xlim([-50.,50.])
plt.ylim([-1.,1.])
plt.grid()

plt.subplot(235)
plt.pcolormesh(v1_an,v2_an,(P_v1v2_intrp).T, cmap='YlOrBr')
plt.colorbar()
plt.contour(v1_an, v2_an, (P_v1v2_intrp).T, levels=(s1_v1v2,s2_v1v2), linewidth=0.1, cmap='hot')
plt.axvline(x=0, ls='--', lw=1, color='k')
plt.axhline(y=0, ls='--', lw=1, color='k')
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')
plt.xlim([-1.,1.])
plt.ylim([-1.,1.])
plt.grid()
#plt.suptitle('$M = %2.1f~M_\odot,~q = %2.1f$' %(M, Q))
#plt.savefig('%s/TestGRPlots_M%2.1f_q%2.1f.png' %(out_dir, M, Q))
plt.savefig('../test/NEW_gaussian_analysis/fig_1.png')
exit()

print '... plotted final plot' 
#########################################################################################

#debug plots

if debug_plots == True: 
  x1, x2 = dp.debug_plots(m1_i, m2_i, logl_i, 'inspiral', M, Q, M_F, A_F, out_dir)
  x1, x2 = dp.debug_plots(m1_r, m2_r, logl_r, 'ringdown', M, Q, M_F, A_F, out_dir)
  x1, x2 = dp.debug_plots(m1_r, m2_r, logl_imr, 'IMR', M, Q, M_F, A_F, out_dir)
  print '... made debug plots' 

print '... complted. time taken = %2.f secs' %(time.time()-start_time)
