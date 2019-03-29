""" 
Compute the posterior on the parameters deltaE and deltaJ using the posteriors 
from the lalinference runs using inspiral and ringdown templates. deltaE and deltaJ
describe the deviation from GR in the radiated energy and angular momentum. 

Abhirup Ghosh and P. Ajith, 2015-04-08 

"""
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
	

import numpy as np
import loadline as ll
import scipy
import scipy.ndimage.filters as filter
import scipy.signal as ss 
from optparse import OptionParser
import plotsettings 
import bayesian as ba
import calc_dE_dJ_debug_plots as dp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# read inputs from command line 
parser = OptionParser()
parser.add_option("-i", "--insp-post", dest="insp_post",
                  help="file containing the posterior samples from the lalinference inspiral run")
parser.add_option("-r", "--ring-post", dest="ring_post",
                  help="file containing the posterior samples from the lalinference ringdown run")
parser.add_option("-f", "--fit-formula", dest="fit_formula",
                  help="fitting formula to be used for the calculation of final mass/spin [options: 'nospin_pan2011'")
parser.add_option("-d", "--debug-plots", dest="debug_plots", help="")
#parser.add_option("--savefig", dest="savefig", type="string",help="savefig")

(options, args) = parser.parse_args()

insp_post = options.insp_post
ring_post = options.ring_post
fit_formula = options.fit_formula
debug_plots = options.debug_plots
#savefig = options.savefig

# read the posterior samples -- inspiral 
legend = ll.load_line(insp_post, dtype='string')
m1_i, m2_i, logl_i = np.loadtxt(insp_post, skiprows=1, usecols = (list(legend).index('m1'),\
	list(legend).index('m2'), list(legend).index('logl')),unpack=True)

# read the posterior samples -- ringdown  
legend = ll.load_line(ring_post, dtype='string')
m1_r, m2_r, logl_r = np.loadtxt(ring_post, skiprows=1, usecols = (list(legend).index('m1'),\
	list(legend).index('m2'), list(legend).index('logl')),unpack=True)

# compute the total mass and mass ratio of the posterior samples 
mtot_i = m1_i+m2_i 
mtot_r = m1_r+m2_r 
eta_i = m1_i*m2_i/mtot_i**2.
eta_r = m1_r*m2_r/mtot_r**2.

# compute mass and spin of the final black hole 
mf_i, af_i = calc_final_mass_spin(mtot_i, eta_i, fit_formula) 
mf_r, af_r = calc_final_mass_spin(mtot_r, eta_r, fit_formula) 

# compute the limits of the correlation integral 
m_max = max(np.append(abs(mf_i), abs(mf_r)))
a_max = max(np.append(abs(af_i), abs(af_r)))

print 'm_max = ', m_max, ', a_max = ', a_max 

nm_bins = 100
na_bins = 100
mbins = np.linspace(-m_max, m_max, nm_bins)
abins = np.linspace(-a_max, a_max, na_bins)

# resolution of the 2D posteriors in the Mf and af space 
dm = np.mean(np.diff(mbins))
da = np.mean(np.diff(abins))

# compute the minimum and maximum value of component masses sampled by the stochastic sampler 
comp_mass_min = min([min(m1_i), min(m2_i), min(m1_r), min(m2_r)]) 
comp_mass_max = max([max(m1_i), max(m2_i), max(m1_r), max(m2_r)])

print 'comp_mass_min = ', comp_mass_min, 'comp_mass_max = ', comp_mass_max

#prior
m1_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)
m2_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)

mtot_pr =  m1_pr + m2_pr
eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)

mf_pr, af_pr = calc_final_mass_spin(mtot_pr, eta_pr, fit_formula)

P_mfaf_pr, abins_pr, mbins_pr = np.histogram2d(af_pr,mf_pr,bins=100, normed=True)
P_mfaf_interp_objct = scipy.interpolate.RectBivariateSpline(abins_pr[:-1], mbins_pr[:-1], P_mfaf_pr)
H4 = P_mfaf_interp_objct(abins[:-1], mbins[:-1])
#----------------------------------------------------------------------------------------



# compute the 2D posteriors in Mf and af 
P_mfaf_i, abins, mbins = np.histogram2d(af_i, mf_i, bins=(abins, mbins))
P_mfaf_r, abins, mbins = np.histogram2d(af_r, mf_r, bins=(abins, mbins))




#---------------------------------------------------------------------------------
# compute the corrected 2D posteriors in Mf and af
P_mfaf_i = P_mfaf_i/H4
P_mfaf_r = P_mfaf_r/H4
#---------------------------------------------------------------------------------



P_mfaf_i = P_mfaf_i.T 
P_mfaf_r = P_mfaf_r.T

# compute the posteiors in deltaE and deltaJ, where 
# deltaE = mf_i - mf_r 
# deltaJ = af_i/mf_i - af_r/mf_r 
P_deltaE_deltaJ = dm*da*ss.correlate2d(P_mfaf_i, P_mfaf_r, boundary='fill', mode='same')
P_deltaE_deltaJ = filter.gaussian_filter(P_deltaE_deltaJ, sigma=1.0)

# details for contours
#s1 = ba.nsigma_value(P_deltaE_deltaJ)
#s2 = ba.nsigma_value(P_deltaE_deltaJ,2)

#print s1, s2

# plot the posteriors 
plt.figure(figsize=(20,6))
plt.subplot(131)
plt.title('Inspiral')
#plt.hist2d(af_i, mf_i, bins=30, cmap='YlOrBr')
plt.pcolormesh(abins[:-1], mbins[:-1], P_mfaf_i, cmap='YlOrBr')
plt.ylim([min(mf_i), max(mf_i)])
plt.xlim([min(af_i), max(af_i)])
#plt.ylim([mbins[0], mbins[-1]])
#plt.xlim([abins[0], abins[-1]])
#plt.imshow(P_mfaf_i, cmap='YlOrBr', extent=(m_min, m_max, a_min, a_max), aspect='auto', interpolation='gaussian')
plt.ylabel('$M_f~[M_\odot]$')
plt.xlabel('$a_f/M_f$')
plt.grid()
plt.subplot(132)
plt.title('Ringdown')
#plt.hist2d(af_r, mf_r, bins=30, cmap='YlOrBr')
plt.pcolormesh(abins[:-1], mbins[:-1], P_mfaf_r, cmap='YlOrBr')
plt.ylim([min(mf_r), max(mf_r)])
plt.xlim([min(af_r), max(af_r)])
#plt.ylim([mbins[0], mbins[-1]])
#plt.xlim([abins[0], abins[-1]])
#plt.imshow(P_mfaf_r, cmap='YlOrBr', extent=(m_min, m_max, a_min, a_max), aspect='auto', interpolation='gaussian')
plt.ylabel('$M_f~[M_\odot]$')
plt.xlabel('$a_f/M_f$')
plt.grid()
plt.subplot(133)
plt.title('Deviation from GR')
plt.pcolormesh(abins[:-1], mbins[:-1], P_deltaE_deltaJ, cmap='YlOrBr')
#plt.imshow(P_deltaE_deltaJ, cmap='YlOrBr', extent=(m_min, m_max, a_min, a_max), aspect='auto', interpolation='gaussian')
#plt.contour(abins[:-1], mbins[:-1], P_deltaE_deltaJ, levels=(s1,s2), linewidth=0.5, cmap='hot')
plt.ylim(-20,20)
plt.xlim(-0.5, 0.5)
plt.ylabel('$\Delta E~ [M_\odot]$')
plt.xlabel('$\Delta J$')
plt.grid()
plt.tight_layout()
plt.savefig('mf_af_contours_%s_%s_%s.png'%(str(M), str(Q), str(n)), dpi=150)


#debug plots

if debug_plots == 'y':
	x1, x2 = dp.debug_plots(m1_i, m2_i, logl_i, 'inspiral', M, Q, n)
	x1, x2 = dp.debug_plots(m1_r, m2_r, logl_r, 'ringdown', M, Q, n)
	
