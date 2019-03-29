import numpy as np
import matplotlib.pyplot as plt
import bayesian as ba
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import standard_gwtransf as gw
import commands
import loadline as ll
import plotsettings 
import scipy.ndimage.filters as filter
import scipy
from scipy import interpolate
import os

date = '2015-05-06'
approximant = 'IMRPhenomB'

phase_list = ['inspiral', 'ringdown']


M_list = [30., 50., 100.]
Q_list = [1., 2., 4.]
n_isco_list = [1.0, 1.0]#[1.0, 0.5]
#n_qnm_list = [1.0]
        
nbins = 30
color = np.genfromtxt('color.dat', dtype='str')
# Output arguments
out_dir = '../plots/%s_%s/contour_corrected'%(date, approximant)
os.system('mkdir -p %s'%(out_dir))


for M in M_list:
  for Q in Q_list:
    if Q == 1.:
        f_qnm = 0.0879
    elif Q == 2.:
        f_qnm = 0.0829
    elif Q == 4.:
        f_qnm = 0.0744

    plt.figure(figsize=(6,6))

    for ((i, phase), n_isco) in zip(enumerate(phase_list), n_isco_list):
      #for n_isco in n_isco_list:		
	M1 = M/(1.+Q)
        M2 = M*Q/(1.+Q)
        ETA = gw.eta_from_comp(M1,M2)

        MF = M*(1. + (np.sqrt(8./9.)-1.)*ETA - 0.4333*(ETA**2.) - 0.4392*(ETA**3.))
        AF_MF = ETA*np.sqrt(12.) - 3.871*(ETA**2.) + 4.028*(ETA**3.)

        data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase, str(M), str(Q), str(n_isco))
	legend = ll.load_line(data_file, dtype='string')

        mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)
		
        mf = mtotal*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.))
        af_mf = (eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.))*(mf**2.)

	print min(mf), max(mf)
	print min(af_mf), max(af_mf)
	nbins = 50


        #POSTERIOR
        H1, xedges1, yedges1 = np.histogram2d(mf,af_mf,bins=nbins, normed=True)
        Hmasked1 = np.log(np.ma.masked_where(H1==0,H1))
        xcenters1 = xedges1[:-1] + 0.5*(xedges1[1:] - xedges1[:-1])
        ycenters1 = yedges1[:-1] + 0.5*(yedges1[1:] - yedges1[:-1])

	#PRIOR
        m1_pr = np.random.uniform(0.25*M1, 2.*M2, 100000)
        m2_pr = np.random.uniform(0.25*M1, 2.*M2, 100000)

        mtotal_pr =  m1_pr + m2_pr
        eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)

        mf_pr = mtotal_pr*(1. + (np.sqrt(8./9.)-1.)*eta_pr - 0.4333*(eta_pr**2.) - 0.4392*(eta_pr**3.))
        af_mf_pr = eta_pr*np.sqrt(12.) - 3.871*(eta_pr**2.) + 4.028*(eta_pr**3.)

        H2, xedges2, yedges2 = np.histogram2d(mf_pr,af_mf_pr,bins=75, normed=True)
        Hmasked2 = np.log(np.ma.masked_where(H2==0,H2))
        xcenters2 = xedges2[:-1] + 0.5*(xedges2[1:] - xedges2[:-1])
        ycenters2 = yedges2[:-1] + 0.5*(yedges2[1:] - yedges2[:-1])

	#INTERPOLATED PRIOR
        H3 = scipy.interpolate.RectBivariateSpline(xcenters2, ycenters2, H2)

        H4 = H3(xcenters1, ycenters1)
        Hmasked4 = np.log(np.ma.masked_where(H4==0,H4))

	#CORRECTED PRIOR
	H5 = H1/H4
	H5 = filter.gaussian_filter(H5, sigma=1.0)	

	s1 = ba.nsigma_value(H5)
	s2 = ba.nsigma_value(H5)

        # Mask zeros
        Hmasked5 = np.ma.masked_where(H5==0,H5) # Mask pixels with a value o

        plt.contourf(ycenters1, xcenters1, H5, levels=(s1,s2), colors=color[i], alpha=0.3)
	plt.axhline(y=MF, ls='--', lw=0.5)
        plt.axvline(x=AF_MF, ls='--', lw=0.5)
	plt.xlabel('$a_f/M_f$')
        plt.ylabel('$M_f$')
	plt.legend(frameon=False, loc='best')
        plt.hold(True)
	print phase
    #plt.suptitle('$f_\mathrm{ISCO}= %2.1f$ Hz, $f_\mathrm{QNM} = %2.1f$\n$M = %2.1f~M_\odot,~q = %2.1f$ \n Contours plotted: 1.0 $f_{isco}$ and 0.25 $f_{qnm}$' %((1./6.)**(3./2)/(LAL_PI*M*LAL_MTSUN_SI), f_qnm/(M*LAL_MTSUN_SI), M, Q))
    plt.savefig('%s/contour_overlap_%s_%s.png'%(out_dir, M, Q))
