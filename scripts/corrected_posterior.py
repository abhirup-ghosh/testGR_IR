import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import compute_rad_energy_from_waveforms as crefw
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import os
import plotsettings
import commands
import standard_gwtransf as gw
import my_plotter as mp
import loadline as ll
import scipy
from scipy import interpolate

date = '2015-06-05'
approximant = 'IMRPhenomB'

phase_list = ['inspiral', 'ringdown']
M_list = [30., 50., 100.]
Q_list = [1., 2., 4.]
n_isco_list = [1.0]



# Loading Data

for phase in phase_list:

  # Output arguments
  out_dir = '../plots/%s_%s/%s/corrected_posteriors'%(date, approximant, phase)
  os.system('mkdir -p %s'%(out_dir))

  for M in M_list:
    for Q in Q_list:
        for n_isco in n_isco_list:
                M1 = M/(1.+Q)
                M2 = M*Q/(1.+Q)
                ETA = (M1*M2)/((M1+M2)**2.)

		MF = M*(1. + (np.sqrt(8./9.)-1.)*ETA - 0.4333*(ETA**2.) - 0.4392*(ETA**3.))
        	AF = ETA*np.sqrt(12.) - 3.871*(ETA**2.) + 4.028*(ETA**3.)


        	data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase, str(M), str(Q), str(n_isco))
        	legend = ll.load_line(data_file, dtype='string')

        	mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)

        	mf = mtotal*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.))
       		af = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)

		
		plt.figure(figsize=(16,16))


		nbins = 50


		#POSTERIOR
		P_mf_af, afbins, mfbins = np.histogram2d(af,mf,bins=nbins, normed=True)
		P_mf_af_maskd = np.log(np.ma.masked_where(P_mf_af==0,P_mf_af))
		#xcenters1 = xedges1[:-1] + 0.5*(xedges1[1:] - xedges1[:-1])
		#ycenters1 = yedges1[:-1] + 0.5*(yedges1[1:] - yedges1[:-1])

		plt.subplot(231)
		plt.pcolormesh(afbins, mfbins, P_mf_af_maskd.T)
		#plt.pcolormesh(afbins, mfbins, P_mf_af.T)
		plt.colorbar()
		plt.axhline(y=MF, ls='--', lw=0.5)
      	 	plt.axvline(x=AF, ls='--', lw=0.5)
		plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
		plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
		plt.grid()		
                plt.title('posterior')
		

		#PRIOR
		comp_mass_min = min([min(m1), min(m2)])
		comp_mass_max = max([max(m1), max(m2)])

		m1_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)
		m2_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)

                mtotal_pr =  m1_pr + m2_pr
                eta_pr = (m1_pr*m2_pr)/((m1_pr+m2_pr)**2.)

                mf_pr = mtotal_pr*(1. + (np.sqrt(8./9.)-1.)*eta_pr - 0.4333*(eta_pr**2.) - 0.4392*(eta_pr**3.))
                af_pr = eta_pr*np.sqrt(12.) - 3.871*(eta_pr**2.) + 4.028*(eta_pr**3.)

                P_mf_af_pr, afbins_pr, mfbins_pr = np.histogram2d(af_pr,mf_pr,bins=75, normed=True)
		P_mf_af_pr_maskd = np.log(np.ma.masked_where(P_mf_af_pr==0,P_mf_af_pr))
		#xcenters2 = xedges2[:-1] + 0.5*(xedges2[1:] - xedges2[:-1])
                #ycenters2 = yedges2[:-1] + 0.5*(yedges2[1:] - yedges2[:-1])
		
		plt.subplot(235)
		plt.pcolormesh(afbins_pr, mfbins_pr, P_mf_af_pr_maskd.T)
		#plt.pcolormesh(afbins_pr, mfbins_pr, P_mf_af_pr.T)
		plt.colorbar()
		plt.axhline(y=MF, ls='--', lw=0.5)
        	plt.axvline(x=AF, ls='--', lw=0.5)
		plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
		plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
                plt.title('uninterpolated prior')
		

		#LIKELIHOOD
		plt.subplot(233)
		plt.scatter(af, mf, c=logl, norm = matplotlib.colors.Normalize(vmax = max(logl), vmin = max(logl)-1), edgecolor='')
		plt.colorbar()
		plt.axhline(y=MF, ls='--', lw=0.5)
        	plt.axvline(x=AF, ls='--', lw=0.5)
		plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
		plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
		plt.grid()		
                plt.title('likelihood')

		

		#INTERPOLATED PRIOR
		P_mf_af_pr_intrp_obj = scipy.interpolate.RectBivariateSpline(afbins_pr[:-1], mfbins_pr[:-1], P_mf_af_pr)
		P_mf_af_pr_intrp = P_mf_af_pr_intrp_obj(afbins[:-1], mfbins[:-1])
		P_mf_af_pr_intrp_maskd = np.log(np.ma.masked_where(P_mf_af_pr_intrp < np.exp(P_mf_af_pr_maskd.min()) ,P_mf_af_pr_intrp))	


		plt.subplot(232)
		plt.pcolormesh(afbins, mfbins, P_mf_af_pr_intrp_maskd.T)
		#plt.pcolormesh(afbins, mfbins, P_mf_af_pr_intrp.T)
                plt.colorbar()
                plt.axhline(y=MF, ls='--', lw=0.5)
                plt.axvline(x=AF, ls='--', lw=0.5)
                plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
                plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
		plt.grid()		
                plt.title('interpolated prior')


		#CORRECTED posterior 
		P_mf_af_corr = P_mf_af/P_mf_af_pr_intrp	
		P_mf_af_corr_maskd = np.log(np.ma.masked_where(P_mf_af_corr==0,P_mf_af_corr))	

		plt.subplot(234)
		plt.pcolormesh(afbins,mfbins, P_mf_af_corr_maskd.T)
		#plt.pcolormesh(afbins,mfbins, P_mf_af_corr.T)
		plt.colorbar()
		plt.axhline(y=MF, ls='--', lw=0.5)
        	plt.axvline(x=AF, ls='--', lw=0.5)
		plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
		plt.xlabel('$a_f/M_f$')
	        plt.ylabel('$M_f$')
        	plt.title('corrected posterior')
		plt.grid()		

		plt.suptitle('M = %s; q = %s; n = %s'%(str(M),str(Q),str(n_isco)))
		plt.tight_layout()
		plt.savefig('%s/posteriors_%s_%s_%s.png'%(out_dir, str(M), str(Q), str(n_isco)))
	        plt.close()

		"""
		plt.figure(figsize=(16,5))
		plt.subplot(121)
                plt.pcolormesh(afbins_pr, mfbins_pr, P_mf_af_pr_maskd.T)
                #plt.pcolormesh(afbins_pr, mfbins_pr, P_mf_af_pr.T)
                plt.colorbar()
                plt.axhline(y=MF, ls='--', lw=0.5)
                plt.axvline(x=AF, ls='--', lw=0.5)
                plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
                plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
                plt.title('uninterpolated prior')

		plt.subplot(122)
                plt.pcolormesh(afbins, mfbins, P_mf_af_pr_intrp_maskd.T)
                #plt.pcolormesh(afbins, mfbins, P_mf_af_pr_intrp.T)
                plt.colorbar()
                plt.axhline(y=MF, ls='--', lw=0.5)
                plt.axvline(x=AF, ls='--', lw=0.5)
                plt.xlim([min(af), max(af)])
                plt.ylim([min(mf), max(mf)])
                plt.xlabel('$a_f/M_f$')
                plt.ylabel('$M_f$')
                plt.grid()
                plt.title('interpolated prior')
		"""
#plt.show()
