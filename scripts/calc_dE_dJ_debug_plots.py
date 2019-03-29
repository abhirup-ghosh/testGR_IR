import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import loadline as ll
import matplotlib
from matplotlib import cm
import plotsettings
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize


def debug_plots(m1, m2, logl, phase, M, Q, M_F, A_F, out_dir):
        """
        Calculate the 1) debug plots in mtot-eta, mtot-q and af-mf (consisting of likelihood, un-interpolated prior, lalinference posterior, corrected posterior) and,
        2) un-interpolated prior, and interpolated prior
        """
	mtot = m1 + m2
	eta = m1*m2/mtot**2
	q = m2/m1
	mf = mtot*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.))
        af = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)


	#posteriors
	P_eta_mtot, eta_bins, mtot_bins = np.histogram2d(eta, mtot, bins=50, normed=True)
	P_q_mtot, q_bins, mtot_bins = np.histogram2d(q, mtot, bins=50, normed=True)
	P_mf_af, mf_bins, af_bins = np.histogram2d(mf, af, bins=50, normed=True)

	#prior
	comp_mass_min = min([min(m1), min(m2)])
	comp_mass_max = max([max(m1), max(m2)])

	m1_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)
	m2_pr = np.random.uniform(comp_mass_min, comp_mass_max, 100000)
	
	mtot_pr = m1_pr + m2_pr
        eta_pr = m1_pr*m2_pr/mtot_pr**2
        q_pr = m2_pr/m1_pr
        mf_pr = mtot_pr*(1. + (np.sqrt(8./9.)-1.)*eta_pr - 0.4333*(eta_pr**2.) - 0.4392*(eta_pr**3.))
        af_pr = eta_pr*np.sqrt(12.) - 3.871*(eta_pr**2.) + 4.028*(eta_pr**3.)


	P_eta_mtot_pr, eta_bins_pr, mtot_bins_pr = np.histogram2d(eta_pr, mtot_pr, bins=50, normed=True)
        P_q_mtot_pr, q_bins_pr, mtot_bins_pr = np.histogram2d(q_pr, mtot_pr, bins=50, normed=True)
        P_mf_af_pr, mf_bins_pr, af_bins_pr = np.histogram2d(mf_pr, af_pr, bins=50, normed=True)
	
	#interpolated prior
	P_eta_mtot_pr_intrp_objct = scipy.interpolate.RectBivariateSpline(eta_bins_pr[:-1], mtot_bins_pr[:-1], P_eta_mtot_pr)	
	P_q_mtot_pr_intrp_objct = scipy.interpolate.RectBivariateSpline(q_bins_pr[:-1], mtot_bins_pr[:-1], P_q_mtot_pr)	
	P_mf_af_pr_intrp_objct = scipy.interpolate.RectBivariateSpline(mf_bins_pr[:-1], af_bins_pr[:-1], P_mf_af_pr)	

	P_eta_mtot_pr_intrp = P_eta_mtot_pr_intrp_objct(eta_bins[:-1], mtot_bins[:-1])
	P_q_mtot_pr_intrp = P_q_mtot_pr_intrp_objct(q_bins[:-1], mtot_bins[:-1])
	P_mf_af_pr_intrp = P_mf_af_pr_intrp_objct(mf_bins[:-1], af_bins[:-1])

	#corrected posteriors
	P_eta_mtot_corr = P_eta_mtot/P_eta_mtot_pr_intrp
	P_q_mtot_corr = P_q_mtot/P_q_mtot_pr_intrp
	P_mf_af_corr = P_mf_af/P_mf_af_pr_intrp

	#plotting
	plt.figure(figsize=(20,16))
	"""
	plt.subplot(341)
	plt.title('Raw Posteriors')
	plt.pcolormesh(eta_bins, mtot_bins, np.log(P_eta_mtot.T), cmap='YlOrBr')
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(eta_bins), max(eta_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
	
	plt.subplot(342)
	plt.title('Interpolated Prior')
	plt.pcolormesh(eta_bins, mtot_bins, np.log(P_eta_mtot_pr_intrp.T), cmap='YlOrBr')
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(eta_bins), max(eta_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
	
	plt.subplot(343)
	plt.title('Likelihood')
	#plt.scatter(eta, mtot, c=logl, norm = matplotlib.colors.Normalize(vmax = max(logl), vmin = max(logl)-1), edgecolor='', cmap='YlOrBr')
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(eta), max(eta)]); plt.ylim([min(mtot), max(mtot)]);
	plt.grid()
	plt.subplot(344)
	
	plt.title('Corrected Posteriors')
	plt.pcolormesh(eta_bins, mtot_bins, np.log(P_eta_mtot_corr.T), cmap='YlOrBr')
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(eta_bins), max(eta_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
	
	plt.subplot(345)
        plt.pcolormesh(q_bins, mtot_bins, np.log(P_q_mtot.T), cmap='YlOrBr')
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(q_bins), max(q_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
        
	plt.subplot(346)
        plt.pcolormesh(q_bins, mtot_bins, np.log(P_q_mtot_pr_intrp.T), cmap='YlOrBr')
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(q_bins), max(q_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
        
	plt.subplot(347)
        #plt.scatter(q, mtot, c=logl, norm = matplotlib.colors.Normalize(vmax = max(logl), vmin = max(logl)-1), edgecolor='', cmap='YlOrBr')
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(q), max(q)]); plt.ylim([min(mtot), max(mtot)]);
	plt.grid()
        
	plt.subplot(348)
        plt.pcolormesh(q_bins, mtot_bins, np.log(P_q_mtot_corr.T), cmap='YlOrBr')
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.xlim([min(q_bins), max(q_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.grid()
	"""
	plt.subplot(349)
        plt.pcolormesh(mf_bins, af_bins, P_mf_af.T, norm=matplotlib.colors.Normalize(vmin=(P_mf_af.T).min(), vmax=(P_mf_af.T).max()), cmap='YlOrBr')
	plt.axvline(x=M_F, ls='--', lw=1, color='k')
	plt.axhline(y=A_F, ls='--', lw=1, color='k')
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.ylim([min(af_bins), max(af_bins)]); plt.xlim([min(mf_bins), max(mf_bins)]);
        
	plt.subplot(3,4,10)
        plt.pcolormesh(mf_bins, af_bins, P_mf_af_pr_intrp.T, norm=matplotlib.colors.Normalize(vmin=(P_mf_af_pr_intrp.T).min(), vmax=(P_mf_af_pr_intrp.T).max()), cmap='YlOrBr')
	plt.axvline(x=M_F, ls='--', lw=1, color='k')
	plt.axhline(y=A_F, ls='--', lw=1, color='k')
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.ylim([min(af_bins), max(af_bins)]); plt.xlim([min(mf_bins), max(mf_bins)]);
	plt.grid()
        
	plt.subplot(3,4,11)		
        plt.scatter(mf, af, c=logl, norm = matplotlib.colors.Normalize(vmax = max(logl), vmin = max(logl)-1), edgecolor='', cmap='YlOrBr')
	plt.axvline(x=M_F, ls='--', lw=1, color='k')
        plt.axhline(y=A_F, ls='--', lw=1, color='k')
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.ylim([min(af), max(af)]); plt.xlim([min(mf), max(mf)]);
	plt.grid()
        
	plt.subplot(3,4,12)
        plt.pcolormesh(mf_bins, af_bins, P_mf_af_corr.T, norm=matplotlib.colors.Normalize(vmin=(P_mf_af_corr.T).min(), vmax=(P_mf_af_corr.T).max()), cmap='YlOrBr')
	plt.axvline(x=M_F, ls='--', lw=1, color='k')
	plt.axhline(y=A_F, ls='--', lw=1, color='k')
	plt.ylim([min(af_bins), max(af_bins)]); plt.xlim([min(mf_bins), max(mf_bins)]);
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.grid()
	plt.tight_layout()		
	x1 = plt.savefig(out_dir+'/debug_plots_%s_%2.1f_%2.1f.png'%(phase, M, Q))

	


	plt.figure(figsize=(6,8))
	plt.subplot(321)
	plt.title('Uninterpolated Prior')
	plt.pcolormesh(eta_bins_pr, mtot_bins_pr, P_eta_mtot_pr.T, cmap='YlOrBr')
	plt.xlim([min(eta_bins), max(eta_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.grid()
	
	plt.subplot(322)
	plt.title('Interpolated Prior')
	plt.pcolormesh(eta_bins, mtot_bins, P_eta_mtot_pr_intrp.T, cmap='YlOrBr')
	plt.xlim([min(eta_bins), max(eta_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.xlabel('$\eta$'); plt.ylabel('$M_{tot}$');
	plt.grid()
	
	plt.subplot(323)
        plt.pcolormesh(q_bins_pr, mtot_bins_pr, P_q_mtot_pr.T, cmap='YlOrBr')
	plt.xlim([min(q_bins), max(q_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.grid()
        
	plt.subplot(324)
        plt.pcolormesh(q_bins, mtot_bins, P_q_mtot_pr_intrp.T, cmap='YlOrBr')
	plt.xlim([min(q_bins), max(q_bins)]); plt.ylim([min(mtot_bins), max(mtot_bins)]);
	plt.xlabel('$q$'); plt.ylabel('$M_{tot}$');
	plt.grid()
	
	plt.subplot(325)
        plt.pcolormesh(mf_bins_pr, af_bins_pr, P_mf_af_pr.T, cmap='YlOrBr')
	plt.ylim([min(af_bins), max(af_bins)]); plt.xlim([min(mf_bins), max(mf_bins)]);
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.grid()
        
	plt.subplot(326)
        plt.pcolormesh(mf_bins, af_bins, P_mf_af_pr_intrp.T, cmap='YlOrBr')
	plt.ylim([min(af_bins), max(af_bins)]); plt.xlim([min(mf_bins), max(mf_bins)]);
	plt.ylabel('$a_f/M_f$'); plt.xlabel('$M_f$');
	plt.grid()
	plt.tight_layout()
	x2 = plt.savefig(out_dir+'/prior_plots_%s_%2.1f_%2.1f.png'%(phase, M, Q))
	
	return x1, x2
