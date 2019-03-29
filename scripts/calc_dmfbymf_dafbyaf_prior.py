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

def calc_P_Mfmean_afmean(Mf_vec, af_vec, Mf_min, Mf_max):
	""" compute the distribution of mean_Mf and mean_af """

	P = np.zeros((len(Mf_vec), len(af_vec)))
	for iMf, Mf in enumerate(Mf_vec):
		dx = 2.*Mf 
		for iaf, af in enumerate(af_vec): 
			dy = 2.*af 
			if dx > Mf_max + Mf_min and dx < 2*Mf_max and 1 < dy < 2:
				P[iMf, iaf] = (-2 + dy)*(dx - 2*Mf_max)
			elif dx > 2*Mf_min and dx <= Mf_max + Mf_min and 1 < dy < 2:
				P[iMf, iaf] = -(-2 + dy)*(dx - 2*Mf_min)
			elif dx > Mf_max + Mf_min and dx < 2*Mf_max and 0 < dy <= 1:
				P[iMf, iaf] = -dy*(dx - 2*Mf_max)
			elif dx > 2*Mf_min and dx <= Mf_max + Mf_min and 0 < dy <= 1:
				P[iMf, iaf] = dy*(dx - 2*Mf_min)
			else:
				P[iMf, iaf] = 0. 
	return P 

def calc_P_dMf_daf(dMf_vec, daf_vec, Mf_min, Mf_max):
	""" compute the distribution of delta_Mf and delta_af """

	P = np.zeros((len(dMf_vec), len(daf_vec)))
	for iMf, dx in enumerate(dMf_vec):
		for iaf, dy in enumerate(daf_vec): 

			if dx > 0 and Mf_max > dx + Mf_min and 0 < dy < 1: 
				P[iMf, iaf] =	-(-1 + dy)*(-dx + Mf_max - Mf_min)
			elif dx + Mf_max > Mf_min and dx <= 0 and 0 < dy < 1: 
				P[iMf, iaf] = -(-1 + dy)*(dx + Mf_max - Mf_min)
			elif dx > 0 and Mf_max > dx + Mf_min and -1 < dy <= 0:
				P[iMf, iaf] = (1 + dy)*(-dx + Mf_max - Mf_min)
			elif dx + Mf_max > Mf_min and dx <= 0 and -1 < dy <= 0:
				P[iMf, iaf] = (1 + dy)*(dx + Mf_max - Mf_min)
			else: 
				P[iMf, iaf] = 0. 

	return P 

if __name__ == '__main__':

	Mf_min, Mf_max = 1., 600. 
	dMfbyMf_lim = 2. 
	dafbyaf_lim = 2. 
	Nbins = 101 
	Mf_vec = np.linspace(-2*Mf_max, 2*Mf_max, Nbins)
	af_vec = np.linspace(-2., 2., Nbins) 


	# evaluate the posteriors of the mean and difference 
	P_Mfmean_afmean = calc_P_Mfmean_afmean(Mf_vec, af_vec, Mf_min, Mf_max)
	P_dMf_daf = calc_P_dMf_daf(Mf_vec, af_vec, Mf_min, Mf_max)
	
	# normalize the posteriors 
	dM = (Mf_vec[-1]-Mf_vec[0])/Nbins
	da = (af_vec[-1]-af_vec[0])/Nbins
	P_Mfmean_afmean /= np.sum(P_Mfmean_afmean)*dM*da 
	P_dMf_daf /= np.sum(P_dMf_daf)*dM*da

	# create interpolation objects 
	P_Mfmean_afmean_interp_obj = interp.interp2d(Mf_vec, af_vec, P_Mfmean_afmean, fill_value=0., bounds_error=False)
	P_dMf_daf_interp_obj = interp.interp2d(Mf_vec, af_vec, P_dMf_daf, fill_value=0., bounds_error=False)

	# evaluate the interpolation objects, find the difference with the original data 		
	P_Mfmean_afmean_interp = P_Mfmean_afmean_interp_obj(Mf_vec, af_vec)	
	P_dMf_daf_interp = P_dMf_daf_interp_obj(Mf_vec, af_vec)	

	# defining limits of delta_Mf/Mf and delta_af/af
	dMfbyMf_vec = np.linspace(-dMfbyMf_lim, dMfbyMf_lim, Nbins)
	dafbyaf_vec = np.linspace(-dafbyaf_lim, dafbyaf_lim, Nbins)

	# compute the P(dMf/Mf, daf/af) by evaluating the integral
	diff_dMfbyMf = np.mean(np.diff(dMfbyMf_vec))
	diff_dafbyaf = np.mean(np.diff(dafbyaf_vec))
	P_dMfbyMf_dafbyaf = np.zeros(shape=(Nbins,Nbins))

	# compute the posterior on the fractional deviation parameters (delta_Mf/Mf, delta_af/af). 
	# Approximate the integral in Eq.(6) of the document LIGO-P1500185-v5 by a discrete sum
	for i, v2 in enumerate(dafbyaf_vec):
			for j, v1 in enumerate(dMfbyMf_vec):
				P_dMfbyMf_dafbyaf[i,j] = tgr.calc_sum(Mf_vec, af_vec, v1, v2, P_dMf_daf_interp_obj, P_Mfmean_afmean_interp_obj)

	# normalization 
	P_dMfbyMf_dafbyaf /= np.sum(P_dMfbyMf_dafbyaf) * diff_dMfbyMf * diff_dafbyaf

	# create interpolation objects, evaluate the interpolation objects, find the difference with the original data 		
	P_dMfbyMf_dafbyaf_interp_obj = interp.interp2d(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMf_dafbyaf, fill_value=0., bounds_error=False)
	P_dMfbyMf_dafbyaf_interp = P_dMfbyMf_dafbyaf_interp_obj(dMfbyMf_vec, dafbyaf_vec)	

	print '... maximum errors due to interpolation are: %3.2e %3.2e %3.2e' %(np.max(abs(P_Mfmean_afmean-P_Mfmean_afmean_interp)), np.max(abs(P_dMf_daf-P_dMf_daf_interp)), np.max(abs(P_dMfbyMf_dafbyaf-P_dMfbyMf_dafbyaf_interp)))

	# save teh interpolation object 
	outfile = 'SemiAnalyticPrior_dMfbyMf_dafbyaf_Mfmin_%.1f_Mfmax_%.1f_afmin0_afmax1_dMfbyMflim%.1f_dafbyaflim%.1f'%(Mf_min, Mf_max, dMfbyMf_lim, dafbyaf_lim)
	f = gzip.open(outfile+".pklz",'wb')
	pickle.dump(P_dMfbyMf_dafbyaf_interp_obj, f)

	plt.figure(figsize=(14,8))
	plt.subplot(231)
	plt.pcolormesh(Mf_vec, af_vec, P_Mfmean_afmean, cmap='YlOrBr')
	plt.contour(Mf_vec, af_vec, P_Mfmean_afmean, cmap='YlOrBr')
	plt.xlabel('$M_f$')
	plt.ylabel('$a_f$')
	plt.grid()
	plt.colorbar()
	plt.title('Mean')

	plt.subplot(234)
	plt.pcolormesh(Mf_vec, af_vec, P_Mfmean_afmean_interp, cmap='YlOrBr')
	plt.contour(Mf_vec, af_vec, P_Mfmean_afmean_interp, cmap='YlOrBr')
	plt.xlabel('$M_f$')
	plt.ylabel('$a_f$')
	plt.grid()
	plt.colorbar()

	plt.subplot(232)
	plt.pcolormesh(Mf_vec, af_vec, P_dMf_daf, cmap='YlOrBr')
	plt.contour(Mf_vec, af_vec, P_dMf_daf, cmap='YlOrBr')
	plt.xlabel('$\Delta M_f$')
	plt.ylabel('$\Delta a_f$')
	plt.grid()
	plt.colorbar()
	plt.title('Difference')

	plt.subplot(235)
	plt.pcolormesh(Mf_vec, af_vec, P_dMf_daf_interp, cmap='YlOrBr')
	plt.contour(Mf_vec, af_vec, P_dMf_daf_interp, cmap='YlOrBr')
	plt.xlabel('$M_f$')
	plt.ylabel('$a_f$')
	plt.grid()
	plt.colorbar()

	plt.subplot(233)
	plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMf_dafbyaf, cmap='YlOrBr')
	plt.contour(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMf_dafbyaf, cmap='YlOrBr')
	plt.xlabel('$\Delta M_f/\\bar{M_f}$')
	plt.ylabel('$\Delta a_f/\\bar{M_f}$')
	plt.grid()
	plt.colorbar()
	plt.title('$\epsilon, \sigma$')

	plt.subplot(236)
	plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMf_dafbyaf_interp, cmap='YlOrBr')
	plt.contour(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMf_dafbyaf_interp, cmap='YlOrBr')
	plt.xlabel('$\Delta M_f/\\bar{M_f}$')
	plt.ylabel('$\Delta a_f/\\bar{M_f}$')
	plt.grid()
	plt.colorbar()
	plt.tight_layout()
	plt.savefig('%s.png' %outfile, dpi=300)
	plt.close()

