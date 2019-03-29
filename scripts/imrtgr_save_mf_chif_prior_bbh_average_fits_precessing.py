import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import nrutils as nr
from scipy import interpolate as interp
import pickle, gzip

N_samples = int(1e7)
Mf_lim = 500
af_lim = 1
N_bins = 400
comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max = 5., 170., 0., 0.99


m1_pr = np.random.uniform(comp_mass_min, comp_mass_max, N_samples)
m2_pr = np.random.uniform(comp_mass_min, comp_mass_max, N_samples)

chi1_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_samples)
chi2_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_samples)
tilt1_pr = np.arccos(np.random.uniform(-1, 1, N_samples))
tilt2_pr = np.arccos(np.random.uniform(-1, 1, N_samples))
phi12_pr = np.random.uniform(-1, 1, N_samples)*np.pi

Mf_fits = ["UIB2016", "HL2016"]
af_fits = ["UIB2016", "HBR2016", "HL2016"]

Mf_pr = nr.bbh_average_fits_precessing(m1_pr, m2_pr, chi1_pr, chi2_pr, tilt1_pr, tilt2_pr, 0, "Mf", Mf_fits)
af_pr = nr.bbh_average_fits_precessing(m1_pr, m2_pr, chi1_pr, chi2_pr, tilt1_pr, tilt2_pr, phi12_pr, "af", af_fits)

Mf_bins = np.linspace(-Mf_lim, Mf_lim, N_bins)
af_bins = np.linspace(-af_lim, af_lim, N_bins)

P_Mfaf_pr, Mf_bins, af_bins = np.histogram2d(Mf_pr, af_pr, bins=(Mf_bins, af_bins), normed=True)
P_Mfaf_pr = P_Mfaf_pr.T
print '... calculated the prior'

# create an interpolation object and save it 
Mf_bins = (Mf_bins[:-1] + Mf_bins[1:])/2.
af_bins = (af_bins[:-1] + af_bins[1:])/2.
outfile = 'Prior_Mfaf_bbh_average_fits_precessing_comp_mass_min%2.1f_comp_mass_max%2.1f_comp_spin_min%2.1f_comp_spin_max%2.2f'%(comp_mass_min, comp_mass_max, comp_spin_min, comp_spin_max)
P_Mfaf_pr_interp_obj = interp.interp2d(Mf_bins, af_bins, P_Mfaf_pr, fill_value=0., bounds_error=False)
f = gzip.open(outfile+".pklz",'wb')
pickle.dump(P_Mfaf_pr_interp_obj, f)
print '... saved the interpolation object.'

# read the interpolation object, reconstruct the data from the interpolation object. This is only used for estimating the error due to the interpolation 
f = gzip.open(outfile+".pklz",'rb')
P_Mfaf_pr_interp_obj = pickle.load(f)
P_Mfaf_pr_interp = P_Mfaf_pr_interp_obj(Mf_bins, af_bins)

# difference between the original and interpolated data 
interp_err = abs(P_Mfaf_pr-P_Mfaf_pr_interp)
print '... maximum difference between the original and interpolated data is %e' %np.max(interp_err)

plt.figure(figsize=(15,4))
plt.subplot(131)
plt.pcolormesh(Mf_bins, af_bins, P_Mfaf_pr, cmap='YlOrBr')
plt.contour(Mf_bins, af_bins, P_Mfaf_pr, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.colorbar()
plt.title('Original')

plt.subplot(132)
plt.pcolormesh(Mf_bins, af_bins, P_Mfaf_pr_interp, cmap='YlOrBr')
plt.contour(Mf_bins, af_bins, P_Mfaf_pr_interp, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.colorbar()
plt.title('Interpolation')

plt.subplot(133)
plt.pcolormesh(Mf_bins, af_bins, interp_err, cmap='YlOrBr')
plt.xlabel('$M_f [M_{\odot}]$')
plt.ylabel('$a_f/M_f$')
plt.xlim(-Mf_lim, Mf_lim)
plt.ylim(-af_lim, af_lim)
plt.grid()
plt.clim(np.min(interp_err)-1e-16,np.max(interp_err)+1e-16)
plt.colorbar()
plt.title('Difference')
plt.savefig(outfile+'.png', dpi=300)

# save a zoomed version of the figure 
plt.subplot(131)
plt.xlim(50,150)
plt.ylim(0,1)
plt.subplot(132)
plt.xlim(50,150)
plt.ylim(0,1)
plt.subplot(133)
plt.xlim(50,150)
plt.ylim(0,1)
plt.savefig(outfile+'_zoom.png', dpi=300)
