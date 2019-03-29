# Code to read in LALInference posterior samples and calculate the medians (and other quantities, including histograms) of the posteriors on the ISCO frequency of the final Kerr black hole and on the peak frequency of the 2,2 mode, using the Taracchini et al. fits for the last of these. Note that we work with the redshifted masses throughout.
# NKJ-M, 12.2015, based on earlier code for the radiated mass and peak luminosity

import numpy as np
import bayesian as ba
import confidence as conf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nr_fits as nrf
import tgrplotsettings 
import argparse

# Set number of bins

N_bin = 100

# Conversion factor Msun -> s

Msun2s = 4.92549e-6 # s/Msun

# Set up the parsing

parser = argparse.ArgumentParser(description = 'Calculate the posteriors on the ISCO frequency of the final black hole and the peak frequency of the 2,2 mode from LALInference posterior samples.')
parser.add_argument("datfile", help = "posterior_samples.dat file to read in")
parser.add_argument("tag", help = "tag for output files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
#parser.add_argument("--redshifted_m", help = "use the redshifted masses instead of the source frame masses (if only the redshifted masses are present, the code uses these and prints a warning)", action = "store_true")
#group = parser.add_mutually_exclusive_group()
#group.add_argument("--Erad", help = "only calculate the radiated energy", action = "store_true")
#group.add_argument("--Lpeak", help = "only calculate the peak luminosity", action = "store_true")
args = parser.parse_args()

print(args.tag)
print("%s\n"%(args.outdir))

# read the posterior samples for the masses and spins, taking the component along the angular momentum for precessing spins
data = np.genfromtxt(args.datfile, dtype=None, names=True)
#if args.redshifted_m:
m1, m2 = data['m1'], data['m2']
#else:
#  if ('m1_source' in data.dtype.names) and ('m2_source' in data.dtype.names):
#    m1, m2 = data['m1_source'], data['m2_source']
#  else:
#    print("WARNING! Source frame masses not available, using redshifted masses instead.")
#    m1, m2 = data['m1'], data['m2']
if ('a1z' in data.dtype.names) and ('a2z' in data.dtype.names):
  chi1, chi2 = data['a1z'], data['a2z']
else:
  chi1, chi2 = np.zeros(len(m1)), np.zeros(len(m2))
#if ('costilt1' in data.dtype.names) and ('costilt2' in data.dtype.names):
#  chi1, chi2 = chi1*data['costilt1'], chi2*data['costilt2']

# Compute the quantities for the ISCO frequency (using the Healy et al. fit to find the final mass and spin)

Mf, af = nrf.bbh_final_mass_and_spin_non_precessing(m1, m2, chi1, chi2)

fISCO = nrf.calc_isco_freq(af)/(Mf*Msun2s)

P_fISCO, fISCO_bins = np.histogram(fISCO, bins=N_bin, normed=True)
P_fISCO = P_fISCO.T

fISCO_median = np.median(fISCO)

fISCO_maP = ba.edg2cen(fISCO_bins)[np.nanargmax(P_fISCO)]

fISCO_conf = conf.confidence_1d(P_fISCO, fISCO_bins)

fISCO_68min, fISCO_68max = fISCO_conf.edges(0.68)
fISCO_90min, fISCO_90max = fISCO_conf.edges(0.90)
fISCO_95min, fISCO_95max = fISCO_conf.edges(0.95)

print("fISCO values (all in Hz):")
print("median: %.3f"%(fISCO_median))
print("maP: %.3f"%(fISCO_maP))
print("0.68 region: %.3f, %.3f"%(fISCO_68min, fISCO_68max))
print("0.90 region: %.3f, %.3f"%(fISCO_90min, fISCO_90max))
print("0.95 region: %.3f, %.3f\n"%(fISCO_95min, fISCO_95max))

# Compute the peak frequency

fpeak = nrf.bbh_22_peak_frequency_non_precessing_Taracchinietal(m1, m2, chi1, chi2)

P_fpeak, fpeak_bins = np.histogram(fpeak, bins=N_bin, normed=True)
P_fpeak = P_fpeak.T

fpeak_median = np.median(fpeak)

fpeak_maP = ba.edg2cen(fpeak_bins)[np.nanargmax(P_fpeak)]

fpeak_conf = conf.confidence_1d(P_fpeak, fpeak_bins)

fpeak_68min, fpeak_68max = fpeak_conf.edges(0.68)
fpeak_90min, fpeak_90max = fpeak_conf.edges(0.90)
fpeak_95min, fpeak_95max = fpeak_conf.edges(0.95)

print("fpeak values (all in Hz):")
print("median: %.3f"%(fpeak_median))
print("maP: %.3f"%(fpeak_maP))
print("0.68 region: %.3f, %.3f"%(fpeak_68min, fpeak_68max))
print("0.90 region: %.3f, %.3f"%(fpeak_90min, fpeak_90max))
print("0.95 region: %.3f, %.3f\n"%(fpeak_95min, fpeak_95max))

# Make plots

plt.figure(figsize=(16,8))
plt.subplot2grid((2, 4), (0, 0))
plt.hist(m1, bins=N_bin, normed=True)
plt.xlabel('$m_1~[M_\odot]$')
plt.ylabel('$P(m_1)$')

plt.subplot2grid((2, 4), (0, 1))
plt.hist(m2, bins=N_bin, normed=True)
plt.xlabel('$m_2~[M_\odot]$')
plt.ylabel('$P(m_2)$')

plt.subplot2grid((2, 4), (0, 2))
plt.hist(chi1, bins=N_bin, normed=True)
plt.xlabel('$\chi_1$')
plt.ylabel('$P(\chi_1)$')

plt.subplot2grid((2, 4), (0, 3))
plt.hist(chi2, bins=N_bin, normed=True)
plt.xlabel('$\chi_2$')
plt.ylabel('$P(\chi_2)$')

plt.subplot2grid((2, 4), (1, 0), colspan = 2)
plt.bar(fISCO_bins[:-1], P_fISCO, width=fISCO_bins[1] - fISCO_bins[0])
plt.xlim(min(fISCO_bins), max(fISCO_bins))
plt.xlabel('$f_\mathrm{ISCO}$ [Hz]')
plt.ylabel('$P(f_\mathrm{ISCO})$')

plt.subplot2grid((2, 4), (1, 2), colspan = 2)
plt.bar(fpeak_bins[:-1], P_fpeak, width=fpeak_bins[1] - fpeak_bins[0])
plt.xlim(min(fpeak_bins), max(fpeak_bins))
plt.xlabel('$f_\mathrm{peak}$ [Hz]')
plt.ylabel('$P(f_\mathrm{peak})$')

plt.tight_layout()
plt.savefig('%s/fISCO_fpeak_debugplots_%s.png' %(args.outdir,args.tag), dpi=200)
