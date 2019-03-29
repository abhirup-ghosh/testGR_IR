import numpy as np
import noiseutils as nu
import matplotlib.pyplot as plt
import lal
import os
import imrtestgr as tgr
import nr_fits as nr
import wavegen as wf

# parameters of the waveform
m1 = 3.
m2 = 3.
spin1x, spin1y, spin1z = 0., 0., 0.
spin2x, spin2y, spin2z = 0., 0., 0.
incl_angle, sky_theta, sky_phi, polarization_angl = 0., 0., 0., 0.,

approx = 'IMRPhenomPv2'
fs = 16384.
N = fs*8.
f_min = 10.
f_max = fs/2.

approx = 'IMRPhenomB' #'EOBNRv2'
dom='freq' # 'time'
amplO, phaseO = 8, 8
fit_formula = 'nonprecspin_Healy2014'

# generate the Fourier transform of h+(f) 
hf = wf.generate_waveforms(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, fs, N, approx, amplO, phaseO, phi0=0., dataFileName=0, waveform_type=0, domain=dom, use_only_22=0, print_waveforms=0)

# unpack the hf 
hf = hf.data.data
hf = abs(hf)

# construct a freq vector
f = np.linspace(0., f_max, len(hf))

# calculating the SNR post-cutoff
f_cut_list = np.linspace(f_min, f_max, 500)
snr_r_list = np.zeros(len(f_cut_list))

for (i, f_cut) in enumerate(f_cut_list):
    snr_r_list[i] = nu.SNR_from_fd(f, hf, f_low=f_cut, f_high=f_max)

# scaling by the full IMR SNR (to see what fraction lies over cutoff)
snr_imr = nu.SNR_from_fd(f, hf, f_low=f_min, f_high=f_max)
snr_r_list = snr_r_list/snr_imr

# calculating the Mf, af, Kerr ISCO freq, dominant QNM freq 
Mf, af = tgr.calc_final_mass_spin(m1, m2, spin1z, spin2z, fit_formula)
f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)
f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)
print 'Mf = %3.2f af = %3.2f f_isco_Kerr = %2.1f Hz f_qnm = %2.1f Hz ' %(Mf, af, f_isco_Kerr, f_qnm)

# plotting
plt.figure()
plt.plot(f_cut_list, snr_r_list)
plt.axvline(x=f_isco_Kerr, color='g', ls='--', label='$f_{Kerr ISCO}$')
plt.axvline(x=f_qnm, color='r', ls='--', label='$f_{QNM}$')
plt.xlabel('flow for ringdown (Hz)')
plt.ylabel('$SNR_r/ SNR_{IMR}$')
plt.title('(%.1f + %.1f) $M_{\odot}$ system'%(m1, m2))
plt.legend(loc='best')
plt.grid()
plt.xlim([0, 2048])
plt.savefig('snr_falloff_%s_m1_%.2f_m2_%.2f_s1z_%.2f_s2z_%.2f_fmin_%.2fHz.png'%(approx, m1, m2, spin1z, spin2z, f_min))
