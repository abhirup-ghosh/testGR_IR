""" 
Routines for performing overlap calculations. 

P. Ajith, 2013-02-05 
$Id: wavegen.py 232 2014-01-31 14:14:57Z ajith_p $
"""

import lalinspiral as lalinsp 
import lal, os 
import lalsimulation as lalsim
import numpy as np
from lal import MSUN_SI, MTSUN_SI, PC_SI, C_SI
from scipy import signal
#import spinm2SphHarm as sSH


""" taper time domain data. h is a numpy array (Ref. Eq. (3.35) of gr-qc/0001023) """
def taper_waveform(h):
	h_temp = h
	peakind = np.array(signal.argrelextrema(abs(h_temp), np.greater)).flatten()
	idx_peak2 = peakind[1]		# index of second extremum
	startind = np.flatnonzero(h_temp)[0]		# index of first non-zero data point
	
	# taper from start to second extremum 
	n = idx_peak2 - startind
	# do the taper using formula Eq. (3.35) of gr-qc/0001023.
	h_temp[startind] = 0
	for i in range(startind+1, startind+n-2):
		z = (n - 1.)/(i-startind) + (n-1.)/(i-startind - (n-1.))
		h_temp[i] = h_temp[i]*1./(np.exp(z) + 1)
	return h_temp


""" Generate h data from .txt files (currently only works for quadrupole mode) """
def generate_nrdata_text(incl_angle, phi, dataFileName):

	t_geom, hp_txt, hc_txt = np.loadtxt(dataFileName, usecols = (0,1,2), unpack=True, comments='#')  # load the text file 
	hp_22, hc_22 = hp_txt/0.6308, hc_txt/0.6308
	dt_geom = t_geom[1] - t_geom[0]

	h_22 = hp_22 - 1j*hc_22	
	y_22 = sSH.spinm2_sphharm(2,2,incl_angle,phi)
	y_22_minus = sSH.spinm2_sphharm(2,-2,incl_angle,phi)
	h_complex = h_22 * y_22 + (-1)**2 * h_22.conjugate() * y_22_minus
	
	hp_nr_geom = h_complex.real					# plus polarization data in geometric units, scaled by (M/D)
	hc_nr_geom = - h_complex.imag					# cross polarization data in geometric units, scaled by (M/D)

	return hp_nr_geom, hc_nr_geom, dt_geom	


""" Generate h data from .mat files """
def generate_nrdata(incl_angle, phi, dataFileName, waveform_type, use_only_22):

	Spinmode_hlm = sio.loadmat('%s'%dataFileName)  # load the mat file 
	def gwpolzonspere(theta,phi,lvec,mvec):

	   if use_only_22 == 0:
     	   	h = np.zeros(len( Spinmode_hlm[waveform_type]['h'][0][0][0][0][0]))		# initialize the h vector
		for i_mode in range(len(lvec)):							
			l,m = int(lvec[i_mode]),int(mvec[i_mode])			# the l,m values corresponding to i_mode		
			y_lm = sSH.spinm2_sphharm(l,m,theta,phi)	# spherical harmonic functions for mode l,m at theta,phi
			y_lm_minus = sSH.spinm2_sphharm(l,-m,theta,phi)  # spherical harmonic functions for mode l,-m at theta,phi
			h_lm = Spinmode_hlm[waveform_type]['h'][0][0][0][i_mode][0]		# h_lm for mode l,m
			if m == 0:					
				h = h + h_lm * y_lm
			elif m > 0:					# if m>0, the part corresponding to -m is given by the last term
				h = h + h_lm * y_lm + (-1)**l * h_lm.conjugate() * y_lm_minus

	   else:
		i_22 = np.intersect1d(np.where(mvec==2)[0], np.where(lvec==2)[0])[0]	# index corresponding to the l=2,m=2 mode
		y_22 = sSH.spinm2_sphharm(2,2,theta,phi)
		y_22_minus = sSH.spinm2_sphharm(2,-2,theta,phi)
		h_22 = Spinmode_hlm[waveform_type]['h'][0][0][0][i_22][0]
		h = h_22 * y_22 + (-1)**2 * h_22.conjugate() * y_22_minus		   

	   return h

	Lvec = Spinmode_hlm[waveform_type]['l'][0][0][0]		# vector of l values determined from the mat file
	Mvec = Spinmode_hlm[waveform_type]['m'][0][0][0]		# vector of m values determined from the mat file
	h_nr_complex_geom = gwpolzonspere(incl_angle , phi, Lvec, Mvec)		# complex nr data (hp-i*hc) in geometric units, scaled by (M/D)
	hp_nr_geom = h_nr_complex_geom.real					# plus polarization data in geometric units, scaled by (M/D)
	hc_nr_geom = - h_nr_complex_geom.imag					# cross polarization data in geometric units, scaled by (M/D)
	dt_geom = float(Spinmode_hlm[waveform_type]['dt'][0][0][0][0])		#time step in geometric units, scaled by (M/D)

	return hp_nr_geom, hc_nr_geom, dt_geom	


""" Generate waveforms from LALSimulation """
def generate_waveforms(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, dL, incl_angle, sky_theta, sky_phi, polarization_angl, f_low, f_final, fs, N, approx, amplO, phaseO, phi0=0., dataFileName=0, waveform_type=0, domain='freq', use_only_22=0, print_waveforms=0): 
	# for the case of non-precessing-spin approximants, set the x, y components 
	# of the spins to zero 
	if approx == 'SEOBNRv2' or approx == 'SEOBNRv1' or approx == 'IMRPhenomB'  or approx == 'IMRPhenomC' or approx == 'IMRPhenomFB' or approx == 'IMRPhenomFC' or approx == 'TaylorF2RedSpin' or approx == 'TaylorF2': 
		if spin1x > 5e-4 or  spin1y > 5e-4 or spin2x > 5e-4 or  spin2y > 5e-4:	# this is a hack
			print 'WARNING from %s: Non-zero transverse spins were given, but setting them to zero since the chosen approximant is non-precessing' %os.path.basename(__file__)
		if spin1x != 0 or  spin1y != 0 or spin2x != 0 or  spin2y != 0: 
			spin1x = spin1y = spin2x = spin2y = 0.
	elif approx == 'TaylorT1' or approx == 'EOBNRv2': 
		if spin1x > 5e-4 or  spin1y > 5e-4 or spin2x > 5e-4 or  spin2y > 5e-4 or spin1z > 5e-4 or spin2z > 5e-4:# this is a hack
			print 'WARNING from %s: Non-zero spins were given, but setting them to zero since the chosen approximant is non-spinning' %os.path.basename(__file__)
		if spin1x != 0 or  spin1y != 0 or spin1z != 0 or spin2x != 0 or  spin2y != 0 or spin2z != 0: 
			spin1x = spin1y = spin1z = spin2x = spin2y = spin2z = 0.

	# this is a hack -- to get an epoch variable - FIXME 
	dummy = lalsim.SimInspiralTaylorF2ReducedSpin(0, 1, 10 * MSUN_SI, 10 * MSUN_SI, 0, 40, 41, 1e6 * PC_SI, 0, 0)

	if approx != 'NRPNHybrid' and approx != 'NRPNHybridmat' and approx != 'NRPNHybridtxt': 			# df is defined for NRPRHybrid separately
		# frequency resolution 
		df = fs/N 
		# get the approximant enum for LALSimulation 
		approx_enum = lalsim.GetApproximantFromString(approx)		# approx_enum is not needed for NRPNHybrid

	# generate templates and compute the Fourier transform 
	if domain == 'time':

		# read time-domain hybrid waveforms from mat file 
		if approx ==  'NRPNHybrid' or approx == 'NRPNHybridmat' or approx == 'NRPNHybridtxt': 
	     	   if approx == 'NRPNHybrid' or approx == 'NRPNHybridmat':
			# read the Hybrid data from mat files
			hp_nr_geom, hc_nr_geom, dt_geom = generate_nrdata(incl_angle, phi0, dataFileName, waveform_type, use_only_22) 
		
		   if approx ==  'NRPNHybridtxt':
			# read the Hybrid data from txt files
			hp_nr_geom, hc_nr_geom, dt_geom = generate_nrdata_text(incl_angle, phi0, dataFileName)

	    	   # rescale to get physical waveform for m1 and m2
		   m = m1 + m2			
      	   	   hp_nr = hp_nr_geom * m * MTSUN_SI * C_SI/ (1.0e6 * PC_SI)
		   hc_nr = hc_nr_geom * m * MTSUN_SI *  C_SI / (1.0e6 * PC_SI)
		   dt = dt_geom * m * MTSUN_SI			
		   fs = 1./dt	 	# sample rate

		   # make sure hp_nr and hc_nr have same lengths
		   if len(hp_nr) != len(hc_nr):
			print 'Different lengths for the plus and cross waveform data (%3.2e and %3.2e)' %(len(hp_nr), len(hc_nr))
			exit(-1) 

		   # convert the NR data into REAL8TimeSeries
		   hp = lal.CreateREAL8TimeSeries("hp(t)", dummy.epoch, 0.0, dt, lal.StrainUnit, len(hp_nr))
		   hc = lal.CreateREAL8TimeSeries("hc(t)", dummy.epoch, 0.0, dt, lal.StrainUnit, len(hc_nr))
		   hp.data.data = hp_nr
		   hc.data.data = hc_nr

		   # project the polarizations to the detector: h = Fp hp + Fc hc 
		   h = project_hplus_hcross(hp, hc, sky_theta, sky_phi, polarization_angl)

		elif approx ==  'SpinTaylorT5':

			hp, hc = lalsim.SimInspiralSpinTaylorT5(
				0,				# GW phase at reference freq (rad)
				1./fs,				# sampling interval (s)
				m1 * lal.MSUN_SI,	# mass of companion 1 (kg)
				m2 * lal.MSUN_SI,	# mass of companion 2 (kg)
				f_low,			# start GW frequency (Hz)
				dL*1e6*PC_SI,			# distance of source (m)
				spin1x,				# initial value of S1x
				spin1y,				# initial value of S1y
				spin1z,				# initial value of S1z
				spin2x,				# initial value of S2x
				spin2y,				# initial value of S2y
				spin2z,				# initial value of S2z
				incl_angle,			# incl_angle angle - careful with definition (line of sight to total vs orbital angular momentum)
				phaseO,				# twice PN phase order
				amplO
				)
			
			# project the polarizations to the detector: h = Fp hp + Fc hc 
			h = project_hplus_hcross(hp, hc, sky_theta, sky_phi, polarization_angl)

		elif approx ==  'SpinTaylorT5_v2':
			import spa_complete_waveform as spa
			
			hp_t5, hc_t5, dt =  spa.spintaylorF2_generate(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, incl_angle, fs, N, approx)
			# make sure hp_t5 and hc_t5 have same lengths
			if len(hp_t5) != len(hc_t5):
				print 'Different lengths for the plus and cross waveform data (%3.2e and %3.2e)' %(len(hp_t5), len(hc_t5))
				exit(-1) 

			# convert the NR data into REAL8TimeSeries
			hp = lal.CreateREAL8TimeSeries("hp(t)", dummy.epoch, 0.0, dt, lal.StrainUnit, len(hp_t5))
			hc = lal.CreateREAL8TimeSeries("hc(t)", dummy.epoch, 0.0, dt, lal.StrainUnit, len(hc_t5))
			hp.data.data = hp_t5
			hc.data.data = hc_t5

			# project the polarizations to the detector: h = Fp hp + Fc hc 
			h = project_hplus_hcross(hp, hc, sky_theta, sky_phi, polarization_angl)

		# Generate template waveforms in the time domain from LALSimulation 
		else: 
			hp, hc = lalsim.SimInspiralChooseTDWaveform(
				0.,						# phi_ref (rad) 
				1./fs,					# delta_t (sec) 
				m1 * MSUN_SI, 		# mass 1 (kg) 
				m2 * MSUN_SI,		# mass 2 (kg) 
				spin1x,					# spin1x (S1x/m1^2)  
				spin1y, 				# spin1y (S1y/m1^2)
				spin1z, 				# spin1z (S1z/m1^2)
				spin2x, 				# spin2x (S2x/m2^2)
				spin2y, 				# spin2y (S2y/m2^2)
				spin2z,					# spin2z (S2z/m2^2)
				f_low, 					# f_low (Hz) 
				f_low, 					# f_ref (Hz) 
				dL*1e6*PC_SI, 		# distance (m)
				incl_angle, 			# inclination angle (rad) 
				0.,						# lambda1 (Love number) 
				0.,						# lambda2 (Love number) 
				None, 					# waveform flags
				None, 					# non-GR parameters
				amplO, 					# amplitude PN order 
				phaseO, 				# phase PN order 
				approx_enum 			# approximant 
				)
	

			# project the polarizations to the detector: h = Fp hp + Fc hc 
			h = project_hplus_hcross(hp, hc, sky_theta, sky_phi, polarization_angl)
			

		# taper
		lalsim.SimInspiralREAL8WaveTaper(h.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)

		if print_waveforms:
			t = np.linspace(0, h.data.length/fs, h.data.length) 
			np.savetxt('%s_%2.1fPN_m1_%2.1f_m2_%2.1f_flow_%2.1f.dat.gz' %(approx, phaseO/2., m1, m2, f_low), np.column_stack((t, h.data.data)))

		# if the length of the generated data is larger than the prescribed N, set N to 
		# be the length og the data (after 'ceiling' to the nearest power of 2)
		if int(2**np.ceil(np.log2(h.data.length))) > N:
			N = int(2**np.ceil(np.log2(h.data.length)))

		# zero-pad to the required length 
		h = lal.ResizeREAL8TimeSeries(h, 0, N)
		df = fs/N
			
		# create vector to hold the FFT and FFT plan, do the FFT
		hf = lal.CreateCOMPLEX16FrequencySeries("h(f)", h.epoch, h.f0, df, lal.HertzUnit, int(N/2+1))
		fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)
		lal.REAL8TimeFreqFFT(hf, h, fftplan)

	elif domain == 'freq':

  	   # frequency resolution 
	   df = fs/N

           if approx == 'SpinTaylorF2':
		import spa_complete_waveform as spa
		hpf_spa, hcf_spa =  spa.spintaylorF2_generate(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, incl_angle, fs, int(N/2+1), approx)

		# create vector to hold the FFT and FFT plan, do the FFT
		hpf = lal.CreateCOMPLEX16FrequencySeries("h(f)", dummy.epoch, 0., df, lal.HertzUnit, int(N/2+1))   #hf.f0 set to zero
		hcf = lal.CreateCOMPLEX16FrequencySeries("h(f)", dummy.epoch, 0., df, lal.HertzUnit, int(N/2+1))   #hf.f0 set to zero
		hpf.data.data = hpf_spa
		hcf.data.data = hcf_spa

	   else:

		# Generate template waveforms in the Fourier domain
		hpf, hcf = lalsim.SimInspiralChooseFDWaveform(
			0,						# phi_ref (rad) 
			df,						# delta_f (Hz) 
			m1 * MSUN_SI, 		# mass 1 (kg) 
			m2 * MSUN_SI,		# mass 2 (kg) 
			spin1x,					# spin1x (S1x/m1^2)  
			spin1y, 				# spin1y (S1y/m1^2)
			spin1z, 				# spin1z (S1z/m1^2)
			spin2x, 				# spin2x (S2x/m2^2)
			spin2y, 				# spin2y (S2y/m2^2)
			spin2z,					# spin2z (S2z/m2^2)
			f_low, 					# f_low (Hz) 
			f_final, 				# f_max (Hz) 
			f_low, 					# f_ref (Hz) 
			dL*1e6*PC_SI, 		# distance (m)
			incl_angle, 			# inclination angle (rad) 
			0,						# lambda1 (Love number) 
			0,						# lambda2 (Love number) 
			None, 					# waveform flags
			None, 					# non-GR parameters
			amplO, 					# amplitude PN order 
			phaseO, 				# phase PN order 
			approx_enum 			# approximant 
			)

	   # project the polarizations to the detector: h = Fp hp + Fc hc 
     	   hf = project_hplus_hcross(hpf, hcf, sky_theta, sky_phi, polarization_angl)

	else:
		print '### %s : Unknown approximant %s' %(os.path.basename(__file__), approx)
		exit(-1)
		
	return hf


def project_hplus_hcross(hplus, hcross, theta, phi, psi):

	# compute antenna factors Fplus and Fcross
	Fp = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.cos(2*psi) - np.cos(theta)*np.sin(2*phi)*np.sin(2*psi)
	Fc = 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.sin(2*psi) + np.cos(theta)*np.sin(2*phi)*np.cos(2*psi)

	# allocate memory for the resulting vector 
	if type(hplus) == lal.COMPLEX16FrequencySeries: 
		h = lal.CreateCOMPLEX16FrequencySeries("h(f)", hplus.epoch, hplus.f0, hplus.deltaF, lal.HertzUnit, hplus.data.length)
	elif type(hplus) == lal.REAL8TimeSeries: 
		h = lal.CreateREAL8TimeSeries("h(t)", hplus.epoch, hplus.f0, hplus.deltaT, lal.SecondUnit, hplus.data.length)

	# form strain signal in detector
	h.data.data = Fp*hplus.data.data + Fc*hcross.data.data

	return h
