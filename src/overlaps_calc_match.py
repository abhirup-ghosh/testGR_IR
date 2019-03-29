""" 
Routines for performing overlap calculations. 

P. Ajith, 2013-02-05 
$Id: overlaps.py 234 2014-02-07 11:38:52Z vijayvarma $
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P
import lalinspiral as lalinsp 
import lal, os 
import lalsimulation as lalsim
import lalinspiral.sbank.waveforms as wf
import numpy as np
import wavegen as wg
from scipy.optimize import minimize
import lalinspiral.sbank.psds as psds
import plotsettings 
from scipy.interpolate import interp1d

""" Generate Advanced LIGO High-Power Zero Detuning PSD using the fit given in Ajith (2011) """
def adligo_psd(f, fLow):
	x = f/245.4;
	Sh = 1e-48*(0.0152*x**-4 + 0.2935*x**(9./4.) + 2.7951*x**(3./2.) - 6.5080*x**(3./4.) + 17.7622)
	df = f[1]-f[2]
	#Sh[:int(fLow / df)] = 0.	# for f < fLow, PSD = 0
	return Sh 

""" Whiten and normalize a frequency-domain waveform """
def whiten_and_normalize_waveform(hf, Sh, f_low, f_final, df):

	# whiten the waveforms 
	hf.data.data[:] /= Sh[:hf.data.length]**0.5 

	# limit the frquency band of integration of f_low to f_final 
	hf.data.data[:int(f_low / df)] = 0.
	hf.data.data[int(f_final / df) : hf.data.length] = 0.

	# normalize the waveforms 
	sigmasq = wf.compute_sigmasq(hf.data.data, df)
	hf.data.data[:] /= sigmasq**0.5

	return hf, sigmasq 


""" Compute the match of a frequency domain waveform h with a frequency-domain template x """
#NOTE that the target waveform htilde is assumed to be whitened, limited to 
# [f_low, f_final], normalized and converted to single precision
def compute_match (htilde, xtilde, Sh, f_low, f_final, df, plot_fft=0,tag='no_tag'): 

	# whiten and normalize the template 
	xtilde, sigmasq_x = whiten_and_normalize_waveform(xtilde, Sh, f_low, f_final, df)

	# down-convert to singla precision
	xtilde = wf.FrequencySeries_to_COMPLEX8FrequencySeries(xtilde)

	# compute matches 
	workspace_cache = lalinsp.CreateSBankWorkspaceCache()
	match = lalinsp.InspiralSBankComputeMatch(htilde, xtilde, workspace_cache)
	if plot_fft:
		fh = np.linspace(0, df*htilde.data.length, htilde.data.length)
		fx = np.linspace(0, df*xtilde.data.length, xtilde.data.length)
		phi_h = np.unwrap(np.angle(htilde.data.data))
		phi_x = np.unwrap(np.angle(xtilde.data.data))
		phi_h-phi_h[0]
		phi_x-phi_x[0]

		P.figure(figsize=(18,5))
		P.subplot(131) 
		P.loglog(fh, abs(htilde.data.data), lw=2, alpha=1, color='red', label='target')
		P.loglog(fx, abs(xtilde.data.data), lw=1, alpha=1, color='black', label='template')
		P.ylabel('$A(f)$')
		P.xlabel('$f$ [Hz]')
		P.legend(loc=3)
		P.xlim(f_low, f_final)

		P.subplot(132) 
		P.semilogx(fh, phi_h, lw=2, alpha=1, color='red', label='target')
		P.semilogx(fh, phi_x, lw=1, alpha=1, color='k', label='template')
		P.ylabel('$\Psi(f)$')
		P.xlabel('$f$ [Hz]')
		P.xlim(f_low, f_final)
		P.title('match = %5.4f, norm = %5.4e' %(match, sigmasq_x**0.5))

		P.subplot(133) 
		P.semilogx(fh, phi_h-phi_x, lw=2, alpha=1, color='cyan', label='diff')
		P.ylabel('$\Delta \Psi(f)$')
		P.xlabel('$f$ [Hz]')
		P.xlim(f_low, f_final)

		P.savefig('fft_whitened_%s.png'%tag)
		P.close()

	# del workspace_cache
	return match, sigmasq_x
 
""" Find the individual masses of the binary from the total mass and \eta """
def individual_masses(m,eta):
   
   # the individual masses of the binary
   m1 = (m/2.)*(1.+ np.sqrt(1.-4.*eta)) 
   m2 = (m/2.)*(1.- np.sqrt(1.-4.*eta))

   return m1, m2

###################################################### Fitting Factor for Hybrid templates ####################################################################################################

""" Find the maximum match between NR target waveform and analytical templates by a 3d maximization """
def calcff_mchirp_eta_chi(m_targ, eta_targ, spin1z_targ, spin2z_targ, incl_angle,\
	 phi0, dataFileName, sample_rate, mchirp_start, eta_start, chi_start, m_templ_min, m_templ_max, eta_templ_min, eta_templ_max, chi_templ_min,\
	 chi_templ_max, maximization_method, f_min, f_max, f_min_wave, f_max_wave, psd,  approx_targ, approx_templ, amplO_templ, amplO_targ,\
	 phaseO_templ, phaseO_targ, spin1x, spin1y, spin2x, spin2y, sky_theta, sky_phi, polarization_angl, waveform_type, domain_targ, domain_templ,N,\
	 only_quadrupole, only_overlap, print_steps, fp, maximization_parameters, chi_definition):

	################################ generate the target waveform and the psd #############################################################
		
	# the individual masses of the target binary
	m1_targ, m2_targ = individual_masses(m_targ,eta_targ)

	# target waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain   
	hf = wg.generate_waveforms(m1_targ, m2_targ, spin1x, spin1y, spin1z_targ, spin2x, spin2y, spin2z_targ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sample_rate, N, approx_targ, amplO_targ, phaseO_targ, phi0, dataFileName,waveform_type, domain=domain_targ, use_only_22=only_quadrupole)

	# generate psd
	f = np.linspace(0, sample_rate/2., hf.data.length)   		# frequency axis
	if psd == 'AdvLIGO-ZeroDet-HighP':	# use the fit to AdvLigo Zero detuned high power given in Ajith (2011)
		Sh = adligo_psd(f, f_min)					# PSD
	elif psd == 'WhiteNoise': 
		Sh = np.ones(len(f))
	else: 
		Sh = psds.noise_models[psd](f)
	
	# whiten and normalize the target 
	hf, sigmasq_h = whiten_and_normalize_waveform(hf, Sh, f_min, f_max, hf.deltaF)
	
	# down-convert to singla precision
	hf = wf.FrequencySeries_to_COMPLEX8FrequencySeries(hf)

	if maximization_parameters == 'mchirp_eta_chi':
		first_param_start = mchirp_start
	elif maximization_parameters == 'm_eta_chi':
		m_start = mchirp_start/eta_start**(3./5)	
		first_param_start = m_start

################# Optimization over mchirp and \eta  if only one value of chi_templ is given ###############################################################
	
	## maximize the match over mchirp and \eta using a 2-d maximization algorithm #############################################################  
	if chi_templ_min == chi_templ_max:			
		only_2d_maximization = 1 			# only 2d maximization
		# function to be maximized. Note that maximizing fun_match is same as minimizing -(fun_match).
		fun_match = lambda x: -(optimize_mchirp_eta_chi(x[0], x[1], chi_templ_min,  hf, Sh, spin1x, spin1y, spin2x, spin2y, m_templ_min,\
		 m_templ_max, eta_templ_min, eta_templ_max, chi_templ_min, chi_templ_max, sky_theta, sky_phi, polarization_angl, incl_angle, f_min_wave,\
		 f_max_wave, f_min, f_max, sample_rate, N, approx_templ, domain_templ, amplO_templ, phaseO_templ, only_2d_maximization, print_steps, fp, maximization_parameters, chi_definition))

		# maximize the match
		res = minimize(fun_match, (first_param_start,eta_start), method=maximization_method) # , options={'xtol': 1e-6, 'ftol': 1e-6})

		chi_mm = chi_templ_min	  				# chi corresponding to maximum match (mm). In this case chi is kept fixed.
		
	################## Optimization over mchirp, \eta and chi if the lower and upper limits of chi_templ is given ######################################################

	## maximize the match over mchirp, \eta and \chi using a 3-d maximization algorithm #############################################################  
	else:
		only_2d_maximization = 0			# 3d maximization
		# function to be maximized. Note that maximizing fun_match is same as minimizing -(fun_match).
		fun_match = lambda x: -(optimize_mchirp_eta_chi(x[0], x[1], x[2],  hf, Sh, spin1x, spin1y, spin2x, spin2y, m_templ_min, m_templ_max,\
		 eta_templ_min, eta_templ_max, chi_templ_min, chi_templ_max, sky_theta, sky_phi, polarization_angl, incl_angle, f_min_wave, f_max_wave, f_min,\
		 f_max, sample_rate, N, approx_templ, domain_templ, amplO_templ, phaseO_templ, only_2d_maximization, print_steps, fp, maximization_parameters, chi_definition))

		# maximize the match
		res = minimize(fun_match, (first_param_start,eta_start,chi_start), method=maximization_method)#, options={'xtol': 1e-14, 'disp': True,'maxfev':500,'ftol': 1e-14})
		chi_mm = float(res.x[2])	  			# chi corresponding to maximum match (mm)

	## get the max match value
	max_match= -res.fun 			# The fitting factor, note the minus sign
	eta_mm = float(res.x[1]) 					# \eta corresponding to maximum match (mm)
	
	if maximization_parameters == 'mchirp_eta_chi':
		mchirp_mm = float(res.x[0])					# Chirp mass corresponding to maximum match (mm)
		m_mm= mchirp_mm/eta_mm**(3./5)  				# Total mass corresponding to maximum match (mm)
	elif maximization_parameters == 'm_eta_chi':
		m_mm = float(res.x[0])					# Total mass corresponding to maximum match (mm)
		mchirp_mm= m_mm * eta_mm**(3./5)  			# Chirp mass corresponding to maximum match (mm)

	if print_steps == 1:
		fp.close()

	return mchirp_mm, m_mm, eta_mm, chi_mm, max_match


""" Find the maximum match between NR target waveform and analytical templates over a grid of template parameters near the Nelder-Mead best match parameters """
# NOTE only suitable for 2-D maximization at the moment
def calcff_grid(m_targ, eta_targ, m_best_nelder, eta_best_nelder, spin1z_targ, spin2z_targ, incl_angle,\
	 phi0, dataFileName, sample_rate, m_templ_min, m_templ_max, eta_templ_min, eta_templ_max, chi_templ_min,\
	 chi_templ_max, maximization_method, f_min, f_max, f_min_wave, f_max_wave, psd,  approx_targ, approx_templ, amplO_templ, amplO_targ,\
	 phaseO_templ, phaseO_targ, spin1x, spin1y, spin2x, spin2y, sky_theta, sky_phi, polarization_angl, waveform_type, domain_targ, domain_templ,N,\
	 only_quadrupole, print_steps, fp, maximization_parameters, chi_definition):

	################################ generate the target waveform and the psd #############################################################
	
	# the individual masses of the target binary
	m1_targ, m2_targ = individual_masses(m_targ,eta_targ)

	# target waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain   
	hf = wg.generate_waveforms(m1_targ, m2_targ, spin1x, spin1y, spin1z_targ, spin2x, spin2y, spin2z_targ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sample_rate, N, approx_targ, amplO_targ, phaseO_targ, phi0, dataFileName,waveform_type, domain=domain_targ, use_only_22=only_quadrupole)

	# generate psd
	f = np.linspace(0, sample_rate/2., hf.data.length)   		# frequency axis
	if psd == 'AdvLIGO-ZeroDet-HighP':	# use the fit to AdvLigo Zero detuned high power given in Ajith (2011)
		Sh = adligo_psd(f, f_min)					# PSD
	elif psd == 'WhiteNoise': 
		Sh = np.ones(len(f))
	else: 
		Sh = psds.noise_models[psd](f)

	# whiten and normalize the target 
	hf, sigmasq_h = whiten_and_normalize_waveform(hf, Sh, f_min, f_max, hf.deltaF)
	
	# down-convert to singla precision
	hf = wf.FrequencySeries_to_COMPLEX8FrequencySeries(hf)

	############################# Find the overlap using the same parameters as the target #################################################################

################# Maximization over mchirp and \eta  if only one value of chi_templ is given ###############################################################
	
	## maximize the match over mchirp and \eta using a 2-d maximization #############################################################  
	if chi_templ_min == chi_templ_max:
		gridwidth_m = m_best_nelder - m_targ
		gridwidth_eta = eta_best_nelder - eta_targ
		m_grid_list = np.linspace(m_targ - gridwidth_m*0.1, m_targ + 1.1 * gridwidth_m, 10)
		eta_grid_list = np.linspace(eta_targ - gridwidth_eta*0.1, eta_targ + 1.1 * gridwidth_eta, 10)

		only_2d_maximization = 0 			# only 2d maximization
		# function to be maximized. Note that maximizing fun_match is same as minimizing -(fun_match).
		match_best_grid = 0
		m_best_grid = 0
		eta_best_grid = 0
		chi_best_grid = chi_templ_min
		for m_templ in m_grid_list:
			for eta_templ in eta_grid_list:
				if maximization_parameters == 'mchirp_eta_chi':
					first_param = m_templ*eta_templ**(3./5)
				elif maximization_parameters == 'm_eta_chi':
					first_param = m_templ
				## match 
				match = optimize_mchirp_eta_chi(first_param, eta_templ, chi_templ_min,  hf, Sh, spin1x, spin1y, spin2x, spin2y, m_templ_min,\
		 m_templ_max, eta_templ_min, eta_templ_max, chi_templ_min, chi_templ_max, sky_theta, sky_phi, polarization_angl, incl_angle, f_min_wave,\
		 f_max_wave, f_min, f_max, sample_rate, N, approx_templ, domain_templ, amplO_templ, phaseO_templ, only_2d_maximization, print_steps, fp, maximization_parameters, chi_definition)
				if match > match_best_grid:
					match_best_grid = match
					m_best_grid = m_templ
					eta_best_grid = m_templ
		mchirp_best_grid = m_best_grid * eta_best_grid**(3./5)
		
	if print_steps == 1:
		fp.close()

	return mchirp_best_grid, m_best_grid, eta_best_grid, chi_best_grid, match_best_grid


""" Generate templates for a given mass,\eta and \chi and find their match with target waveforms """
def optimize_mchirp_eta_chi(first_param, eta, chi,  hf, Sh, spin1x, spin1y, spin2x, spin2y, m_templ_min, m_templ_max, eta_templ_min, eta_templ_max,\
	 chi_templ_min, chi_templ_max, sky_theta, sky_phi, polarization_angl, incl_angle, f_min_wave, f_max_wave, f_min, f_max, sample_rate, N,\
	 approx_templ, domain_templ, amplO_templ, phaseO_templ, only_2d_maximization, print_steps, fp, maximization_parameters, chi_definition):

	if maximization_parameters == 'mchirp_eta_chi':
		mchirp = first_param
		m=mchirp/eta**(3./5)							# The total mass of the template waveform
	elif maximization_parameters == 'm_eta_chi':
		m = first_param

	## return 0 if the maximization goes beyond the bounds
	if only_2d_maximization == 1:				# 1 for 2d maximization 
		# apply some bounds on m and eta 
		if m > m_templ_max or m < m_templ_min or eta > eta_templ_max or eta < eta_templ_min: 
			return 0

	elif only_2d_maximization == 0:				# 0 for 3d maximization
		# apply some bounds on m, eta and chi
		if m>m_templ_max or m<m_templ_min or eta > eta_templ_max or eta < eta_templ_min or chi > chi_templ_max or chi < chi_templ_min: 
			return 0

	# generate templates for the mass m, symmetric mass ratio \eta and reduced spin, \chi and find their match with target waveforms 
	# the individual masses of the binary
	m1,m2 = individual_masses(m,eta)

	# Spins along z-axis for template waveform, note that these depend
	if chi_definition == 'IMRPhenomB_spin_parameter':
		spin1z_templ = spin2z_templ = chi
	elif chi_definition == 'reduced_spin':
		spin1z_templ = spin2z_templ = chi/(1.-76.*eta/113.)

	# template waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain 
	xf = wg.generate_waveforms(m1, m2, spin1x, spin1y, spin1z_templ, spin2x, spin2y, spin2z_templ, incl_angle, sky_theta, sky_phi, polarization_angl, \
		f_min_wave, f_max_wave, sample_rate, N, approx_templ, amplO_templ, phaseO_templ, domain=domain_templ)

	# make sure that the two vectors have the same length, freq resolution etc 
	if abs(xf.deltaF- hf.deltaF) < 1e-14:   #FIXME this is a hack 
		df = xf.deltaF
	else: 
		print 'Different frequency sampling for the target and template waveforms (%3.2e and %3.2e)' %(xf.deltaF, hf.deltaF)
		exit(-1) 
	if xf.data.length != hf.data.length:
		print 'Different lengths for the target and template waveform vectors (%3.2e and %3.2e)' %(xf.data.length, hf.data.length)
		exit(-1)		
	# now compute the match -- note that the frequency cutoffs passed here will be 
	# used for setting the limit of the match integral
	match = compute_match (hf, xf, Sh, f_min, f_max, df,  plot_fft=0) 

	if print_steps == 1:
		fp.write('%10.9e\t%10.9e\t%10.9e\t%10.9e\t%10.9e\n'%(m,mchirp,eta,chi,match))
	return match


"""Plot whitened FFTs of target and templates """
def calcmatch_plotFFT(m_targ, eta_targ, m_templ, eta_templ, chi_templ, spin1z_targ, spin2z_targ, incl_angle, phi0, dataFileName, \
sample_rate, f_min, f_max, f_min_wave, f_max_wave, psd,  approx_targ, approx_templ, amplO_templ, amplO_targ, phaseO_templ, \
phaseO_targ, spin1x, spin1y, spin2x, spin2y, sky_theta, sky_phi, polarization_angl, waveform_type, domain_targ, domain_templ,\
N, only_quadrupole, plot_fft, chi_definition, only_optimal_snr, tag):

	################################ generate the target waveform and the psd #############################################################
		
	# the individual masses of the target binary
	m1_targ, m2_targ = individual_masses(m_targ,eta_targ)
	m1_templ, m2_templ = individual_masses(m_templ,eta_templ)

	# target waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain   
	hf = wg.generate_waveforms(m1_targ, m2_targ, spin1x, spin1y, spin1z_targ, spin2x, spin2y, spin2z_targ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sample_rate, N, approx_targ, amplO_targ, phaseO_targ, phi0, dataFileName,waveform_type, domain=domain_targ, use_only_22=only_quadrupole)

	# generate psd
	f = np.linspace(0, sample_rate/2., hf.data.length)   		# frequency axis
	if psd == 'AdvLIGO-ZeroDet-HighP':	# use the fit to AdvLigo Zero detuned high power given in Ajith (2011)
		Sh = adligo_psd(f, f_min)					# PSD
	elif psd == 'WhiteNoise': 
		Sh = np.ones(len(f))
	else: 
		Sh = psds.noise_models[psd](f)
	
	# whiten and normalize the target 
	hf, sigmasq_h = whiten_and_normalize_waveform(hf, Sh, f_min, f_max, hf.deltaF)
	snr_optimal = sigmasq_h**0.5		# optimal signal to noise ratio for the target waveform

	if only_optimal_snr == 1:
		match = 0
		return match, snr_optimal
	
	# down-convert to singla precision
	hf = wf.FrequencySeries_to_COMPLEX8FrequencySeries(hf)

	if chi_definition == 'IMRPhenomB_spin_parameter':
		spin1z_templ = spin2z_templ = chi_templ
	elif chi_definition == 'reduced_spin':
		spin1z_templ = spin2z_templ = chi_templ/(1.-76.*eta_targ/113.)

	# template waveform in Fourier domain
	xf = wg.generate_waveforms(m1_templ, m2_templ, spin1x, spin1y, spin1z_templ, spin2x, spin2y, spin2z_templ, incl_angle,sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sample_rate, N, approx_templ, amplO_templ, phaseO_templ, domain=domain_templ)
 	
	# make sure that the two vectors have the same length, freq resolution etc
	if abs(xf.deltaF- hf.deltaF) < 1e-14:   #FIXME this is a hack 
		df = xf.deltaF
	else: 
		print 'Different frequency sampling for the template and target waveforms (%3.2e and %3.2e)' %(xf.deltaF, hf.deltaF)
		exit(-1) 

	if xf.data.length != hf.data.length:
		print 'Different lengths for the template and target waveform vectors (%3.2e and %3.2e)' %(xf.data.length, hf.data.length)
		exit(-1)		
 
	# now compute the overlap (faithfulness) 
	match = compute_match (hf, xf, Sh, f_min, f_max, df,  plot_fft, tag) 

	return match, snr_optimal



""" Set some parameters based on the input arguments """
def input_arguments(opts):
	
	# the default values for the arguments, so if any argument is not given, the corresponding default value is used
	m_targ= 12.							# Target Mass
 	eta_targ= 0.25							# Target Symmetric mass ration
	chi_templ_min=chi_templ_max=0.0					# Template Reduced Spin Bounds. Unless specified, maximization over \chi is not done and default chi=0.0
	spin1z_targ = 0.						# Target spins of individual compact objects
	spin2z_targ = 0.
	maximization_method='Nelder-Mead'				# Method for optimizing (maximizing) match.
	maximization_parameters = 'mchirp_eta_chi'			# parameters over which match is maximized
	sample_rate = 4096.						# sample rate (not used for hybrid waveforms)
	f_min = 20. 							# low freq cutoff of the detector (Hz) 
	f_max = sample_rate/2. 						# upper frequency cutoff for match calculatin (Hz) 
	psd = 'AdvLIGO-ZeroDet-HighP'					# detector noise 
	spin1x = spin1y = spin2x = spin2y = 0 				# Spins
	approx_targ = 'NRPNHybrid'					# Target waveform
	approx_templ = 'IMRPhenomB'					# Template waveform
	amplO_templ = amplO_targ = 0					# Amplitude PN order
	phaseO_templ = phaseO_targ = 7 					# Phase PN order
	sky_theta = 0							# sky location and orientation of the binary in the detector frame. For optimally 
	sky_phi = 0							# located and oriented binary, set everything to zero 
	polarization_angl = 0
	incl_angle = 0							# inclination angle
	resultsDir = 'NONE'
	plot_fft = 0
	N_startpoints=1							# Number of starting points
	N_cores = 1							# Number of parallel processes	
	error_mchirp=0.1  						# Relative change in starting mchirp and \eta for template.
	error_eta= 0.1							# The maximum of match is expected to occur near mchirp_targ, so we vary \eta by a much larger fraction
	error_chi = 0.1	   						# The maximum is expected to occur near chi = chi_targ_raduced
	m_templ_min = 70.0						# bounds on the total mass and symmetric mass ratio, \eta . 
	m_templ_max = 500.0 						# Note that though we use mchirp as a parameter for maximization, we keep bounds on the total mass,m.
	eta_templ_min = 0.10 
	eta_templ_max = 0.25
	phi0 = 0.0							# phase angle
	tag_file = 'HybData_TaylorT1_3.5PN_BAM_q1.00_spin1[0.00,0.00,0.00]_spin2[0.00,0.00,0.00]'   #name of .mat file		
	print_heading = 0						# heading showing differnt paramter names in the print file
	data_file_location = '/home/ajith/working/cbc/phenom_hh/data/Hybrid_pnOffset/' #location of .mat files
	waveform_type = 'None'						# Type of waveforms, e.g. 'Hyb','PN','NR', use 'None' if using others
	chi_definition = 'IMRPhenomB_spin_parameter'			# Formula for chi: 'reduced_spin' or 'IMRPhenomB_spin_parameter'
	outFileTag = 'output_file'					# Name of the file to which results are written
	print_steps = 0							# print intermediate values of match if value is 1
	print_steps_tag = 'print_intermediate_file'			# Name of the file to which intermediate match values are written
	domain_targ = 'time'						# Domain of the target waveform
	domain_templ = 'freq'						# Domain of the template waveform
	N_multiply= 2							# Factor to multiply length of Hybridwaveform
	only_quadrupole = 0						# use only quadrupole mode if value is 1	
	only_overlap = 0						# calculate only overlap if value is 1	
	only_optimal_snr = 0						# calculate only optimal_snr if value is 1	
	grid_search = 0							# do a grid based search after Nelder-Mead if value is 1
	second_Nelder = 0						# do a second Nelder-Mead search if value is 1

	# if a particular argument is given by the user, store it in the appropriate variable
	for opt, arg in opts:
		if opt == '--m_targ':
		  m_targ=float(arg)
		elif opt == '--eta_targ':
		  eta_targ=float(arg)
		elif opt == '--chi_templ':								#If chi_templ is specified, maximization over chi is not done.
		  chi_templ_min=chi_templ_max= float(arg)						#To optimize over chi_templ specify chi_min and chi_max
		elif opt == '--maximization_method':
		  maximization_method=arg
		elif opt == '--f_min':
		  f_min=float(arg)
		elif opt == '--f_max':
		  f_max=float(arg)
		elif opt == '--sample_rate':
		  sample_rate = float(arg)
		elif opt == '--psd':
		  psd= arg
		elif opt == '--approx_targ':
		  approx_targ= arg
		elif opt == '--approx_templ':
		  approx_templ= arg
		elif opt == '--amplO_templ':
		  amplO_templ=int(arg)
		elif opt == '--amplO_targ':
		  amplO_targ= int(arg)
		elif opt == '--phaseO_templ':
		  phaseO_templ= int(arg)
		elif opt == '--phaseO_targ':
		  phaseO_targ= int(arg)
		elif opt == '--spin1x':
		  spin1x=float(arg)
		elif opt == '--spin1y':
		  spin1y= float(arg)
		elif opt == '--spin2x':
		  spin2x= float(arg)
		elif opt == '--spin2y':
		  spin2y= float(arg)
		elif opt == '--sky_theta':
		  sky_theta= float(arg)
		elif opt == '--sky_phi':
		  sky_phi= float(arg)
		elif opt == '--polarization_angl':
		  polarization_angl= float(arg)
		elif opt == '--incl_angle':
		  incl_angle= float(arg)
		elif opt == '--resultsDir':
		   resultsDir= str(arg)
		elif opt == '--plot_fft':
		   plot_fft= int(arg)
		elif opt == '--N_startpoints':
		   N_startpoints=int(arg)
		elif opt == '--N_cores':
		   N_cores=int(arg)
		elif opt == '--error_mchirp':
		   error_mchirp= float(arg)
		elif opt == '--error_eta':
		   error_eta = float(arg)
		elif opt == '--error_chi':
		   error_chi = float(arg)
		elif opt == '--spin1z_targ':
		  spin1z_targ= float(arg)	
		elif opt == '--spin2z_targ':
		   spin2z_targ=  float(arg)
		elif opt == '--m_templ_min':
		  m_templ_min= float(arg)
		elif opt == '--m_templ_max':
		  m_templ_max= float(arg)
		elif opt == '--eta_templ_min':
		  eta_templ_min= float(arg)
		elif opt == '--eta_templ_max':
		  eta_templ_max= float(arg)
		elif opt == '--chi_templ_min':
		  chi_templ_min= float(arg)
		elif opt == '--chi_templ_max':
		  chi_templ_max= float(arg)
		elif opt == '--phi0':
		  phi0= float(arg)
		elif opt == '--tag_file':
		  tag_file= str(arg)
		elif opt == '--print_heading':
		   print_heading = int(arg)
		elif opt == '--data_file_location':
		  data_file_location= str(arg)
		elif opt == '--waveform_type':
		   waveform_type = str(arg)
		elif opt == '--outFileTag':
		   outFileTag = str(arg)
		elif opt == '--print_steps':
		   print_steps = int(arg)
		elif opt == '--print_steps_tag':
		   print_steps_tag = str(arg)
		elif opt == '--domain_templ':
		   domain_templ = str(arg)
		elif opt == '--domain_targ':
		   domain_targ = str(arg)
		elif opt == '--N_multiply':
		   N_multiply = int(arg)
		elif opt == '--only_quadrupole':
		   only_quadrupole = int(arg)
		elif opt == '--only_overlap':
		   only_overlap = int(arg)
		elif opt == '--only_optimal_snr':
		   only_optimal_snr = int(arg)
		elif opt == '--grid_search':
		   grid_search = int(arg)
		elif opt == '--maximization_parameters':
		   maximization_parameters = str(arg)
		elif opt == '--chi_definition':
		   chi_definition = str(arg)
		elif opt == '--second_Nelder':
		   second_Nelder = int(arg)

	return m_targ, eta_targ, maximization_method, f_min, f_max, psd, sample_rate, approx_targ, approx_templ, amplO_templ, amplO_targ,phaseO_templ,\
		 phaseO_targ, spin1x, spin1y, spin2x, spin2y, spin1z_targ, spin2z_targ, sky_theta, sky_phi, polarization_angl, resultsDir,\
		 plot_fft, N_startpoints, N_cores, error_mchirp, error_eta, error_chi, m_templ_min, m_templ_max, eta_templ_min, eta_templ_max,\
		 chi_templ_min, chi_templ_max,incl_angle,phi0, tag_file, print_heading, data_file_location, waveform_type, outFileTag, print_steps, print_steps_tag,\
		 domain_targ, domain_templ, N_multiply, only_quadrupole, only_overlap, grid_search, maximization_parameters, chi_definition, second_Nelder, only_optimal_snr

""" compute the match between two (one target and one template) waveforms. Note that spin variables are vectors (Sx, Sy, Sz)  """
def calc_match_waveforms(m1_targ, m2_targ, spin1_targ, spin2_targ, m1_templ, m2_templ, spin1_templ, spin2_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, approx_templ, amplO_templ, phaseO_templ, domain_targ, domain_templ, det_name, plot_fft, printW, targ_wave_file_name=0, use22only=0): 

	# unpack spin components 
	spin1x_targ, spin1y_targ, spin1z_targ = spin1_targ 
	spin2x_targ, spin2y_targ, spin2z_targ = spin2_targ 
	spin1x_templ, spin1y_templ, spin1z_templ = spin1_templ 
	spin2x_templ, spin2y_templ, spin2z_templ = spin2_templ 

	# target waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain 
	hf = wg.generate_waveforms(m1_targ, m2_targ, spin1x_targ, spin1y_targ, spin1z_targ, spin2x_targ, spin2y_targ, spin2z_targ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, 0., targ_wave_file_name, 0, domain_targ, use_only_22=use22only, print_waveforms=printW)	
	print '... generated target waveforms (N = %d, df = %f)' %(hf.data.length, hf.deltaF)

	# redefine N and sample_rate so that the tempalte waveform will have the same df. This is probably a hack FIXME 
	sampl_rate = (hf.data.length-1)*hf.deltaF*2. 
	N = 2*(hf.data.length-1)

	# template waveform in the Fourier domain (either by computing the FFT of the time-domain waveform or by generating the waveform directly in the Fourier domain 
	xf = wg.generate_waveforms(m1_templ, m2_templ, spin1x_templ, spin1y_templ, spin1z_templ, spin2x_templ, spin2y_templ, spin2z_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min_wave, f_max_wave, sampl_rate, N, approx_templ, amplO_templ, phaseO_templ, 0., 0, 0, domain_templ, use_only_22=use22only, print_waveforms=printW)	
	print '... generated templates (N = %d, df = %f)' %(xf.data.length, xf.deltaF)

	# make sure that the two vectors have the same length, freq resolution etc 
	if xf.deltaF == hf.deltaF:
		df = hf.deltaF
	else: 
		print 'Different frequency sampling for the target and template waveforms (%e and %e)' %(xf.deltaF, hf.deltaF)
		exit(-1) 
		
	if xf.data.length != hf.data.length:
		print 'Different data length for the target and template waveforms (%d and %d)' %(xf.data.length, hf.data.length)
		exit(-1) 

	# generate psd
	f = np.linspace(0, sampl_rate/2., xf.data.length)
	if det_name == 'WhiteNoise': 
		Sh = np.ones(len(f))
	elif det_name == 'aLIGO_EARLY_HIGH':
		os.system('. ${HOME}/.profile')
		home = os.getenv('HOME')
		psd_path = os.path.join(home, 'src/lalsuite/lalsimulation/src')
		psd = os.path.join(psd_path, 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt')
		f_init, Sh_init = np.loadtxt(psd, unpack=True)
		Sh_interp_object = interp1d(f_init, Sh_init, bounds_error=False, fill_value=0)
		Sh = Sh_interp_object(f)
	else: 
		Sh = psds.noise_models[det_name](f)
	print '... generated psd'

	# whiten and normalize the target waveform  
	hf, sigmasq_h = whiten_and_normalize_waveform(hf, Sh, f_min, f_max, df)
	snr_optimal = sigmasq_h**0.5		# optimal signal to noise ratio for the target waveform

	# down-convert target waveform to single precision
	hf = wf.FrequencySeries_to_COMPLEX8FrequencySeries(hf)

	# now compute the match -- note that the frequency cutoffs passed here will be 
	# used for setting the limit of the match integral 
	tag = '%s_m1_%2.1f_m2_%2.1f_chi1_%3.3f_chi2_%3.3f_Vs_%s_m1_%2.1f_m2_%2.1f_chi1_%3.3f_chi2_%3.3f_%s' %(approx_targ, m1_targ, m2_targ, spin1z_targ, spin2z_targ, approx_templ, m1_templ, m2_templ, spin1z_templ, spin2z_templ, det_name)

	match, sigmasq_x = compute_match (hf, xf, Sh, f_min, f_max, df, plot_fft, tag) 
	#print '#', m1_templ, m2_templ, spin1z_templ, spin2z_templ, snr_optimal, np.sqrt(sigmasq_x), f_min, f_max, int(f_min/df), int(f_max/df), match 
	return match 


""" distribute starting points on a 3D-ellipsoid """
def generate_ellipsoid_points(x_center, y_center, z_center, x_semimajor, y_semimajor,  z_semimajor, N_startpoints):

	# intialize the vectors
	phi_vec = []
	costheta_vec = []


	num_ellip_points = N_startpoints - 1		# first point reserved for the center of the ellipsoid
	# starting points at the north pole #1
	if num_ellip_points > 0: 
		phi_vec.append(0.)
		costheta_vec.append(1.)


	# starting points at the south pole #2
	if num_ellip_points > 1: 
		phi_vec.append(0.)
		costheta_vec.append(-1.)


	# starting points on the equator #6
	if num_ellip_points > 2: 
		n_points_equator = min(num_ellip_points-2, 4)
		phi_vec.extend(np.linspace(0.,2*np.pi,n_points_equator, endpoint=False))
		costheta_vec.extend(0*np.ones(n_points_equator))

	
	# starting points at cos(theta) = 0.5 #9
	if num_ellip_points > 6: 
		n_points_latitude = min(num_ellip_points-6, 3)
		phi_vec.extend(np.linspace(0.,2*np.pi,n_points_latitude, endpoint=False) )
		costheta_vec.extend(0.5* np.ones(n_points_latitude))

	# starting points at cos(theta) = -0.5 #12
	if num_ellip_points > 9: 
		n_points_latitude = min(num_ellip_points-6, 3)
		phi_vec.extend(np.linspace(0.,2*np.pi,n_points_latitude, endpoint=False) )
		costheta_vec.extend(-0.5* np.ones(n_points_latitude))


	theta_vec = np.arccos(costheta_vec)

	# initialze points array
	ellipsoid_points = np.zeros([3,N_startpoints])    	# The first set is reserved for x_center,y_center,z_center
		
	for k in range(N_startpoints-1):
	    
	    theta = theta_vec[k]
	    phi = phi_vec[k]	
	    x = x_semimajor*np.cos(phi)*np.sin(theta)		# distance of point along x-axis from center
	    y = y_semimajor*np.sin(phi)*np.sin(theta)		# distance of point along y-axis from center
	    z =	z_semimajor*np.cos(theta)			# distance of point along z-axis from center
	    ellipsoid_points[0,k+1]=x+x_center			# The x value of point on the surface of ellipsoid
	    ellipsoid_points[1,k+1]=y+y_center			# The y value of point on the surface of ellipsoid
	    ellipsoid_points[2,k+1]=z+z_center			# The z value of point on the surface of ellipsoid
		
	ellipsoid_points[0,0]=x_center				# The first set for x_center,y_center,z_center
	ellipsoid_points[1,0]=y_center
	ellipsoid_points[2,0]=z_center

	return ellipsoid_points 



""" distribute starting points on a 2D-ellipse """
#def generate_ellipse_points(x_center, y_center, z_center, x_semimajor, y_semimajor,  z_semimajor, N_startpoints, rot_angle=np.pi/4.):
def generate_ellipse_points(x_center, y_center, z_center, error_x, error_z, N_startpoints, rot_angle=0):

	# intialize the vectors
	phi_vec = []
	costheta_vec = []
	num_ellip_points = N_startpoints - 1		# first point reserved for the center of the ellipse
	angle = np.linspace(0,2*np.pi,num_ellip_points,  endpoint=False)

	# distribute the points on an ellipse
	x_err_vals = error_x*np.cos(angle)
	z_err_vals = error_z*np.sin(angle)

	## rotate the ellipse if required
	if rot_angle != 0:
		## use temp varialbes to avoid overwriting
		x_err_temp = x_err_vals * np.cos(rot_angle) + z_err_vals * np.sin(rot_angle)
		z_err_temp = -x_err_vals * np.sin(rot_angle) + z_err_vals * np.cos(rot_angle)
		x_err_vals = x_err_temp
		z_err_vals = z_err_temp

	# initialze points array
	ellipse_points = np.zeros([3,N_startpoints])    

	# The first set for x_center,y_center,z_center
	ellipse_points[0,0]=x_center				
	ellipse_points[1,0]=y_center
	ellipse_points[2,0]=z_center

	# distribute the points on an ellipse
	ellipse_points[0,1:]=x_center*(np.ones(len(x_err_vals)) + x_err_vals)
	ellipse_points[1,1:]=y_center*np.ones(len(x_err_vals)) 				# chi_start is always chi_center in this case
	ellipse_points[2,1:]=z_center*(np.ones(len(x_err_vals)) + z_err_vals) 

	return ellipse_points 
