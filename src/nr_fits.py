""" 
Various fitting formulas provided by numerical relativity 

P. Ajith, 2015-04-09
N. K. Johnson-McDaniel, 2015-10-09 (added simple recoil fit) 
"""

import numpy as np 
import scipy.optimize as so
#from scipy.integrate import dblquad, nquad # Comment out before committing...
#from skmonaco import mcquad, mcmiser
#Removing these packages, since they're not installed on some machines. This means that the bbh_recoil_simple_averaged* functions will not work, except for the analytically averaged version in the ...system_frame one.

def calc_isco_radius(a): 
	""" 
	Calculate the ISCO radius of a Kerr BH as a function of the Kerr parameter 
	
	a : Kerr parameter (numpy array) 
	"""

	# Ref. Eq. (2.5) of Ori, Thorne Phys Rev D 62 124022 (2000)
	z1 = 1.+(1.-a**2.)**(1./3)*((1.+a)**(1./3) + (1.-a)**(1./3))
	z2 = np.sqrt(3.*a**2 + z1**2)
	a_sign = np.sign(a)
	return 3+z2 - np.sqrt((3.-z1)*(3.+z1+2.*z2))*a_sign

def calc_isco_freq(a): 
	""" 
	Calculate the ISCO frequency of a Kerr BH as a function of the Kerr parameter 
	
	a : Kerr parameter (numpy array) 
	"""
	
	r_isco = calc_isco_radius(a)
	u_isco = r_isco**-0.5	
	v_isco = u_isco*(1.-a*u_isco**3.+a**2.*u_isco**6)**(1./3.)
	return v_isco**3./np.pi
	
def _final_spin_diff(a_f, eta, delta_m, S, Delta): 
	""" Internal function: the final spin is determined by minimizing this function """
	
	# calculate ISCO radius 
	r_isco = calc_isco_radius(a_f)
	
	# angular momentum at ISCO -- Eq.(2.8) of Ori, Thorne Phys Rev D 62 124022 (2000)
	J_isco = (3*np.sqrt(r_isco)-2*a_f)*2./np.sqrt(3*r_isco)

	# fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004 (2014) 
	# [forth order fits]
	L0  = 0.686710
	L1  = 0.613247 
	L2a = -0.145427
	L2b = -0.115689
	L2c = -0.005254
	L2d = 0.801838 
	L3a = -0.073839
	L3b = 0.004759 
	L3c = -0.078377
	L3d = 1.585809 
	L4a = -0.003050
	L4b = -0.002968
	L4c = 0.004364 
	L4d = -0.047204
	L4e = -0.053099
	L4f = 0.953458 
	L4g = -0.067998
	L4h = 0.001629 
	L4i = -0.066693

	a_f_new = (4.*eta)**2.*(L0  +  L1*S +  L2a*Delta*delta_m + L2b*S**2. + L2c*Delta**2 \
		+ L2d*delta_m**2. + L3a*Delta*S*delta_m + L3b*S*Delta**2. + L3c*S**3. \
		+ L3d*S*delta_m**2. + L4a*Delta*S**2*delta_m + L4b*Delta**3.*delta_m \
		+ L4c*Delta**4. + L4d*S**4. + L4e*Delta**2.*S**2. + L4f*delta_m**4 + L4g*Delta*delta_m**3. \
		+ L4h*Delta**2.*delta_m**2. + L4i*S**2.*delta_m**2.) \
		+ S*(1. + 8.*eta)*delta_m**4. + eta*J_isco*delta_m**6.

	return np.ravel(abs(a_f-a_f_new))


def bbh_final_mass_and_spin_non_precessing(m1, m2, chi1, chi2): 
	""" 
	Calculate the mass and spin of the final BH resulting from the 
	merger of two black holes with non-precessing spins

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_non_precessing)(m1, m2, chi1, chi2)
	
	# binary parameters 
	m = m1+m2  
	q = m1/m2 
	eta = q/(1.+q)**2. 
	delta_m = (m1-m2)/m

	S1 = chi1*m1**2				# spin angular momentum 1 
	S2 = chi2*m2**2				# spin angular momentum 2 
	S = (S1+S2)/m**2 			# symmetric spin (dimensionless -- called \tilde{S} in the paper) 
	Delta = (S2/m2-S1/m1)/m		# antisymmetric spin (dimensionless -- called tilde{Delta} in the paper

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# compute the final spin 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# res = so.minimize_scalar(_final_spin_diff, bounds=(-0.99999, 0.99999), args=(eta, delta_m, S, Delta), method='Bounded', tol=1e-6, options={'maxiter':100, 'disp':False})
	# a_f = res.x

	x, cov_x = so.leastsq(_final_spin_diff, 0., args=(eta, delta_m, S, Delta))
	a_f = x[0]

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# now compute the final mass  
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	r_isco = calc_isco_radius(a_f)

	# fitting coefficients - Table X1 of Healy et al Phys Rev D 90, 104004 (2014) 
	# [forth order fits]
	M0  = 0.951507 
	K1  = -0.051379
	K2a = -0.004804
	K2b = -0.054522
	K2c = -0.000022
	K2d = 1.995246 
	K3a = 0.007064 
	K3b = -0.017599
	K3c = -0.119175
	K3d = 0.025000 
	K4a = -0.068981
	K4b = -0.011383
	K4c = -0.002284
	K4d = -0.165658
	K4e = 0.019403 
	K4f = 2.980990 
	K4g = 0.020250 
	K4h = -0.004091
	K4i = 0.078441 

	# binding energy at ISCO -- Eq.(2.7) of Ori, Thorne Phys Rev D 62 124022 (2000)
	E_isco = (1. - 2./r_isco + a_f/r_isco**1.5)/np.sqrt(1. - 3./r_isco + 2.*a_f/r_isco**1.5)

	# final mass -- Eq. (14) of Healy et al Phys Rev D 90, 104004 (2014)
	mf = (4.*eta)**2*(M0 + K1*S + K2a*Delta*delta_m + K2b*S**2 + K2c*Delta**2 + K2d*delta_m**2 \
		+ K3a*Delta*S*delta_m + K3b*S*Delta**2 + K3c*S**3 + K3d*S*delta_m**2 \
		+ K4a*Delta*S**2*delta_m + K4b*Delta**3*delta_m + K4c*Delta**4 + K4d*S**4 \
		+ K4e*Delta**2*S**2 + K4f*delta_m**4 + K4g*Delta*delta_m**3 + K4h*Delta**2*delta_m**2 \
		+ K4i*S**2*delta_m**2) + (1+eta*(E_isco+11.))*delta_m**6.
	
	return mf*m, a_f  

def bbh_final_mass_and_spin_HLZ_extension_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Healy, Lousto, and Zlochower fit for the aligned part
	and including the contribution from the in-plane spin components

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_HLZ_extension_precessing)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# First compute the final mass and parallel component of the final spin using the aligned components of the initial spins
	Mf, afpara = bbh_final_mass_and_spin_non_precessing(m1, m2, chi1*np.cos(tilt1), chi2*np.cos(tilt2))

	# Now compute the magnitude of the in-plane final dimensionful spin, first computing the magnitudes of the initial in-plane spins

	S1perpmag = m1*m1*chi1*np.sin(tilt1)
	S2perpmag = m2*m2*chi2*np.sin(tilt2)

	Sperpmag = (S1perpmag*S1perpmag + S2perpmag*S2perpmag + 2.*S1perpmag*S2perpmag*np.cos(phi12))**0.5

	# Put everything together

	return Mf, (afpara*afpara + Sperpmag*Sperpmag/Mf**4.)**0.5

def qnmfreqs_berti(a, l, m, n):
	"""
	compute the (complex) QNM frequencies for different 
	overtones using the fits provided by Berti, Cardoso, Will 
	
	a     	: Kerr parameter of the BH 
	l, m, n	: indices of spheroidal harnonics modes 
	"""

	# load the data file containing the fits (Berti et al,  gr-qc/0512160)
	lVec, mVec, nVec, f1Vec, f2Vec, f3Vec, q1Vec, q2Vec, q3Vec = np.loadtxt('../src/Berti_QNMfitcoeffsWEB.dat', unpack=True)

	idx = np.logical_and(np.logical_and(lVec == l, mVec == m), nVec == n)

	# evaluate the Berti et al fits to the complex frequencies 
	if len(lVec[idx]) == 1:

		 f1 = f1Vec[idx]
		 f2 = f2Vec[idx]
		 f3 = f3Vec[idx]
		 q1 = q1Vec[idx]
		 q2 = q2Vec[idx]
		 q3 = q3Vec[idx]

		 omega = (f1 + f2*(1.-a)**f3) # fit to omega
		 Q = q1 + q2*(1.-a)**q3       # fit to quality factor 
		 tau = 2.*Q/omega

		 return Q, omega - 1j/tau  # complex frequency 

	elif len(lVec[idx]) < 1:
			print '# No data matching this l, m, n combination (l = #d m = #d n = #d)' %(l, m, n)
			exit()
	else:
			print '# More than on fit point corresponding to this l, m, n combination'
			exit()

def calc_fqnm_dominant_mode(af): 
	"""
	calculate the (dimensionless) freq of the dominant QNM of a Kerr BH with spin af 
	"""

	Q, Omega = qnmfreqs_berti(af, 2, 2, 0)
	return Q, np.real(Omega)/(2*np.pi) 

def bbh_recoil_simple_aligned(q, chi1para, chi2para):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with aligned spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013)

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1para: the component of the dimensionless spin of the *lighter* black hole along the angular momentum (z)
	chi2para: the component of the dimensionless spin of the *heavier* black hole along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_aligned)(q, chi1para, chi2para)
	
	# Check that the spins are at most extremal

	# Check that the spins are at most extremal

	if abs(chi1para) > 1:
		print("WARNING: aligned component of spin 1 has a magnitude of %3f\n > 1" %abs(chi1para))
	if abs(chi2para) > 1:
		print("WARNING: aligned component of spin 2 has a magnitude of %3f\n > 1" %abs(chi2para))

	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2para, chi1para = chi1para, chi2para

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# Define the relevant constants from the fit
	Am = 12000.
	Bm = -0.93
	H = 6900. # The quoted error bar is 500 
	cxi = -0.819152 # cos of xi

	# Compute the two contributions to the recoil. Here we take H_S = 0, since there's not a value given in the paper.

	Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
	Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)

	return (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2.)**0.5

def bbh_recoil_simple_averaged_over_angles(q, chi1perpx, chi1perpy, chi1para, chi2perpx, chi2perpy, chi2para):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with generic spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013), averaged over the angles present in the in-plane spins portion of the fit

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1perpx, chi1perpy, chi1para: the components of the dimensionless spin of the *lighter* black hole in the orbital plane (x and y) and along the angular momentum (z)
	chi2perpx, chi2perpy, chi2para: the components of the dimensionless spin of the *heavier* black hole in the orbital plane (x and y) and along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_averaged_over_angles)(q, chi1perpx, chi1perpy, chi1para, chi2perpx, chi2perpy, chi2para)
	
	# Check that the spins are at most extremal

	if (chi1perpx**2. + chi1perpy**2. + chi1para**2.) > 1:
		print("WARNING: spin 1 has a magnitude of %3f\n > 1" %(chi1perpx**2. + chi1perpy**2. + chi1para**2.)**0.5)
	if (chi2perpx**2. + chi2perpy**2. + chi2para**2.) > 1:
		print("WARNING: spin 2 has a magnitude of %3f\n > 1" %(chi2perpx**2. + chi2perpy**2. + chi2para**2.)**0.5)


	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2perpx, chi1perpx = chi1perpx, chi2perpx
	   chi2perpy, chi1perpy = chi1perpy, chi2perpy
	   chi2para, chi1para = chi1para, chi2para


	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# Define the relevant constants from the fit
	Am = 12000.
	Bm = -0.93
	H = 6900. # The quoted error bar is 500 
	K = 59000. # The quoted error bar is 1000
	KS = -4.254
	cxi = -0.819152 # cos of xi

	# Compute the three contributions to the recoil, starting with two norms needed for the parallel component. Here we take H_S = 0, since there's not a value given in the paper.

	norm1 = ((chi2perpx - q*chi1perpx)**2. + (chi2perpy - q*chi1perpy)**2.)**0.5
	norm2 = ((chi2perpx + q**2.*chi1perpx)**2. + (chi2perpy + q**2.*chi1perpy)**2.)**0.5

	Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
	Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)
	Vpara1 = K*(eta**2./(1.+q))*norm1
	Vpara2 = K*(eta**2./(1.+q))*KS*((1.-q)/(1.+q)**2.)*norm2

	Vmagavg = nquad(lambda psiD, psiS: (1./(4.*np.pi**2.))*(Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + (Vpara1*np.cos(psiD) + Vpara2*np.cos(psiS))**2.)**0.5, [[0, 2.*np.pi], [0, 2.*np.pi]])

	# Note that psiD := \phi_\Delta - \phi_1 and psiS := \phi_S - \phi_2, where the expressions on the right are the ones from the paper

	# Check that the error estimate for the integration is appropriately small, and complain if it is not

	if Vmagavg[1]/Vmagavg[0] > 1.e-6:
		print("WARNING: Fractional error is larger than anticipated, at: %3f\n" %(Vmagavg[1]/Vmagavg[0]))

	return Vmagavg[0]

def bbh_recoil_simple_incl_angles_system_frame(q, chi1, chi2, ctilt1, ctilt2, phi12, psiD, psiS):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with generic spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013), averaged over the angles present in the in-plane spins portion of the fit

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1, ctilt1: magnitude of the dimensionless spin of the *lighter* black hole and the cosine of its angle with the orbital angular momentum
	chi2, ctilt2: magnitude of the dimensionless spin of the *heavier* black hole and the cosine of its angle with the orbital angular momentum
	phi12: angle between the in-plane spins
	phiD: angle \phi_\Delta - \phi_1 used in the in-plane spin portion of the recoil fit
	phiS: angle \phi_S - \phi_2 used in the in-plane spin portion of the recoil fit
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_incl_angles_system_frame)(q, chi1, chi2, ctilt1, ctilt2, phi12, psiD, psiS)
	
	# Check that the spins are at most extremal

	if chi1 > 1:
		print("WARNING: spin 1 has a magnitude of %3f\n > 1" %chi1)
	if chi2 > 1:
		print("WARNING: spin 2 has a magnitude of %3f\n > 1" %chi2)


	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2, chi1 = chi1, chi2
	   ctilt2, ctilt1 = ctilt1, ctilt2

	# Computing the spins parallel and perpendicular to the orbital angular momentum
	chi1para = chi1*ctilt1
	chi2para = chi2*ctilt2
	chi1perp = chi1*(1 - ctilt1**2.)**0.5
	chi2perp = chi2*(1 - ctilt2**2.)**0.5

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# Define the relevant constants from the fit
	Am = 12000.
	Bm = -0.93
	H = 6900. # The quoted error bar is 500 
	K = 59000. # The quoted error bar is 1000
	KS = -4.254
	cxi = -0.819152 # cos of xi

	# Compute the three contributions to the recoil, starting with two norms needed for the parallel component. Here we take H_S = 0, since there's not a value given in the paper.

	cp12 = np.cos(phi12)

	norm1 = (chi2perp**2. - 2.*q*cp12*chi1perp*chi2perp + q**2.*chi1perp**2.)**0.5
	norm2 = (chi2perp**2. + 2.*q**2.*cp12*chi1perp*chi2perp + q**4.*chi1perp**2.)**0.5

	Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
	Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)
	Vpara1 = K*(eta**2./(1.+q))*norm1
	Vpara2 = K*(eta**2./(1.+q))*KS*((1.-q)/(1.+q)**2.)*norm2

	Vmagavg = (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + (Vpara1*np.cos(psiD) + Vpara2*np.cos(psiS))**2.)**0.5 
	
	return Vmagavg

def bbh_recoil_simple_averaged_over_angles_system_frame(q, chi1, chi2, ctilt1, ctilt2, phi12):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with generic spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013), averaged over the angles present in the in-plane spins portion of the fit

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1, ctilt1: magnitude of the dimensionless spin of the *lighter* black hole and the cosine of its angle with the orbital angular momentum
	chi2, ctilt2: magnitude of the dimensionless spin of the *heavier* black hole and the cosine of its angle with the orbital angular momentum
	phi12: angle between the in-plane spins
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_averaged_over_angles_system_frame)(q, chi1, chi2, ctilt1, ctilt2, phi12)
	
	# Check that the spins are at most extremal

	if chi1 > 1:
		print("WARNING: spin 1 has a magnitude of %3f\n > 1" %chi1)
	if chi2 > 1:
		print("WARNING: spin 2 has a magnitude of %3f\n > 1" %chi2)


	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2, chi1 = chi1, chi2
	   ctilt2, ctilt1 = ctilt1, ctilt2

	# Computing the spins parallel and perpendicular to the orbital angular momentum
	chi1para = chi1*ctilt1
	chi2para = chi2*ctilt2
	chi1perp = chi1*(1 - ctilt1**2.)**0.5
	chi2perp = chi2*(1 - ctilt2**2.)**0.5

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# Define the relevant constants from the fit
	Am = 12000.
	Bm = -0.93
	H = 6900. # The quoted error bar is 500 
	K = 59000. # The quoted error bar is 1000
	KS = -4.254
	cxi = -0.819152 # cos of xi

	# Compute the three contributions to the recoil, starting with two norms needed for the parallel component. Here we take H_S = 0, since there's not a value given in the paper.

	cp12 = np.cos(phi12)

	norm1 = (chi2perp**2. - 2.*q*cp12*chi1perp*chi2perp + q**2.*chi1perp**2.)**0.5
	norm2 = (chi2perp**2. + 2.*q**2.*cp12*chi1perp*chi2perp + q**4.*chi1perp**2.)**0.5

	Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
	Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)
	Vpara1 = K*(eta**2./(1.+q))*norm1
	Vpara2 = K*(eta**2./(1.+q))*KS*((1.-q)/(1.+q)**2.)*norm2

	#Vmagavg = nquad(lambda psiD, psiS: (1./(4.*np.pi**2.))*(Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + (Vpara1*np.cos(psiD) + Vpara2*np.cos(psiS))**2.)**0.5, [[0, 2.*np.pi], [0, 2.*np.pi]])

	#Vmagavg = dblquad(lambda psiD, psiS: (1./(4.*np.pi**2.))*(Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + (Vpara1*np.cos(psiD) + Vpara2*np.cos(psiS))**2.)**0.5, 0, 2.*np.pi, lambda phiS: 0, lambda phiD: 2.*np.pi)

	# Note that psiD := \phi_\Delta - \phi_1 and psiS := \phi_S - \phi_2, where the expressions on the right are the ones from the paper

	# Check that the error estimate for the integration is appropriately small, and complain if it is not

	#if Vmagavg[1]/Vmagavg[0] > 1.e-6:
	#	print("WARNING: Fractional error is larger than anticipated, at: %3f\n" %(Vmagavg[1]/Vmagavg[0]))

	#return Vmagavg[0]

	Vmagavg = (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + (Vpara1**2. + Vpara2**2.)/2)**0.5 # Version where we analytically average over the angles in Vpara**2 before taking the square root
	
	return Vmagavg

def bbh_recoil_simple_averaged_over_perp_spins(q, chi1para, chi2para):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with generic spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013), averaged over the in-plane spins, constrained by the input aligned component (as well as the angles present in the in-plane spin portion of the fit)

	N.B.: Does not work for exactly extremal spins, due to division by zero. (This could be fixed by giving appropriate cases...)

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1perpx, chi1perpy, chi1para: the components of the dimensionless spin of the *lighter* black hole in the orbital plane (x and y) and along the angular momentum (z)
	chi2perpx, chi2perpy, chi2para: the components of the dimensionless spin of the *heavier* black hole in the orbital plane (x and y) and along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_averaged_over_perp_spins)(q, chi1para, chi2para)
	
	# Check that the spins are at most extremal

	if abs(chi1para) > 1:
		print("WARNING: aligned component of spin 1 has a magnitude of %3f\n > 1" %abs(chi1para))
	if abs(chi2para) > 1:
		print("WARNING: aligned component of spin 2 has a magnitude of %3f\n > 1" %abs(chi2para))

	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2para, chi1para = chi1para, chi2para

	# Compute the maximum possible magnitude of the in-plane spin, given the aligned component
	   
	chi1perpmax = (1. - chi1para**2.)**0.5
	chi2perpmax = (1. - chi2para**2.)**0.5

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


	def Vmagintegrand(args):

		chi1perp, theta1, chi2perp, theta2, psiD, psiS = args

		# Define the relevant constants from the fit
		Am = 12000.
		Bm = -0.93
		H = 6900. # The quoted error bar is 500 
		K = 59000. # The quoted error bar is 1000
		KS = -4.254
		cxi = -0.819152 # cos of xi

		# Compute the three contributions to the recoil, starting with two norms needed for the parallel component. Here we take H_S = 0, since there's not a value given in the paper.		

		chi1perpx = chi1perp*np.cos(theta1)
		chi1perpy = chi1perp*np.sin(theta1)

		chi2perpx = chi2perp*np.cos(theta2)
		chi2perpy = chi2perp*np.sin(theta2)

		norm1 = ((chi2perpx - q*chi1perpx)**2. + (chi2perpy - q*chi1perpy)**2.)**0.5
		norm2 = ((chi2perpx + q**2.*chi1perpx)**2. + (chi2perpy + q**2.*chi1perpy)**2.)**0.5

		Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
		Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)
		Vpara = K*(eta**2./(1.+q))*(norm1*np.cos(psiD) + KS*((1.-q)/(1.+q)**2.)*norm2*np.cos(psiS))

		# Note that psiD := \phi_\Delta - \phi_1 and psiS := \phi_S - \phi_2, where the expressions on the right are the ones from the paper

		return (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + Vpara**2.)**0.5*chi1perp*chi2perp

	Vmagavg = mcquad(Vmagintegrand, npoints=20000, xl=[0., 0., 0., 0., 0., 0.], xu=[chi1perpmax, 2.*np.pi, chi2perpmax, 2.*np.pi, 2.*np.pi, 2.*np.pi])

	#Vmagavg = mcmiser(Vmag, npoints=500000, xl=[0., 0., 0., 0., 0., 0.], xu=[chi1perpmax, 2.*np.pi, chi2perpmax, 2.*np.pi, 2.*np.pi, 2.*np.pi])

	# Check that the error estimate for the integration is appropriately small, and complain if it is not

	if Vmagavg[1]/Vmagavg[0] > 2.e-2:
		print("WARNING: Fractional error is larger than anticipated, at: %3f\n" %(Vmagavg[1]/Vmagavg[0]))

	return Vmagavg[0]/(4*np.pi**4*chi1perpmax**2.*chi2perpmax**2.) #, Vmagavg[1]//(4*np.pi**4*chi1perpmax**2.*chi2perpmax**2.)


def bbh_recoil_simple_averaged_over_perp_spins_weighted(q, chi1para, chi2para, sigmachi):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with generic spins using the RIT
	group's simple fit given in, e.g., Eqs. (1) and (2) of Lousto and Zlochower, PRD 87, 084027 (2013), averaged over the in-plane spins, constrained by the input aligned component (as well as the angles present in the in-plane spin portion of the fit)

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1perpx, chi1perpy, chi1para: the components of the dimensionless spin of the *lighter* black hole in the orbital plane (x and y) and along the angular momentum (z)
	chi2perpx, chi2perpy, chi2para: the components of the dimensionless spin of the *heavier* black hole in the orbital plane (x and y) and along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_averaged_over_perp_spins_weighted)(q, chi1para, chi2para, sigmachi)
	
	# Check that the spins are at most extremal

	if abs(chi1para) > 1:
		print("WARNING: aligned component of spin 1 has a magnitude of %3f\n > 1" %abs(chi1para))
	if abs(chi2para) > 1:
		print("WARNING: aligned component of spin 2 has a magnitude of %3f\n > 1" %abs(chi2para))

	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2para, chi1para = chi1para, chi2para

	# Compute the maximum possible magnitude of the in-plane spin, given the aligned component
	   
	chi1perpmax = (1. - chi1para**2.)**0.5
	chi2perpmax = (1. - chi2para**2.)**0.5

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


	def Vmagintegrand(args):

		chi1perp, theta1, chi2perp, theta2, psiD, psiS = args

		# Define the relevant constants from the fit
		Am = 12000.
		Bm = -0.93
		H = 6900. # The quoted error bar is 500 
		K = 59000. # The quoted error bar is 1000
		KS = -4.254
		cxi = -0.819152 # cos of xi

		# Compute the three contributions to the recoil, starting with two norms needed for the parallel component. Here we take H_S = 0, since there's not a value given in the paper.		

		chi1perpx = chi1perp*np.cos(theta1)
		chi1perpy = chi1perp*np.sin(theta1)

		chi2perpx = chi2perp*np.cos(theta2)
		chi2perpy = chi2perp*np.sin(theta2)

		norm1 = ((chi2perpx - q*chi1perpx)**2. + (chi2perpy - q*chi1perpy)**2.)**0.5
		norm2 = ((chi2perpx + q**2.*chi1perpx)**2. + (chi2perpy + q**2.*chi1perpy)**2.)**0.5

		Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)
		Vperp = H*(eta**2./(1.+q))*(chi2para - q*chi1para)
		Vpara = K*(eta**2./(1.+q))*(norm1*np.cos(psiD) + KS*((1.-q)/(1.+q)**2.)*norm2*np.cos(psiS))

		# Note that psiD := \phi_\Delta - \phi_1 and psiS := \phi_S - \phi_2, where the expressions on the right are the ones from the paper

		return (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2. + Vpara**2.)**0.5*np.exp(-(chi1perp/sigmachi)**2./2.)*np.exp(-(chi2perp/sigmachi)**2./2.)*chi1perp*chi2perp

	Vmagavg = mcquad(Vmagintegrand, npoints=20000, xl=[0., 0., 0., 0., 0., 0.], xu=[chi1perpmax, 2.*np.pi, chi2perpmax, 2.*np.pi, 2.*np.pi, 2.*np.pi])

	#Vmagavg = mcmiser(Vmag, npoints=500000, xl=[0., 0., 0., 0., 0., 0.], xu=[chi1perpmax, 2.*np.pi, chi2perpmax, 2.*np.pi, 2.*np.pi, 2.*np.pi])

	# Check that the error estimate for the integration is appropriately small, and complain if it is not

	if Vmagavg[1]/Vmagavg[0] > 4.e-2:
		print("WARNING: Fractional error is larger than anticipated, at: %3f\n" %(Vmagavg[1]/Vmagavg[0]))

	return Vmagavg[0]/(16*np.pi**4*(1. - np.exp(-(chi1perpmax/sigmachi)**2./2.))*(1. - np.exp(-(chi2perpmax/sigmachi)**2./2.))*sigmachi**4.) #, Vmagavg[1]/(16*np.pi**4*(1. - np.exp(-(chi1perpmax/sigmachi)**2./2.))*(1. - np.exp(-(chi2perpmax/sigmachi)**2./2.))*sigmachi**4.)


def bbh_recoil_HLZ_aligned(q, chi1para, chi2para):
	""" 
	Calculate the magnitude of the recoil of the final BH resulting from the 
	merger of two black holes with aligned spins using the RIT
	group's more accurate fit given in Eqs. (2A), (12), (17), and (18) (and parameter values from Table VI) of Healy, Lousto, and Zlochower, PRD 90, 104004 (2013). (Note that the nonspinning part of the fit is the same as in bbh_recoil_simple_aligned.)

	q: mass ratio (here m1/m2--the code automatically converts to the q < 1 convention required by the fit)
	chi1para: the component of the dimensionless spin of the *lighter* black hole along the angular momentum (z)
	chi2para: the component of the dimensionless spin of the *heavier* black hole along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_recoil_simple_aligned)(q, chi1para, chi2para)
	
	# Check that the spins are at most extremal

	# Check that the spins are at most extremal

	if abs(chi1para) > 1:
		print("WARNING: aligned component of spin 1 has a magnitude of %3f\n > 1" %abs(chi1para))
	if abs(chi2para) > 1:
		print("WARNING: aligned component of spin 2 has a magnitude of %3f\n > 1" %abs(chi2para))

	# First compute eta and then convert the mass ratio to the \leq 1 convention used in the fit, if it's not already in it. Also reverse the spins to correspond to this convention, if necessary
	eta = q/(1.+q)**2. 
	if q > 1.:
	   q = 1./q
	   chi2para, chi1para = chi1para, chi2para
	dm = (q -1.)/(q + 1.) # dm is \delta m

	dm2 = dm*dm
	dm3 = dm2*dm
	
	# Now compute the spin combinations used in the fit
	Stpara = (chi2para + chi1para*q**2.)/(1. + q)**2.
	Deltatpara = (chi2para - q*chi1para)/(1. + q)

	Stpara2 = Stpara*Stpara
	Stpara3 = Stpara2*Stpara

	Deltatpara2 = Deltatpara*Deltatpara
	Deltatpara3 = Deltatpara2*Deltatpara

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# Compute the magnitude of the recoil 
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# Define the relevant constants from the fit
	Am = 12000.
	Bm = -0.93
	H = 7367.250029 # \pm 66.122336
	H2a = -1.626094 # \pm 0.053888
	H2b = -0.578177 # \pm 0.055790
	H3a = -0.717370 # \pm 0.077605
	H3b = -2.244229 # \pm 0.137982
	H3c = -1.221517 # \pm 0.176699
	H3d = -0.002325 # \pm 0.021612
	H3e = -1.064708 # \pm 0.133021
	H4a = -0.579599 # \pm 0.297351
	H4b = -0.455986 # \pm 0.302432
	H4c =  0.010963 # \pm 0.174289
	H4d =  1.542924 # \pm 0.274459
	H4e = -4.735367 # \pm 0.430869 
	H4f = -0.284062 # \pm 0.174087

	a   =  2.611988 # \pm 0.028327
	b   =  1.383778 # \pm 0.092915
	c   =  0.549758 # \pm 0.113300

	# Compute the cosine of xi

	cxi = np.cos(a + b*Stpara + c*dm*Deltatpara)

	# Compute the two contributions to the recoil. Here we take H_S = 0, since there's not a value given in the paper.

	Vm = Am*eta**2.*((1.-q)/(1.+q))*(1 + Bm*eta)

	VperpH2 = H2a*Stpara*dm + H2b*Deltatpara*Stpara
	VperpH3 = H3a*Deltatpara2*dm + H3b*Stpara2*dm + H3c*Deltatpara*Stpara2 + H3d*Deltatpara3 + H3e*Deltatpara*dm2
	VperpH4 = H4a*Stpara*Deltatpara2*dm + H4b*Stpara3*dm + H4c*Stpara*dm3 + H4d*Deltatpara*Stpara*dm2 + H4e*Deltatpara*Stpara3 + H4f*Stpara*Deltatpara3

	Vperp = H*eta**2*(Deltatpara + VperpH2 + VperpH3 + VperpH4)

	return (Vm**2. + 2.*Vm*Vperp*cxi + Vperp**2.)**0.5

def recoil_sys_frame_test(q, chi1perpx, chi1perpy, chi1para, chi2perpx, chi2perpy, chi2para):

	chi1 = (chi1perpx**2. + chi1perpy**2. + chi1para**2.)**0.5
	chi2 = (chi2perpx**2. + chi2perpy**2. + chi2para**2.)**0.5

	ctilt1 = chi1para/chi1
	ctilt2 = chi2para/chi2

	phi12 = np.arccos((chi1perpx*chi2perpx + chi1perpy*chi2perpy)/(chi1perpx**2. + chi1perpy**2.)**0.5/(chi2perpx**2. + chi2perpy**2.)**0.5)

	vr_orig = bbh_recoil_simple_averaged_over_angles(q, chi1perpx, chi1perpy, chi1para, chi2perpx, chi2perpy, chi2para)
	vr_new = bbh_recoil_simple_averaged_over_angles_system_frame(q, chi1, chi2, ctilt1, ctilt2, phi12)

	return vr_orig, vr_new

def bbh_aligned_Lpeak22_LL(q, chi1para, chi2para):
	""" 
	Calculate the peak luminosity (currently just due to the 2,2 and 2,-2 modes) in geometrized units of a binary black hole with aligned spins using the fit made by Lionel London

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak22_LL)(q, chi1para, chi2para)
	
	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	# Evaluate Lionel's fit for rh22dot
	
	l = [ 	0.371321246247, 1.066104836526, 0.165120075182, 0.081229510674]
	rh22dot = eta*( l[1]*eta + l[2]*chi_eff + l[3]*chi_eff*chi_eff + l[0] )

	# Evaluate Eq. (3.8) of Ruiz et al. GRG 40, 2467 (2008), arXiv:0707.4654v3 to compute the peak luminosity; currently just for the 2,2 and 2,-2 modes (hence the overall factor of 2)

	return 2.*(rh22dot**2.) / (16.0*np.pi)


def bbh_aligned_Lpeak_lmax5_LL(q, chi1para, chi2para):
	""" 
	Calculate the peak luminosity (using the modes up to \ell = 5) in geometrized units of a binary black hole with aligned spins using the fit made by Lionel London

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak_lmax5_LL)(q, chi1para, chi2para)
	
	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	# Evaluate Lionel's fit for the square root of the luminosity
	
	l = [ 0.105937217783, 0.091694475457, 0.032953100645, 0.019520946419 ]
	L05 = eta*( l[1]*eta + l[2]*chi_eff + l[3]*chi_eff*chi_eff + l[0] )

	# Return the luminosity

	return L05*L05

def bbh_aligned_Lpeak22_SH(q, chi1para, chi2para):
	""" 
	Calculate the peak luminosity (currently just due to the 2,2 and 2,-2 modes) in geometrized units of a binary black hole with aligned spins using the fit made by Sascha Husa

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak22_SH)(q, chi1para, chi2para)
	
	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	chi_eff2 = chi_eff*chi_eff
	chi_eff3 = chi_eff2*chi_eff
	chi_eff4 = chi_eff3*chi_eff

	# Evaluate Sascha's fit for N22/sqrt(8pi)
	
	N22norm = 0.378321*eta*eta*eta + 0.721946*(chi1para - chi2para)*(1. - 4.*eta)**0.704535*eta**3.34223 + eta*(0.0870396 + 0.0403944*chi_eff + 0.0404659*chi_eff2 + 0.0180177*chi_eff3 - 0.013243*chi_eff4) + eta*eta*(0.0649694 - 0.0489852*chi_eff - 0.117612*chi_eff2 - 0.0369386*chi_eff3 + 0.0749239*chi_eff4)

	return N22norm*N22norm

def bbh_aligned_Lpeak_6mode_SH(q, chi1para, chi2para):
	""" 
	Calculate the peak luminosity (using modes 22, 21, 33, 32, 44, and 43) in geometrized units of a binary black hole with aligned spins using the fit made by Sascha Husa

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak_6mode_SH)(q, chi1para, chi2para)
	
	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	dm2 = abs(1. - 4.*eta) # We introduce an explicit absolute value to keep from problems with rounding, since we raise this to a fractional power

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	chi_eff2 = chi_eff*chi_eff
	chi_eff3 = chi_eff2*chi_eff
	chi_eff4 = chi_eff3*chi_eff

	chi_diff = chi1para - chi2para

	# Return Sascha's fit

	return 0.000344859*eta3 + 0.0563281*eta4 + 1.27178*chi_diff*dm2**0.805446*eta**5.27679 + 31.1881*chi_diff*chi_diff*dm2**1.45574*eta**7.21542 + eta*(0.000512702*chi_eff + 0.00158118*chi_eff2 + 0.0010312*chi_eff3 - 0.0000894223*chi_eff4) + eta2*(0.0127984 + 0.00516765*chi_eff - 0.00288145*chi_eff2 - 0.000592951*chi_eff3 + 0.00320978*chi_eff4)

def bbh_aligned_Lpeak_6mode_UIB(q, chi1para, chi2para):
	""" 
	Calculate the peak luminosity (using modes 22, 21, 33, 32, 44, and 43) in geometrized units of a binary black hole with aligned spins using the fit made by the UIB group (Sascha Husa, David Keitel, and Xisco Jimenez Forteza)--this is a newer fit than Sascha's original

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak_6mode_UIB)(q, chi1para, chi2para)
	
	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	eta2 = eta*eta
	eta4 = eta2*eta2

	dm2 = 1. - 4.*eta

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	chi_eff2 = chi_eff*chi_eff
	chi_eff3 = chi_eff2*chi_eff
	chi_eff4 = chi_eff3*chi_eff

	chi_diff = chi1para - chi2para

	# Return the UIB fit

	return (0.012851338846828302 + 0.009567880362125462*chi_eff + 0.010592710363652837*chi_eff2 + 0.006764389333842757*chi_eff3 + 0.00047325493557317513*chi_eff4)*eta2 + (0.05681786589129071 - 0.04056527662647619*chi_eff - 0.1107866473486155*chi_eff2 - 0.04525613383811283*chi_eff3 + 0.03108896270210911*chi_eff4)*eta4 + 0.018687217489025302*dm2**0.4687674712184375*eta**3.03713059782352*chi_diff + 0.0005568407558348465*dm2**0.3771581831078265*eta**1.6648159167967396*chi_diff*chi_diff

def bbh_final_mass_and_spin_non_precessing_Husaetal(m1, m2, chi1, chi2): 
	""" 
	Calculate the mass and spin of the final BH resulting from the 
	merger of two black holes with non-precessing spins using the fits
	used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
	arXiv:1508.07250. Note that Eq. (3.8) gives the radiated energy, not
	the final mass directly

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_non_precessing_Husaetal)(m1, m2, chi1, chi2)
	
	# binary parameters 
	m = m1+m2  
        msq = m*m

	eta = m1*m2/msq
	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	S1 = chi1*m1**2/msq	       		# spin angular momentum 1 (in m = 1 units) 
	S2 = chi2*m2**2/msq		       	# spin angular momentum 2 (in m = 1 units)
	S = S1+S2		    	        # total spin
	Sh = S/(1. - 2.*eta)                    # rescaled total spin

	S2 = S*S
	S3 = S2*S
	S4 = S3*S

	# Expressions copied from LALSimIMRPhenomD_internals.c (except with two notation differences: S is capitalized in chif and s -> Sh in Mf, in addition to the "m*(1. - ...)" to obtain the final mass from the radiated mass in m = 1 units which is calculated in the LAL code)

	Mf = m*(1. - ((0.055974469826360077*eta + 0.5809510763115132*eta2 - 0.9606726679372312*eta3 + 3.352411249771192*eta4)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*Sh))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*Sh))

        chif = 3.4641016151377544*eta - 4.399247300629289*eta2 + 9.397292189321194*eta3 - 13.180949901606242*eta4 + (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*S + (0.1014665242971878*eta - 2.0967746996832157*eta2)*S2 + (-1.3546806617824356*eta + 4.108962025369336*eta2)*S3 + (-0.8676969352555539*eta + 2.064046835273906*eta2)*S4 

	return Mf, chif

def bbh_final_mass_and_spin_non_precessing_Meta_Husaetal(m, eta, chi1, chi2): 
	""" 
	Calculate the mass and spin of the final BH resulting from the 
	merger of two black holes with non-precessing spins using the fits
	used by IMRPhenomD, given in Eqs. (3.6) and (3.8) of Husa et al.
	arXiv:1508.07250. Note that Eq. (3.8) gives the radiated energy, not
	the final mass directly

	m: total mass
	eta: symmetric mass ratio
	chi1, chi2: dimensionless spins of two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(eta, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_non_precessing_Meta_Husaetal)(m, eta, chi1, chi2)
	
	# binary parameters   
        msq = m*m
	m1 = (1. + (1. - 4.*eta)**0.5)*m/2.
	m2 = (1. - (1. - 4.*eta)**0.5)*m/2.

	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta

	S1 = chi1*m1**2/msq	       		# spin angular momentum 1 (in m = 1 units) 
	S2 = chi2*m2**2/msq		       	# spin angular momentum 2 (in m = 1 units)
	S = S1+S2		    	        # total spin
	Sh = S/(1. - 2.*eta)                    # rescaled total spin

	S2 = S*S
	S3 = S2*S
	S4 = S3*S

	# Expressions copied from LALSimIMRPhenomD_internals.c (except with two notation differences: S is capitalized in chif and s -> Sh in Mf, in addition to the "m*(1. - ...)" to obtain the final mass from the radiated mass in m = 1 units which is calculated in the LAL code)

	Mf = m*(1. - ((0.055974469826360077*eta + 0.5809510763115132*eta2 - 0.9606726679372312*eta3 + 3.352411249771192*eta4)*
    (1. + (-0.0030302335878845507 - 2.0066110851351073*eta + 7.7050567802399215*eta2)*Sh))/(1. + (-0.6714403054720589 - 1.4756929437702908*eta + 7.304676214885011*eta2)*Sh))

        chif = 3.4641016151377544*eta - 4.399247300629289*eta2 + 9.397292189321194*eta3 - 13.180949901606242*eta4 + (1 - 0.0850917821418767*eta - 5.837029316602263*eta2)*S + (0.1014665242971878*eta - 2.0967746996832157*eta2)*S2 + (-1.3546806617824356*eta + 4.108962025369336*eta2)*S3 + (-0.8676969352555539*eta + 2.064046835273906*eta2)*S4 

	return Mf, chif

def bbh_final_mass_and_spin_precessing_IMRPhenomPv2(m1, m2, chi1z, chi2z, chi_p): 
	""" 
	Calculate the mass and dimensionless spin of the final BH resulting from the 
	merger of two black holes with precessing spins using the modification
	of the IMRPhenomD fit used by IMRPhenomPv2, given in
	FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH in LALSimIMRPhenomP.c
	plus the final mass fit with just the aligned components of the spins,
	again as in IMRPhenomPv2

	m1, m2: component masses
	chi1z, chi2z: components of the dimensionless spins of the two BHs along the orbital angular momentum
	chi_p: IMRPhenomP in-plane spin parameter
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_mass_and_spin_precessing_IMRPhenomPv2)(m1, m2, chi1z, chi2z, chi_p)

	# Compute the ratio of the more massive black hole to the total mass

	qf = max(m1, m2)/(m1 + m2)

	# Compute the component of the spin parallel to the orbital angular momentum using the IMRPhenomD fit

	mf, chif_par = bbh_final_mass_and_spin_non_precessing_Husaetal(m1, m2, chi1z, chi2z)

	# Compute the in-plane spin with the scaling from the mass ratio

	Splane = qf*qf*chi_p

	return mf, (Splane*Splane + chif_par*chif_par)**0.5

def bbh_final_spin_precessing_Barausse_and_Rezzolla(m1, m2, a1, a2, tilt1, tilt2, phi12): 
	""" 
	Calculate the dimensionless spin of the final BH resulting from the 
	merger of two black holes with precessing spins using the fit from Barausse and Rezzolla ApJL 704, L40 (2009). We base our implementation on the IMRPhenomPv2 one in FinalSpinBarausse2009 in LALSimIMRPhenomP.c.

	m1, m2: component masses (with m1 > m2)
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_precessing_Barausse_and_Rezzolla)(m1, m2, a1, a2, tilt1, tilt2, phi12)

	# Coefficients

	s4 = -0.1229
	s5 = 0.4537
	t0 = -2.8904
	t2 = -3.5171
	t3 = 2.5763

	# Computing angles
	cos_beta_tilde = np.cos(tilt1)
	cos_gamma_tilde = np.cos(tilt2)
	cos_alpha = ((1 - cos_beta_tilde*cos_beta_tilde)*(1 - cos_gamma_tilde*cos_gamma_tilde))**0.5*np.cos(phi12) + cos_beta_tilde*cos_gamma_tilde

	# Definitions
	q = m2/m1
	nu = m1*m2/(m1+m2)**2

	# Shorthands
	nu2 = nu*nu
	q2 = q*q
	q4 = q2*q2
	q2p = 1. + q2
	q2p2 = q2p*q2p
	qp = 1. + q
	qp2 = qp*qp
	a1_2 = a1*a1
	a2_2 = a2*a2

	# Compute the final spin and return it

	l = 2.*3.**0.5 + t2*nu + t3*nu2 + (s4 / q2p2) * (a1_2 + a2_2*q4 + 2.*a1*a2*q2*cos_alpha) + ((s5*nu + t0 + 2.)/q2p) * (a1*cos_beta_tilde + a2*cos_gamma_tilde*q2)
	
	return (1. / qp2) * (a1_2 + a2_2*q4 + 2.*a1*a2*q2*cos_alpha + 2.*(a1*cos_beta_tilde + a2*q2*cos_gamma_tilde)*l*q + l*l*q2)**0.5

def bbh_22_peak_frequency_non_precessing_Taracchinietal(m1, m2, chi1, chi2): 
	""" 
	Calculate the peak frequency in Hz of the 2,2 mode of a binary black
	hole gravitational waveform given the individual masses and
	dimensionless (aligned) spins using the fit given in Eq. (42) of
	Taracchini et al. PRD 86, 024011 (2012).

	m1, m2: individual masses of the two BHs (in solar masses) with m1 >= m2
	chi1, chi2: dimensionless spins of the two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( bbh_22_peak_frequency_non_precessing_Taracchinietal)(m1, m2, chi1, chi2)
	
	# Conversion factor Msun -> s

	Msun2s = 4.92549e-6 # s/Msun

	# compute m and eta

	m = m1 + m2
	eta = m1*m2/(m*m)

	# compute chi  

	chiS = 0.5*(chi1 + chi2)
	chiA = 0.5*(chi1 - chi2)
	chi = chiS + chiA*(1 - 4.*eta)**0.5/(1. - 2.*eta) # From Eq. (32) in Taracchini et al.

	# Give fit expressions from Tables IV and V

	fNR0 = 0.2758 - 0.08898*np.log(1. - chi) # The eta = 0 expression for f^NR
	fNR1o4  = 0.3604 + 0.08242*chi + 0.02794*chi*chi # The eta = 1/4 expression for f^NR
	cb1 = 0.1935 # \bar{c}_1 

	# Put everything together as in Eq. (42)

	fNR = (16.*(fNR1o4 - fNR0) - 4.*cb1)*eta*eta + cb1*eta + fNR0

	# Now scale by the mass and convert to Hz

	return fNR/(m*Msun2s*2.*np.pi)

def bbh_22_peak_frequency_non_precessing_Boheetal(m1, m2, chi1, chi2): 
	""" 
	Calculate the peak frequency in Hz of the 2,2 mode of a binary black
	hole gravitational waveform given the individual masses and
	dimensionless (aligned) spins using the fit given in Appendix A 3 of
	Boh{\'e} et al., arXiv:1611.03703.

	m1, m2: individual masses of the two BHs (in solar masses) with m1 >= m2
	chi1, chi2: dimensionless spins of the two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( bbh_22_peak_frequency_non_precessing_Boheetal)(m1, m2, chi1, chi2)
	
	# Conversion factor Msun -> s

	Msun2s = 4.92549e-6 # s/Msun

	# compute m and eta

	m = m1 + m2
	eta = m1*m2/(m*m)

	# compute chi  

	chiS = 0.5*(chi1 + chi2)
	chiA = 0.5*(chi1 - chi2)
	chi = chiS + chiA*(1 - 4.*eta)**0.5/(1. - 2.*eta) # From Eq. (32) in Taracchini et al.; cf. Eq. (2.8) in Boh{\'e} et al. (noting that m1 >= m2)

	# Give coefficients

	p0TPL = 0.562679
	p1TPL = -0.087062
	p2TPL = 0.001743
	p3TPL = 25.850378
	p4TPL = 25.819795

	p3EQ = 10.262073
	p4EQ = 7.629922

	# Put everything together as in Eqs. (A8-A9)

	A3 = p3EQ + 4.*(p3EQ - p3TPL)*(eta - 0.25)
	A4 = p4EQ + 4.*(p4EQ - p4TPL)*(eta - 0.25)

	fNR = p0TPL + (p1TPL + p2TPL*chi)*np.log(A3 - A4*chi)

	# Now scale by the mass and convert to Hz

	return fNR/(m*Msun2s*2.*np.pi)
          
def bbh_22_peak_dimless_frequency_non_precessing_Taracchinietal(m1, m2, chi1, chi2): 
	""" 
	Calculate the peak dimensionless frequency of the 2,2 mode of a binary black
	hole gravitational waveform given the individual masses and
	dimensionless (aligned) spins using the fit given in Eq. (42) of
	Taracchini et al. PRD 86, 024011 (2012).

	m1, m2: individual masses of the two BHs, with m1 >= m2
	chi1, chi2: dimensionless spins of the two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( bbh_22_peak_dimless_frequency_non_precessing_Taracchinietal)(m1, m2, chi1, chi2)
	
	# compute eta

	m = m1 + m2
	eta = m1*m2/(m*m)

	# compute chi  

	chiS = 0.5*(chi1 + chi2)
	chiA = 0.5*(chi1 - chi2)
	chi = chiS + chiA*(1 - 4.*eta)**0.5/(1. - 2.*eta) # From Eq. (32) in Taracchini et al.

	# Give fit expressions from Tables IV and V

	fNR0 = 0.2758 - 0.08898*np.log(1. - chi) # The eta = 0 expression for f^NR
	fNR1o4  = 0.3604 + 0.08242*chi + 0.02794*chi*chi # The eta = 1/4 expression for f^NR
	cb1 = 0.1935 # \bar{c}_1 

	# Put everything together as in Eq. (42)

	return (16.*(fNR1o4 - fNR0) - 4.*cb1)*eta*eta + cb1*eta + fNR0
 
def bbh_22_peak_amplitude_non_precessing_Taracchinietal(m1, m2, chi1, chi2): 
	""" 
	Calculate the peak GW amplitude of the 2,2 mode of a binary black
	hole gravitational waveform given the individual masses and
	dimensionless (aligned) spins using the fit given in Eq. (42) of
	Taracchini et al. PRD 86, 024011 (2012).

	m1, m2: individual masses of the two BHs (in solar masses) with m1 >= m2
	chi1, chi2: dimensionless spins of the two BHs
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( bbh_22_peak_amplitude_non_precessing_Taracchinietal)(m1, m2, chi1, chi2)
	
	# compute eta

	m = m1 + m2
	eta = m1*m2/(m*m)

	# compute chi  

	chiS = 0.5*(chi1 + chi2)
	chiA = 0.5*(chi1 - chi2)
	chi = chiS + chiA*(1 - 4.*eta)**0.5/(1. - 2.*eta) # From Eq. (32) in Taracchini et al.

	# Give fit expressions from Tables IV and V

	fNR0 = 0. # The eta = 0 expression for f^NR
	fNR1o4  = 0.3961 # The eta = 1/4 expression for f^NR
	cb1 = 0.1355 # \bar{c}_1 

	# Put everything together as in Eq. (42)

	return (16.*(fNR1o4 - fNR0) - 4.*cb1)*eta*eta + cb1*eta + fNR0

def peak_freq_bns_BDN_old(m1, m2, lambda1, lambda2):
	""" 
	Calculate the peak frequency in Hz of (the 2,2 mode of) a binary
	neutron star gravitational waveform given the individual masses and
	dimensionless Love numbers using the fit given in Eq. (4) of
	Bernuzzi, Dietrich, and Nagar, arXiv:1504.01764v1. This is superseded
	by the updated fit given in the published PRL, but we keep it, for
	comparison.

	m1, m2: individual masses of the two NSs (in solar masses)
	lambda1, lambda2: dimensionless Love numbers of the two NSs

	Note that we use the definition lamdbda_i = (2/3)k_2^{(i)}/C_i^5, where 
	C_i = M_i/R_i is the compactness of the ith star 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( peak_freq_bns_BDN)(m1, m2, lambda1, lambda2)

	# Conversion factor Msun -> s

	Msun2s = 4.92549e-6 # s/Msun

	# First put everything in the m1 >= m2 convention

	if m1 < m2:
          m1, m2 = m2, m1
	  lambda1, lambda2 = lambda2, lamnbda1

	# Compute kappa_2^T

	q = m1/m2

	kap2T = (4./3.)*q*(q**3.*lambda1 + lambda2)/(1. + q)**5.

	return 0.36*(1. + 2.62e-2*kap2T - 6.32e-6*kap2T*kap2T)/(1. + 6.18e-2*kap2T)/((m1 + m2)*Msun2s)/(2.*np.pi)


def peak_freq_bns_BDN(m1, m2, lambda1, lambda2):
	""" 
	Calculate the peak frequency in Hz of (the 2,2 mode of) a binary
	neutron star gravitational waveform given the individual masses and
	dimensionless Love numbers using the fit given in Eq. (2) and Table II of
	Bernuzzi, Dietrich, and Nagar PRL 115, 091101(2015).

	m1, m2: individual masses of the two NSs (in solar masses)
	lambda1, lambda2: dimensionless Love numbers of the two NSs

	Note that we use the definition lamdbda_i = (2/3)k_2^{(i)}/C_i^5, where 
	C_i = M_i/R_i is the compactness of the ith star 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize( peak_freq_bns_BDN)(m1, m2, lambda1, lambda2)

	# Conversion factor Msun -> s

	Msun2s = 4.92549e-6 # s/Msun

	# First put everything in the m1 >= m2 convention

	if m1 < m2:
          m1, m2 = m2, m1
	  lambda1, lambda2 = lambda2, lamnbda1

	# Compute kappa_2^T

	q = m1/m2

	kap2T = (4./3.)*q*(q**3.*lambda1 + lambda2)/(1. + q)**5.

	return 0.3596*(1. + 2.4384e-2*kap2T - 1.7167e-5*kap2T*kap2T)/(1. + 6.8865e-2*kap2T)/((m1 + m2)*Msun2s)/(2.*np.pi)

def bbh_aligned_Lpeak_6mode_SHXJDK(q, chi1para, chi2para):
	"""
	Calculate the peak luminosity (using modes 22, 21, 33, 32, 44, and 43) in geometrized units of a binary black hole with aligned spins using the fit made by Sascha Husa, Xisco Jimenez Forteza, David Keitel
	using 5th order in chieff

	q: mass ratio (here m1/m2)
	chi1para: the component of the dimensionless spin of m1 along the angular momentum (z)
	chi2para: the component of the dimensionless spin of m2 along the angular momentum (z)
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(q, '__len__'):
		return np.vectorize(bbh_aligned_Lpeak_6mode_SHXJDK)(q, chi1para, chi2para)

	# Calculate eta and the effective spin

	eta = q/(1.+q)**2.

	eta2 = eta*eta
	eta3 = eta2*eta
	eta4 = eta3*eta
	dm2 = 1. - 4.*eta

	chi_eff = (q*chi1para + chi2para)/(1. + q)

	chi_eff2 = chi_eff*chi_eff
	chi_eff3 = chi_eff2*chi_eff
	chi_eff4 = chi_eff3*chi_eff
	chi_eff5 = chi_eff4*chi_eff

	chi_diff = chi1para - chi2para
	chi_diff2 = chi_diff*chi_diff

	# Return best fit (from ../UIBfits/PeakLuminosityUIBCformFit_S5.txt)

        return (0.012851338846828302 + 0.007822265919928252*chi_eff + 0.010221856361035788*chi_eff2 + 0.015805535732661396*chi_eff3 + 0.0011356206806770043*chi_eff4 - 0.009868152529667197*chi_eff5)*eta2 + (0.05681786589129071 - 0.0017473702709303457*chi_eff - 0.10150706091341818*chi_eff2 - 0.2349153289253309*chi_eff3 + 0.015657737820040145*chi_eff4 + 0.19556893194885075*chi_eff5)*eta4 + 0.026161288241420833*dm2**0.541825641769908*eta**3.1629576945611757*chi_diff + 0.0007771032100485481*dm2**0.4499151697918658*eta**1.7800346166040835*chi_diff2

def bbh_final_mass_and_spin_ZL_precessing(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
        """ 
        Calculate the spin of the final BH resulting from the 
        merger of two black holes with precessing spins using
        the Zlochower and Lousto fit from PRD 92, 024022 (2015)

        m1, m2: component masses
        chi1, chi2: dimensionless spins of two BHs
        tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
        phi12: angle between in-plane spin components 
        """
        # Vectorize the function if arrays are provided as input
        if hasattr(m1, '__len__'):
                return np.vectorize(bbh_final_mass_and_spin_ZL_precessing)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

        # Define mass combinations

        m = m1 + m2
        msq = m*m

        eta = m1*m2/msq

        F = 4.*eta
        F2 = F*F

        dM = (m1 - m2)/m

        dM2 = dM*dM
        dM3 = dM2*dM
        dM4 = dM2*dM2
        dM6 = dM4*dM2

        # Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

        S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
        S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

        # Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

        S = (S1 + S2)/msq
        Delta = (S2/m2 - S1/m1)/m
        S0 = S + 0.5*dM*Delta

        # Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

        Sz = S[2]
        Sp2 = S[0]*S[0] + S[1]*S[1]
        Sp = Sp2**0.5
        Ssq = Sp2 + Sz*Sz

        Deltaz = Delta[2]
        Deltaz2 = Deltaz*Deltaz
        Dp2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
        Dp = Dp2**0.5   
        
        S0z = S0[2]
        S0z2 = S0z*S0z
        S0p2 = S0[0]*S0[0] + S0[1]*S0[1]
        S0p = S0p2**0.5

        # Compute ISCO energy and angular momentum corresponding to Sz

        risco = calc_isco_radius(Sz)

        Lisco = (2./3.**0.5)*(2. - 2.*Sz/risco**0.5) # Incorrect version (3. -> 2. in first term of second parentheses) from RIT notebook

        #Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5) # Correct version, but not what this fit seems to be designed for...

        Eisco = (1. - 2./3./risco)**0.5 # cf. the expression above Eq. (2.20) in Bardeen, Press, and Teukolsky ApJ 178, 347 (1972)

        # Compute fractional radiated mass and square of final spin

        aHU = 0.686403 + 0.613203*S0z - 0.107373*S0z2 - 0.0784152*S0z2*S0z - 0.079896*S0z2*S0z2 # Eq. (36) in Zlochower and Lousto

        a2 = F2*(dM2*S0p2*(0.869791 + 4.65308*S0z) + S0p2*(0.8401 - 0.3277*S0z - 0.6088*S0z2) + Dp2*(-0.0209 - 0.0381*S0z + 0.0429*S0z2)) + F2*(-0.00522711*Deltaz2 - 0.556242*Deltaz*dM - 5.55112e-17*dM2 -1.22018*Deltaz*dM3 + 1.14178*dM4 - 6.93889e-18*Deltaz2*S0z - 1.61133*Deltaz*dM*S0z + 1.70469*dM2*S0z + aHU*aHU) + dM6*(Ssq + 12.*Ssq*eta + 2.*Sz*eta*Lisco + F2*(3.*Ssq + Sz*Lisco))

        Eradfrac = F2*(0.0025829 + Deltaz2*(0.000743 + 0.000124*Deltaz2) - 0.0166606*Deltaz*dM - 0.0100117*Deltaz2*Deltaz*dM - 1.73472e-18*dM2 + 0.0370949*dM4 - 0.0187652*dM2*Dp2 - 0.243822*Deltaz2*S0p2 + dM2*S0p2*(-0.0750609 - 0.286384*S0z) + 1.73472e-18*Deltaz2*S0z - 0.0911093*Deltaz*dM*S0z - 0.146479*Deltaz*dM*S0z2 - 0.07730790000000000/(-1.693959000000000 + 2*S0z) + Dp2*(0.0043 + 0.005*S0z - 0.009*S0z2) + S0p2*(0.0356 + 0.096*S0z + 0.1217*S0z2)) + dM6*eta*(1. - Eisco)

        return m*(1. - Eradfrac), np.maximum(a2,0.)**0.5

def bbh_final_mass_and_spin_ZL_precessing_test(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
        """ 
        Calculate the spin of the final BH resulting from the 
        merger of two black holes with precessing spins using
        the Zlochower and Lousto fit from PRD 92, 024022 (2015),
        using the corrected Lisco

        m1, m2: component masses
        chi1, chi2: dimensionless spins of two BHs
        tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
        phi12: angle between in-plane spin components 
        """
        # Vectorize the function if arrays are provided as input
        if hasattr(m1, '__len__'):
                return np.vectorize(bbh_final_mass_and_spin_ZL_precessing_test)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

        # Define mass combinations

        m = m1 + m2
        msq = m*m

        eta = m1*m2/msq

        F = 4.*eta
        F2 = F*F

        dM = (m1 - m2)/m

        dM2 = dM*dM
        dM3 = dM2*dM
        dM4 = dM2*dM2
        dM6 = dM4*dM2

        # Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

        S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
        S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

        # Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

        S = (S1 + S2)/msq
        Delta = (S2/m2 - S1/m1)/m
        S0 = S + 0.5*dM*Delta

        # Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

        Sz = S[2]
        Sp2 = S[0]*S[0] + S[1]*S[1]
        Sp = Sp2**0.5
        Ssq = Sp2 + Sz*Sz

        Deltaz = Delta[2]
        Deltaz2 = Deltaz*Deltaz
        Dp2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
        Dp = Dp2**0.5   
        
        S0z = S0[2]
        S0z2 = S0z*S0z
        S0p2 = S0[0]*S0[0] + S0[1]*S0[1]
        S0p = S0p2**0.5

        # Compute ISCO energy and angular momentum corresponding to Sz

        risco = calc_isco_radius(Sz)

        #Lisco = (2./3.**0.5)*(2. - 2.*Sz/risco**0.5) # Incorrect version (3. -> 2. in first term of second parentheses) from RIT notebook

        Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5) # Correct version, but not what this fit seems to be designed for...

        Eisco = (1. - 2./3./risco)**0.5 # cf. the expression above Eq. (2.20) in Bardeen, Press, and Teukolsky ApJ 178, 347 (1972)

        # Compute fractional radiated mass and square of final spin

        aHU = 0.686403 + 0.613203*S0z - 0.107373*S0z2 - 0.0784152*S0z2*S0z - 0.079896*S0z2*S0z2 # Eq. (36) in Zlochower and Lousto

        a2 = F2*(dM2*S0p2*(0.869791 + 4.65308*S0z) + S0p2*(0.8401 - 0.3277*S0z - 0.6088*S0z2) + Dp2*(-0.0209 - 0.0381*S0z + 0.0429*S0z2)) + F2*(-0.00522711*Deltaz2 - 0.556242*Deltaz*dM - 5.55112e-17*dM2 -1.22018*Deltaz*dM3 + 1.14178*dM4 - 6.93889e-18*Deltaz2*S0z - 1.61133*Deltaz*dM*S0z + 1.70469*dM2*S0z + aHU*aHU) + dM6*(Ssq + 12.*Ssq*eta + 2.*Sz*eta*Lisco + F2*(3.*Ssq + Sz*Lisco))

        Eradfrac = F2*(0.0025829 + Deltaz2*(0.000743 + 0.000124*Deltaz2) - 0.0166606*Deltaz*dM - 0.0100117*Deltaz2*Deltaz*dM - 1.73472e-18*dM2 + 0.0370949*dM4 - 0.0187652*dM2*Dp2 - 0.243822*Deltaz2*S0p2 + dM2*S0p2*(-0.0750609 - 0.286384*S0z) + 1.73472e-18*Deltaz2*S0z - 0.0911093*Deltaz*dM*S0z - 0.146479*Deltaz*dM*S0z2 - 0.07730790000000000/(-1.693959000000000 + 2*S0z) + Dp2*(0.0043 + 0.005*S0z - 0.009*S0z2) + S0p2*(0.0356 + 0.096*S0z + 0.1217*S0z2)) + dM6*eta*(1. - Eisco)

        return m*(1. - Eradfrac), np.maximum(a2,0.)**0.5

def bbh_final_spin_ZL_precessing_new_AC(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Zlochower and Lousto fit from PRD 92, 024022 (2015),
	using the first corrected fit from Yosef

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_ZL_precessing_new_AC)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# Define mass combinations

	m = m1 + m2
	msq = m*m

	eta = m1*m2/msq

	F = 4.*eta
	F2 = F*F

	dM = (m1 - m2)/m

	dM2 = dM*dM
	dM3 = dM2*dM
	dM4 = dM2*dM2
	dM6 = dM4*dM2

	# Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

	S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
	S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

	# Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

	S = (S1 + S2)/msq
	Delta = (S2/m2 - S1/m1)/m
	S0 = S + 0.5*dM*Delta

	# Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

	Sz = S[2]
	Sz2 = Sz*Sz
	Sp2 = S[0]*S[0] + S[1]*S[1]
	Sp = Sp2**0.5

	deltaz = Delta[2]
	deltaz2 = deltaz*deltaz
	deltap2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
	deltap = deltap2**0.5	

	# Compute ISCO energy and angular momentum corresponding to Sz

	risco = calc_isco_radius(Sz)

	Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5)

	aalign = (0.68671 - 0.005254*deltaz2 + 0.004364*deltaz2*deltaz2 - 0.145427*deltaz*dM - 0.006018*deltaz2*deltaz*dM + 0.801838*dM2 + 0.001629*deltaz2*dM2 - 0.067998*deltaz*dM3 + 0.953458*dM4 + 0.613247*Sz + 0.004759*deltaz2*Sz - 0.073839*deltaz*dM*Sz + 1.58581*dM2*Sz - 0.115689*Sz2 - 0.053099*deltaz2*Sz2 - 0.066693*dM2*Sz2 - 0.078377*Sz2*Sz - 0.047204*Sz2*Sz2)*F2 + dM4*Sz*(1 + 2.*F) + dM6*eta*Lisco # \alpha_\text{align}

	Ac = F2*(-2.80392*deltap2*deltaz2 + 5.13991*deltaz2*Sp2 + dM2*Sp2*(3.00684 - 7.32272*Sz) + deltap2*dM2*(-2.00948 + 5.08985*Sz) + Sp2*(0.8401 - 0.3277*Sz - 0.6088*Sz2) + deltap2*(-0.0209 - 0.0381*Sz + 0.0429*Sz2)) + dM6*(Sp2 + 12.*Sp2*eta) + aalign*aalign

	return Ac**0.5 #np.minimum(Ac**0.5, 1.)

def bbh_final_spin_ZL_precessing_new_APC(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Zlochower and Lousto fit from PRD 92, 024022 (2015),
	using the second corrected fit from Yosef

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_ZL_precessing_new_APC)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# Define mass combinations

	m = m1 + m2
	msq = m*m

	eta = m1*m2/msq

	F = 4.*eta
	F2 = F*F

	dM = (m1 - m2)/m

	dM2 = dM*dM
	dM3 = dM2*dM
	dM4 = dM2*dM2
	dM6 = dM4*dM2

	# Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

	S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
	S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

	# Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

	S = (S1 + S2)/msq
	Delta = (S2/m2 - S1/m1)/m
	S0 = S + 0.5*dM*Delta

	# Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

	Sz = S[2]
	S0z = S0[2]
	S0z2 = S0z*S0z
	Sp2 = S[0]*S[0] + S[1]*S[1]
	Sp = Sp2**0.5

	deltaz = Delta[2]
	deltaz2 = deltaz*deltaz
	deltap2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
	deltap = deltap2**0.5	

	S0p = Sp + 0.5*dM*deltap
	S0p2 = S0p*S0p

	# Compute ISCO energy and angular momentum corresponding to Sz

	risco = calc_isco_radius(Sz)

	Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5)

	aalign = F2*(0.686403 - 0.51655*deltaz*dM + 0.809512*dM2 + 0.691714*dM4 + 1.00321*deltaz*dM4 + 0.613203*S0z + 1.58707*dM2*S0z - 0.107373*S0z2 - 0.0784152*S0z2*S0z - 0.079896*S0z2*S0z2) + dM4*(1. + 2.*F)*Sz + dM6*eta*Lisco # \alpha_\text{align}

	Ac = F2*(0.547173*deltap2*dM2 + 2.18869*dM2*S0p2 + S0p2*(0.8401 - 0.3277*S0z - 0.6088*S0z2) + deltap2*(-0.0209 - 0.0381*S0z + 0.0429*S0z2)) + dM6*(Sp2 + 12.*Sp2*eta) + aalign*aalign

	return Ac**0.5 #np.minimum(Ac**0.5, 1.)

def bbh_final_spin_ZL_precessing_new_APCnew(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Zlochower and Lousto fit from PRD 92, 024022 (2015),
	using the second corrected fit from Yosef

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_ZL_precessing_new_APCnew)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# Define mass combinations

	m = m1 + m2
	msq = m*m

	eta = m1*m2/msq

	F = 4.*eta
	F2 = F*F

	dM = (m1 - m2)/m

	dM2 = dM*dM
	dM3 = dM2*dM
	dM4 = dM2*dM2
	dM6 = dM4*dM2

	# Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

	S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
	S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

	# Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

	S = (S1 + S2)/msq
	Delta = (S2/m2 - S1/m1)/m
	S0 = S + 0.5*dM*Delta

	# Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

	Sz = S[2]
	S0z = S0[2]
	S0z2 = S0z*S0z
	Sp2 = S[0]*S[0] + S[1]*S[1]
	Sp = Sp2**0.5

	deltaz = Delta[2]
	deltaz2 = deltaz*deltaz
	deltap2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
	deltap = deltap2**0.5	

	# Compute ISCO energy and angular momentum corresponding to Sz

	risco = calc_isco_radius(Sz)

	Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5)

	aalign = F2*(0.686403 - 0.51655*deltaz*dM + 0.809512*dM2 + 0.691714*dM4 + 1.00321*deltaz*dM4 + 0.613203*S0z + 1.58707*dM2*S0z - 0.107373*S0z2 - 0.0784152*S0z2*S0z - 0.079896*S0z2*S0z2) + dM4*(1. + 2.*F)*Sz + dM6*eta*Lisco # \alpha_\text{align}

	Ac = F2*(-2.55016*deltap2*dM2 + deltap2*(-0.0209 - 0.0381*S0z + 0.0429*S0z2) + 4.33119*dM2*Sp2 + (0.8401 - 0.3277*S0z - 0.6088*S0z2)*Sp2)  + dM6*(Sp2 + 12.*Sp2*eta) + aalign*aalign

	return Ac**0.5 #np.minimum(Ac**0.5, 1.)

def bbh_final_spin_ZL_precessing_new_ALT(m1, m2, chi1, chi2, tilt1, tilt2, phi12):
	""" 
	Calculate the spin of the final BH resulting from the 
	merger of two black holes with precessing spins using
	the Zlochower and Lousto fit from PRD 92, 024022 (2015),
	using the third corrected fit from Yosef

	m1, m2: component masses
	chi1, chi2: dimensionless spins of two BHs
	tilt1, tilt2: tilt angles of the spins from the orbital angular momentum
	phi12: angle between in-plane spin components 
	"""
	# Vectorize the function if arrays are provided as input
	if hasattr(m1, '__len__'):
		return np.vectorize(bbh_final_spin_ZL_precessing_new_ALT)(m1, m2, chi1, chi2, tilt1, tilt2, phi12)

	# Define mass combinations

	m = m1 + m2
	msq = m*m

	eta = m1*m2/msq

	F = 4.*eta
	F2 = F*F

	dM = (m1 - m2)/m

	dM2 = dM*dM
	dM3 = dM2*dM
	dM4 = dM2*dM2
	dM6 = dM4*dM2

	# Compute the dimensionful spin vectors, aligning the in-plane portion of S1 along the x-axis, without loss of generality

	S1 = m1*m1*chi1*np.array([np.sin(tilt1), 0., np.cos(tilt1)])
	S2 = m2*m2*chi2*np.array([np.sin(tilt2)*np.cos(phi12), np.sin(tilt2)*np.sin(phi12), np.cos(tilt2)])

	# Compute the vector dimensionless spin quantities for the binary [Eqs. (6)--(8) in Zlochower and Lousto]; note that these are the tilded quantities, though we do not note this in the name

	S = (S1 + S2)/msq
	Delta = (S2/m2 - S1/m1)/m
	S0 = S + 0.5*dM*Delta

	# Select out parallel and magnitude of perpendicular components of these spin vectors; we follow the notation in the RIT group's Mathematica notebook

	Sz = S[2]
	Sz2 = Sz*Sz
	Sp2 = S[0]*S[0] + S[1]*S[1]
	Sp = Sp2**0.5

	deltaz = Delta[2]
	deltaz2 = deltaz*deltaz
	deltap2 = Delta[0]*Delta[0] + Delta[1]*Delta[1]
	deltap = deltap2**0.5

	# Compute ISCO energy and angular momentum corresponding to Sz

	risco = calc_isco_radius(Sz)

	Lisco = (2./3.**0.5)*(3. - 2.*Sz/risco**0.5)

	aalign = dM4*(1 + 2.*F)*Sz + F2*(0.686403 - 0.00564042*deltaz2 - 0.148208*deltaz*dM - 0.0136799*deltaz2*deltaz*dM + 0.815626*dM2 + 0.836887*dM4 + 0.613203*Sz - 0.0946628*deltaz*dM*Sz + 1.61408*dM2*Sz - 0.107373*Sz2 - 0.0845327*deltaz2*Sz2 - 0.0914532*dM2*Sz2 - 0.0784152*Sz2*Sz - 0.079896*Sz2*Sz2) + dM6*eta*Lisco # \alpha_\text{align}

	Ac = F2*(-3.09868*deltap2*deltaz2 + 5.99933*deltaz2*Sp2 + dM2*Sp2*(3.40195 - 4.57129*Sz) + deltap2*dM2*(-2.21562 + 3.71341*Sz) + Sp2*(0.8401 - 0.3277*Sz - 0.6088*Sz2) + deltap2*(-0.0209 - 0.0381*Sz + 0.0429*Sz2)) + dM6*(Sp2 + 12.*Sp2*eta) + aalign*aalign

	return Ac**0.5 #np.minimum(Ac**0.5, 1.)
