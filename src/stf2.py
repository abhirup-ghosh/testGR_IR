import numpy as np 
import matplotlib.pyplot as plt 
import plotsettings
from numpy.core.umath_tests import inner1d
#from lal import MTSUN_SI, PC_SI, C_SI
from numpy import sin, cos, arccos, arctan2, log 
#from lal import PI as pi 
GAMMA = 0.577216
pi = np.pi
MTSUN_SI = 1
PC_SI = 1
C_SI = 1

""" compute the orbital phase of a precessing binary in timedomain in spintaylort5 approximation """
def compute_orbital_phase(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, phi_0, phaseO):

	# Define powers of v to make writing easier(change powers to multiplication)		
	v2 = v*v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v; v7 = v6*v

	# Define different phiase coefficients 
	phiN = phi2 = phi3 = phi4 = phi5 = phi6 = phi7 = 0.

	if phaseO >= 0:
		phiN = 3./(256.*eta)
	if phaseO >= 2:
		phi2 = (4.914021164021164 + (55*eta)/9.)
	if phaseO >= 3:
		phi3 = ((113*chisdL)/3. + (113*chiadL*delta)/3. - (76*chisdL*eta)/3. - 16*np.pi)
	if phaseO >= 4:
		phi4 = (30.103152950995213 - (3595*(chiadL**2))/48. + (1165*chiasqr)/48. - (3595*(chisdL**2))/48. + \
			(1165*chissqr)/48. + (1165*chisdchia*delta)/24. - (3595*chiadL*chisdL*delta)/24. + (27145*eta)/504. +\
			 300*(chiadL**2)*eta - 100*chiasqr*eta - (5*(chisdL**2)*eta)/12. + (35*chissqr*eta)/12. +   (3085*(eta**2))/72.)
	if phaseO >= 5:
		phi5 =	((-5*chiadL*delta*(146597 + 7056*eta))/2268. - (5*(chisdL*(146597 - 135856*eta - 17136*(eta**2)) +\
			 3*(-7729 + 1092*eta)*np.pi))/2268.)*(3*np.log(v) + 1)		
	if phaseO >= 6:
		phi6 = 2467.5541189728633 - (15737765635*eta)/3.048192e6 + (76055*(eta**2))/1728. - (127825*(eta**3))/1296.\
			- (640*(np.pi**2))/3. + (2255*eta*(np.pi**2))/12. -  (6848*GAMMA)/21. - (6848*np.log(4))/21. - (6848*np.log(v))/21. 
	if phaseO >= 7:
		phi7 = (77096675*np.pi)/254016. + (378515*eta*np.pi)/1512. - (74045*(eta**2)*np.pi)/756.

	phi_orb = (phiN/v5)*(1 + phi2*v2 + phi3*v3 + phi4*v4 +  phi5*v5 + phi6*v6 + phi7*v7) 

	# add a constant phase such that phi_orb at v=v0 is phi_0/2. 
	return phi_orb-phi_orb[0]+phi_0/2.


def calc_modulation_functions(alpha, angle_i, incl_angle): 

	# compute the modulation functions 
	Cp = 2*cos(alpha)**2*(1 + cos(angle_i)**2*cos(incl_angle)**2) - (2 + cos(2*angle_i) + cos(2*incl_angle))*sin(alpha)**2 \
		+ 2*sin(angle_i)**2*sin(incl_angle)**2 + cos(alpha)*sin(2*angle_i)*sin(2*incl_angle)
	Sp = 2*sin(alpha)*(-(cos(alpha)*cos(angle_i)*(3 + cos(2*incl_angle))) - 2*cos(incl_angle)*sin(angle_i)*sin(incl_angle))
	Cc = sin(alpha)*(2*cos(alpha)*(3 + cos(2*angle_i))*cos(incl_angle) + 4*cos(angle_i)*sin(angle_i)*sin(incl_angle))
	Sc = 4*(cos(2*alpha)*cos(angle_i)*cos(incl_angle) + cos(alpha)*sin(angle_i)*sin(incl_angle))

	return Cp, Sp, Cc, Sc 

""" compute the GW polarizations in the radiation frame """
def calc_polarizations_in_rad_frame(m, eta, d, incl_angle, v, PhiS, alpha, angle_i): 

	Cp, Sp, Cc, Sc = calc_modulation_functions(alpha, angle_i, incl_angle)

	# compute the amplitude and phase of the polarizations 
	d_sec = d*PC_SI/C_SI
	Ap = -m*MTSUN_SI*v**2.*eta*np.sqrt(Cp**2.+Sp**2.)/d_sec
	Ac = -m*MTSUN_SI*v**2.*eta*np.sqrt(Cc**2.+Sc**2.)/d_sec
	Psip = np.unwrap(2*PhiS-arctan2(Sp, Cp))
	Psic = np.unwrap(2*PhiS-arctan2(Sc, Cc))

	# return the ampl and phaseo of the polzns in radiation frame
	return Ap, Ac, Psip, Psic 

""" compute the re-exapnded dEnergy/flux """
def denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, order):

	# different powers of v 
	v2 = v*v; v3 = v2*v; v4 = v3*v; v5 = v4*v; v6 = v5*v; v7 = v6*v; v9 = v7*v2 

	# initialize the cofficients
	dEbF0 = dEbF2 = dEbF3 = dEbF4 = dEbF5 = dEbF6 = dEbF6L = dEbF7 = 0.

	if order >= 0:
		dEbF0 = -5./(32.*eta)
	if order >= 2:
		dEbF2 = 2.2113095238095237 + (11*eta)/4. 
	if order >= 3:
		dEbF3 = (113*chiadL*delta)/12. + chisdL*(9.416666666666666 - (19.*eta)/3.) - 4.*pi 
	if order >= 4:
		dEbF4 = 3.010315295099521 + (233*chisdchia*delta)/48. - (719.*chiadL*chisdL*delta)/48. + \
			chiasqr*(2.4270833333333335 - 10.*eta) + pow(chisdL,2.)*(-7.489583333333333 - eta/24.) + \
			chissqr*(2.4270833333333335 + (7.*eta)/24.) + (5429.*eta)/1008. + (617*pow(eta,2.))/144. + \
			pow(chiadL,2.)*(-7.489583333333333 + 30.*eta) 
	if order >= 5:
			dEbF5 = chiadL*delta*(72.71676587301587 + (7*eta)/2.) + chisdL*(72.71676587301587 - \
						(1213*eta)/18. - (17*pow(eta,2))/2.) - (7729*pi)/672. + (13*eta*pi)/8. 
	if order >= 6:
		dEbF6 = -115.2253249962622 - (15211*pow(eta,2))/6912. + (25565*pow(eta,3))/5184. + \
					(32*pow(pi,2))/3. + eta*(258.1491854023631 - (451*pow(pi,2))/48.) + (1712*GAMMA)/105. 
		dEbF6L = 1712./105. 
	if order >= 7:
		dEbF7 = (-15419335.*pi)/1.016064e6 - (75703.*eta*pi)/6048. + (14809.*pow(eta,2)*pi)/3024. 

	return (dEbF0/v9)*(1. + dEbF2*v2 + dEbF3*v3 + dEbF4*v4 + dEbF5*v5 + (dEbF6+dEbF6L*log(4.*v))*v6 + dEbF7*v7) 


""" compute the spintaylorf2 approximant """
def spintaylorf2(m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, Lx, Ly, Lz, v, incl_angle, d, dt, phi_0, phaseO, amplO, mkplots):

	# compute total mass and mass ratios 
	q = m2/m1 
	m = m1+m2 
	eta = m1*m2/m**2.
	delta = (m1-m2)/m

	# these quantities show the time evolution of the spin and L vectors 
	L = np.transpose([Lx,Ly,Lz])		
	chi1 = np.transpose([chi1x,chi1y,chi1z])
	chi2 = np.transpose([chi2x,chi2y,chi2z])
	chia = (chi1-chi2)/2.	
	chis = (chi1+chi2)/2.	

	# in the case of non-precessing binaries (Lx ~ Ly ~ 0) the alpha estimation is corrupted 
	# by numerical errors. this is  a temporary hack to avoid that. FIXME 
	zero_idx = abs(L[:,1]) < 1e-13 
	L[zero_idx,1] = 0.
	zero_idx = abs(L[:,0]) < 1e-13 
	L[zero_idx,0] = 0.

	# angles describing the evolution of LN
	alpha = np.unwrap(np.arctan2(L[:,1],L[:,0]))
	angle_i = np.unwrap(np.arccos(L[:,2]))
	
	# compute the inner products between spins and L 
	chisdL = inner1d(chis, L)          	# chis.L
	chiadL = inner1d(chia, L)          	# chia.L
	chissqr = inner1d(chis, chis) 		# chis.chis 	
	chiasqr = inner1d(chia, chia) 		# chia.chia
	chisdchia = inner1d(chis, chia)     # chis.chia 

	# compute the orbital phase 
 	phi_orb = compute_orbital_phase(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, phi_0, phaseO)

	# compute the phase at ascending nodes 
	dN = 1
	dalphadt = np.gradient(alpha, dN)/(dN*dt)  
	dalphadt_cosi = dalphadt*cos(angle_i)
	PhiS = phi_orb-np.cumsum(dalphadt_cosi)*dt 

	# compute the ampl and phase of the polarizations in time domain  
	Ap_t, Ac_t, Psip_t, Psic_t = calc_polarizations_in_rad_frame(m, eta, d, incl_angle, v, PhiS, alpha, angle_i)

	# compute the modulation observed in the radiation frame 
	Cp, Sp, Cc, Sc = calc_modulation_functions(alpha, angle_i, incl_angle)
	xi_p = 2.*dalphadt_cosi+np.gradient(np.unwrap(arctan2(Sp, Cp)), dN)/(dN*dt)
	xi_c = 2.*dalphadt_cosi+np.gradient(np.unwrap(arctan2(Sc, Cc)), dN)/(dN*dt)
	
	d_arctan_SpbyCp = np.gradient(np.unwrap(arctan2(Sp, Cp)), dN)/(dN*dt)
	d_arctan_ScbyCc = np.gradient(np.unwrap(arctan2(Sc, Cc)), dN)/(dN*dt)

	# plot the corrected stationary point againt v
	f = v**3./(pi*m*MTSUN_SI)
	vf_p = v*(1+xi_p/(2*pi*f))**(1./3)
	vf_c = v*(1+xi_c/(2*pi*f))**(1./3)

	f1_p = vf_p**3./(pi*m*MTSUN_SI)
	f1_c = vf_c**3./(pi*m*MTSUN_SI)

	dv = np.gradient(v,dN) 
	dxibydv_p = np.gradient(xi_p,dN)/dv
	dxibydv_c = np.gradient(xi_c,dN)/dv
	dxibydv_p *= m*MTSUN_SI/(12.*v**2)
	dxibydv_c *= m*MTSUN_SI/(12.*v**2)

	plt.figure()
	plt.plot(v, dxibydv_p, 'r', label='$d\\xi_+/dv$')
	plt.plot(v, dxibydv_c, 'b', label='$d\\xi_\\times/dv$')
	plt.xlabel('$v$')
	plt.ylabel('d\\xi/dv')
	plt.legend()
	plt.grid()
	plt.show()
	
	# compute dE/flux 
	dEbyF = denergy_by_flux(v, eta, delta, chiadL, chisdL, chiasqr, chissqr, chisdchia, amplO)

	Afp = -Ap_t*np.sqrt(-dEbyF*pi/3.)*m*MTSUN_SI/(2.*v)
	Afc = -Ac_t*np.sqrt(-dEbyF*pi/3.)*m*MTSUN_SI/(2.*v)

	# compute the frequency domain phase 
	#int_xi_p = np.cumsum(xi_p)*dt
	#int_xi_c = np.cumsum(xi_c)*dt
	
	int_xi_p = 2*np.cumsum(dalphadt_cosi)*dt-np.unwrap(arctan2(Sp, Cp))
	int_xi_c = 2*np.cumsum(dalphadt_cosi)*dt-np.unwrap(arctan2(Sc, Cc))

	plt.figure(figsize=(10,4))
	plt.subplot(231)
	plt.plot(v, alpha, 'c', label='$\\alpha$')
	plt.plot(v, angle_i, 'm', label='$i$')
	plt.plot(v, Sp, 'r', label='$S_+$')
	plt.plot(v, Cp, 'b', label='$C_+$')
	plt.plot(v, Sc, 'r--', lw=1, label='$S_\\times$')
	plt.plot(v, Cc, 'b--', lw=1, label='$C_\\times$')
	plt.legend(frameon=False)
	plt.grid()
	plt.subplot(232)
	plt.plot(v, L[:,1], 'b', label='$L_y$')
	plt.plot(v, L[:,0], 'r', label='$L_x$')
	plt.legend(frameon=False)
	plt.grid()
	plt.subplot(233)
	plt.plot(v, int_xi_p, 'b', lw=2, label='$\int \\xi_+ dt $') 
	plt.plot(v, 2*np.cumsum(dalphadt_cosi)*dt, 'c', lw=1, label='$2 \int d\\alpha/dt \cos i dt $') 
	plt.plot(v, -np.unwrap(arctan2(Sp, Cp)), 'm', lw=1, label='$\\arctan(S_+/C_+)$') 
	plt.grid()
	plt.legend(frameon=False)
	plt.xlabel('$v$')
	plt.subplot(234)
	plt.plot(v, int_xi_c, 'b', lw=2, label='$\int \\xi_\\times dt $') 
	plt.plot(v, 2*np.cumsum(dalphadt_cosi)*dt, 'c', lw=1,  label='$2 \int d\\alpha/dt \cos i dt $') 
	plt.plot(v, -np.unwrap(arctan2(Sc, Cc)), 'm', lw=1, label='$\\arctan(S_\\times/C_\\times)$') 
	plt.grid()
	plt.legend(frameon=False)
	plt.xlabel('$v$')
	plt.subplot(235)
	plt.semilogx(f, xi_p/(2*pi*f), 'r', f, xi_c/(2*pi*f), 'b')
	plt.grid()
	plt.xlim(30, 220)
	plt.xlabel('$f$ [Hz]')
	plt.ylabel('$\\xi_+/2\pi f$ and $\\xi_\\times/2\pi f$')
	plt.subplot(236)
	plt.semilogx(f, d_arctan_SpbyCp, 'r', f, d_arctan_ScbyCc, 'b')
	plt.grid()
	plt.xlim(30, 220)
	plt.xlabel('$f$ [Hz]')
	plt.ylabel('$d \\arctan(S/C) / dt$')
	plt.show()

	int_xi_p -= int_xi_p[0]
	int_xi_c -= int_xi_c[0]

	# plot the evolution of various quantities 
	if mkplots == True: 




		chi1sqr = inner1d(chi1, chi1) 		# chi1.chi1 	
		chi2sqr = inner1d(chi2, chi2) 		# chi2.chi2 	

		plt.figure(figsize=(6,6))
		plt.subplot(221)
		plt.plot(v, chi1[:,0], label='$\chi_{1x}$')
		plt.plot(v, chi1[:,1], label='$\chi_{1x}$', lw=1)
		plt.plot(v, chi1[:,2], label='$\chi_{1z}$', lw=1)
		plt.plot(v, chi2[:,0], label='$\chi_{2x}$')
		plt.plot(v, chi2[:,1], label='$\chi_{2x}$', lw=1)
		plt.plot(v, chi2[:,2], label='$\chi_{2z}$', lw=1)
		plt.plot(v, chi1sqr, '--', label='$\chi_1^2$', lw=1)
		plt.plot(v, chi2sqr, '--', label='$\chi_2^2$', lw=1)
		plt.xlabel('$v$')
		plt.ylabel('spins')
		plt.legend(loc=4, frameon=False)
		plt.tight_layout()
		plt.subplot(222)
		plt.plot(v, L[:,0], label='$L_x$')
		plt.plot(v, L[:,1], label='$L_x$')
		plt.plot(v, L[:,2], label='$L_z$')
		plt.xlabel('$v$')
		plt.ylabel('$L_N$')
		plt.legend(loc=4, frameon=False)
		plt.tight_layout()
		plt.subplot(223)
		plt.plot(v, alpha, label='$\\alpha$')
		plt.plot(v, angle_i, label='$i$')
		plt.xlabel('$v$')
		plt.ylabel('$L$ angles')
		plt.legend(loc=4, frameon=False)
		plt.tight_layout()
		plt.subplot(224)
		plt.plot(v, chisdL, label='$\chi_s \cdot L$')
		plt.plot(v, chiadL, label='$\chi_a \cdot L$')
		plt.plot(v, chisdchia, label='$\chi_s \cdot \chi_a$')
		plt.plot(v, chissqr, label='$|\chi_s|^2$')
		plt.plot(v, chiasqr, label='$|\chi_a|^2$')
		plt.xlabel('$v$')
		plt.ylabel('dot products')
		plt.legend(loc=4, frameon=False)
		plt.tight_layout()
		plt.show()

		plt.figure()
		plt.subplot(121)
		plt.plot(v, xi_p, 'b', v, xi_c, 'r')
		plt.xlabel('$v$')
		plt.ylabel('$\\xi_{+,\\times}$')
		plt.subplot(122)
		plt.plot(v, dxibydv_p, 'b', v, dxibydv_c, 'r')
		plt.xlabel('$v$')
		plt.ylabel('$d\\xi_{+,\\times}/dv/(12 v^2)$')
		plt.show()
	
	return Afp, Afc, Psip_t, Psic_t, phi_orb, PhiS, dEbyF, dxibydv_p, dxibydv_c, int_xi_p, int_xi_c 
