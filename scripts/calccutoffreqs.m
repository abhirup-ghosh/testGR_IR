% calculate the cutoff frequencies for the inspiral-merger-ringdown consistency test 
% 
% P. Ajith, 2015-09-17 

clear 
addpath('../src');
setconstants 

% set binary parameters here -- for the time being we neglect the spin paraemter 
% use the final-mass/final-spin formula for non-spinning binaries 
m1 = 80;
m2 = 20;
chi = -0.05;

% calc. total mass, mass ratio, etc 
m = m1+m2;
q = m1/m2;
eta = q/(1+q)^2; 

% calculate the mass and spin of the final BH 
[mf, af, A, Omega] = finalmassandspin_eobnrv2(m, eta); 

mf
af

% QNM freq (l = 2, m = 2, n = 0 mode) - Berti's formula 
f_QNM = real(Omega(1))/(2*pi*m*LAL_MTSUN_SI)
tau = -1/imag(Omega(1))*m*LAL_MTSUN_SI;
Q = tau*f_QNM*sqrt(2)*pi;
sigma = f_QNM/Q;

% ISCO freq - Schwarzschild 
[flso_Schw, flRing_Schw] = calcflso(m);

% ISCO freq - Kerr. Also QNM freq using old formulas -- for cross checking 
[flso_Kerr] = CalcFlsoKerrBH(mf, af)

