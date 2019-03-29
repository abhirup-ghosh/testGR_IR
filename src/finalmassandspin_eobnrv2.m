function [mf, af, A, Omega] = finalmassandspin_eobnrv2(m, eta)
% 
% FINALMASSANDSPIN_EOBNRV2 - calculate the final mass and spin of a merged BH binary 
% using EOBNRv2 fitting formula. Initial BHs are assumed to be non-spinning. 
% 
% usage: [mf, af, A, Omega] = finalmassandspin_eobnrv2(m, eta)
% 
% m, eta : total mass and symmetric mass ratio of the binary 
% mf, af : mass and spin of the final BH 
% A      : complex amplitude of the first seven overtones (n=0..7) of the l=m=2 mode  A = A_abs exp(i phi)
% Omega  : complex frequency of the first seven overtones (n=0..7) of the l=m=2 mode. Omega = omega - i/tau
% 
% Note: Amplitude and frequency should be rescaled with the total mass of the binary; not the final mass. 
% 
% P. Ajith, 15 March 2012 
% 
% $Id: finalmassandspin_eobnrv2.m 3161 2015-05-05 02:24:38Z ajith $

% final mass and final spin fits [Pan et al 1106.1021; Eq. (29)]
mf = m*(1 + ((8/9.)^0.5 - 1.)*eta - 0.4333*eta.^2 - 0.4392*eta.^3); 
af = 12^0.5*eta - 3.871*eta.^2 + 4.028*eta.^3;		% This actually (af/mf) note that this is in units of mf

% fits to the amplityde  of the first three overtones of the l=m=2 
% amplitude of the ring down (n = 0...7). (fits from data provided
% by the EOBNRv2 code) 
A_abs(1) = 1.9597e-02 + 1.7997e+00*eta^1 + 2.8446e+01*eta^2 + -3.4359e+02*eta^3 + 2.3337e+03*eta^4 + -8.4717e+03*eta^5 + 1.3557e+04*eta^6;
A_abs(2) = 4.9074e-02 + 1.4406e+00*eta^1 + 7.4339e+01*eta^2 + -8.2738e+02*eta^3 + 5.5826e+03*eta^4 + -1.9023e+04*eta^5 + 3.2774e+04*eta^6;
A_abs(3) = 1.3187e-01 + -4.2678e+00*eta^1 + 1.8816e+02*eta^2 + -2.1684e+03*eta^3 + 1.4310e+04*eta^4 + -4.6620e+04*eta^5 + 7.1710e+04*eta^6;
A_abs(4) = 1.6492e-01 + -7.8883e+00*eta^1 + 2.3387e+02*eta^2 + -2.7713e+03*eta^3 + 1.8468e+04*eta^4 + -6.1365e+04*eta^5 + 9.1920e+04*eta^6;
A_abs(5) = 1.5475e-01 + -8.3554e+00*eta^1 + 2.1225e+02*eta^2 + -2.5231e+03*eta^3 + 1.6649e+04*eta^4 + -5.5755e+04*eta^5 + 8.1116e+04*eta^6;
A_abs(6) = 8.2481e-02 + -4.5976e+00*eta^1 + 1.1026e+02*eta^2 + -1.3031e+03*eta^3 + 8.4833e+03*eta^4 + -2.8284e+04*eta^5 + 4.0056e+04*eta^6;
A_abs(7) = 1.4592e-02 + -8.1989e-01*eta^1 + 1.9532e+01*eta^2 + -2.3103e+02*eta^3 + 1.5028e+03*eta^4 + -5.0100e+03*eta^5 + 7.0623e+03*eta^6;
A_abs(8) = 8.6923e-04 + -4.8959e-02*eta^1 + 1.1768e+00*eta^2 + -1.4013e+01*eta^3 + 9.1689e+01*eta^4 + -3.0703e+02*eta^5 + 4.3445e+02*eta^6;

A_arg(1) = -2.4858e+00 + 2.5396e+00*eta^1 + -1.4878e+00*eta^2 + 2.9794e+01*eta^3;
A_arg(2) = 1.4616e+00 + 1.2826e-01*eta^1 + 3.6362e+01*eta^2 + -8.2587e+01*eta^3;
A_arg(3) = -1.4167e+00 + 1.9941e+00*eta^1 + 4.0849e+01*eta^2 + -1.2077e+02*eta^3;
A_arg(4) = 1.7823e+00 + 5.6680e+00*eta^1 + 2.5379e+01*eta^2 + -1.0750e+02*eta^3;
A_arg(5) = -1.3492e+00 + 8.9583e+00*eta^1 + 6.8264e+00*eta^2 + -8.1815e+01*eta^3;
A_arg(6) = 1.8427e+00 + 1.0453e+01*eta^1 + -4.3855e+00*eta^2 + -6.8326e+01*eta^3;
A_arg(7) = -1.1692e+00 + 9.1804e+00*eta^1 + -4.4233e-01*eta^2 + -8.2307e+01*eta^3;
A_arg(8) = 2.0957e+00 + 6.4560e+00*eta^1 + 1.3641e+01*eta^2 + -1.0925e+02*eta^3;

omega(1) = 3.8153e-01 + 2.7783e-01*eta^1 + 1.6219e+00*eta^2;
omega(2) = 3.5493e-01 + 3.3120e-01*eta^1 + 1.6343e+00*eta^2;
omega(3) = 3.0970e-01 + 4.1750e-01*eta^1 + 1.6434e+00*eta^2;
omega(4) = 2.5959e-01 + 4.9876e-01*eta^1 + 1.6314e+00*eta^2;
omega(5) = 2.1274e-01 + 5.7745e-01*eta^1 + 1.5523e+00*eta^2;
omega(6) = 1.7057e-01 + 6.5914e-01*eta^1 + 1.6131e+00*eta^2;
omega(7) = 1.3012e-01 + 8.3389e-01*eta^1 + 1.5422e+00*eta^2;
omega(8) = 9.8427e-02 + 1.0591e+00*eta^1 + 1.1535e+00*eta^2;

tau(1) = 1.1107e+01 + 1.9933e+00*eta^1;
tau(2) = 3.6092e+00 + 9.1964e-01*eta^1;
tau(3) = 2.0697e+00 + 8.2133e-01*eta^1;
tau(4) = 1.4036e+00 + 8.0517e-01*eta^1;
tau(5) = 1.0415e+00 + 7.9676e-01*eta^1;
tau(6) = 8.2230e-01 + 7.6581e-01*eta^1;
tau(7) = 6.8404e-01 + 6.5816e-01*eta^1;
tau(8) = 5.8920e-01 + 5.3011e-01*eta^1;

A = A_abs.*exp(1i*A_arg);		% complex amplitude: A = A_abs exp(i phi) 
Omega = omega-1i./tau;			% complex frequency: Omega = omega - i/tau 



% old fits  
%A_abs(1) = 776.8336*eta^4 + -355.6494*eta^3 + 60.4432*eta^2 + -1.3418*eta + 0.1108;
%A_abs(2) = 2925.2033*eta^4 + -1175.1857*eta^3 + 194.9479*eta^2 + -8.9058*eta + 0.3349;
%A_abs(3) = 4703.3909*eta^4 + -1837.4023*eta^3 + 303.8194*eta^2 + -17.4766*eta + 0.5330;
%
%A_arg(1) = 29.7941*eta^3 + -1.4878*eta^2 + 60.4432*eta + -1.3418;
%A_arg(2) = -82.5869*eta^3 + 36.3616*eta^2 + 194.9479*eta + -8.9058;
%A_arg(3) = -120.7687*eta^3 + 40.8486*eta^2 + 303.8194*eta + -17.4766;



