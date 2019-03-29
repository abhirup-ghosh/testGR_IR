function [Fphi, Frstar] = EOBRRPN(nu, dydt, y, EOBopt, EOBmet, EOBHam)

%EOBRRno Compute radiation reaction with PN-Taylor flux.
%
%   [Fphi, Frstar] = EOBRRPN(nu, dydt, y, EOBopt, EOBmet, EOBHam)
%

Omega   = dydt(1);
Omegao3 = Omega.^(1/3);
hatF    = EOBFluxTaylor(nu,Omegao3); 
Fphi    = - 32.d0/5.d0 * nu * Omegao3.^7 .* hatF;
Frstar  = 0; 