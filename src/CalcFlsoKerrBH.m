function [flso,flRing,fQNM, Q, rlso, vlso, flsoFit] = CalcFlsoKerrBH(M, chi)
% 
% CalcFlsoKerrBH - compute various characteristic frequencies of a test
% particle inspiralling into a Kerr BH. 
% 
% usage:  [flso,flRing,fQNM, Q, rlso, vlso, flsoFit] = CalcFlsoKerrBH(M, chi)
% 
% M         - mass of the big BH (in units of M_sun)
% chi       - dimensionless spin parameter of the big BH
% flso      - GW frequency corresponding to the LSO. 
% flRing    - GW frequency corresponding to light ring 
% fQNM      - dominant quasi-normal-mode freqeuncy 
% Q         - quality factor of the ringing
% 
% P. Ajith, 21.07.06
% 
% $Id: CalcFlsoKerrBH.m 120 2010-10-20 03:01:35Z ajith $

MSOLAR_TIME = 4.92579497077314e-06;

% convert the mass into seconds
M = M*MSOLAR_TIME;

% location of the ISCO 
z1 = 1+(1-chi.^2).^(1/3).*((1+chi).^(1/3) + (1-chi).^(1/3));
z2 = sqrt(3*chi.^2 + z1.^2);

% if abs(chi) > 0
%     spinAlign = chi./abs(chi);
% else 
%     spinAlign = 1;
% end

spinAlign = ones(size(chi));
kerrIdx = find(abs(chi) > 0);
spinAlign(kerrIdx) = chi(kerrIdx)./abs(chi(kerrIdx));

rlso = 3.+z2 - spinAlign.*sqrt((3.-z1).*(3.+z1+2.*z2));

% location of the light ring 
rlRing = 2*(1+cos(2/3*acos(-chi)));

ulso = sqrt(1./rlso);
ulRing = sqrt(1./rlRing);

vlRing = ulRing.*(1-chi.*ulRing.^3+chi.^2.*ulRing.^6).^(1/3);
vlso = ulso.*(1-chi.*ulso.^3+chi.^2.*ulso.^6).^(1/3);
    
% calculate the GW freq corresponding to the LSO
flso = vlso.^3./(pi*M);
    
% calculate the GW freq corresponding to the light ring
flRing = vlRing.^3./(pi*M);

% calculate the ring-down frequency 
fQNM = (1. - 0.63*power(1.-chi, 0.3))/(2.*pi*M);

% calculate the quality factor of the ringing 
Q = 2*power(1.-chi, -0.45);

% simple fit to the ISCO frequency in terms of the spin (Ajith, 13.09.09)
flsoFit = (1 - 4.4547*(1-chi).^0.217 + 3.521*(1-chi).^0.26)/(pi*M);


