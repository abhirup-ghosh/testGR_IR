function [flso, flRing] = calcflso(M)
% 
% CALCFLSO - calculate the GW frequency corresponding to the last-stable
% orbit and the light ring in the Scwarzschild geometry (test particle 
% limit)
% 
% usage: [flso, flRing] = calcflso(M)
% 
% M         - total mass of the binary in solar masses
% flso      - GW frequency corresponding to the LSO. 
% flRing    - GW frequency corresponding to light ring 
% 
% P. Ajith, 21.07.06
% 
% $Id: calcflso.m 120 2010-10-20 03:01:35Z ajith $

MSOLAR_TIME = 4.92579497077314e-06;

% convert the mass into seconds
M = M*MSOLAR_TIME;

% schwarzschild LSO
vlso = sqrt(1/6);
    
% schwarzschild light ring
vlRing = sqrt(1/3);

% calculate the GW freq corresponding to the LSO
flso = vlso.^3./(pi*M);

% calculate the GW freq corresponding to the light ring
flRing = vlRing.^3./(pi*M);


