function [Fphi, Frstar] = EOBRRDIN_mod(nu, dydt, y, EOBopt, EOBmet, EOBHam)

%EOBRRDIN Compute radiation reaction with resummed flux. 
%   
%   Damour,Iyer & Nagar RR for  circular orbits.
%
%   [Fphi, Frstar] = EOBRRDIN(nu, dydt, y, EOBopt, EOBmet, EOBHam)
%
%   Reference(s)
%   Damour, Iyer & Nagar, PRD 79, 064004 (2009) [azimuthal force]
%   Bini & Damour, PRD 86 (2012), 124012        [radial force]


% Compute ddot r (approximate value)
[ddotr dr_dt] = EOBopt.Computeddotr(nu, y,dydt, EOBopt, EOBmet, EOBHam);                


% Compute the RR
Omega   = dydt(1);
[Fphi, Frstar] = EOBRRDINddotr_mod(nu, Omega, y, EOBopt, EOBmet, EOBHam, ddotr);

    