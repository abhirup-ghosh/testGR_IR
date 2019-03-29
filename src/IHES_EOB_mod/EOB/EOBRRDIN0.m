function [Fphi, Frstar] = EOBRRDIN0(nu, dydt, y, EOBopt, EOBmet, EOBHam)

%EOBRRDIN0 Compute radiation reaction with resummed flux for nu->0. 
%
%   Damour,Iyer & Nagar RR for circular orbits.
%
%   [Fphi, Frstar] = EOBRRDIN0(nu, omega, y, EOBopt, EOBmet, EOBHam)
%
%   Reference(s)
%   Damour,Iyer & Nagar ...
%


% Compute ddot r (approximate value) - NOT USED FOR PARTICLE !
%[ddotr dr_dt] = EOBopt.Computeddotr(nu, y,dydt, EOBopt, EOBmet, EOBHam);
ddotr = 0.;    


% Compute the RR
Omega = dydt(1);
%[Fphi, Frstar] = EOBRRDINddotr(nu,Omega,y,EOBopt,EOBmet,EOBHam,ddotr);


% specific code for particle, here below:


% Unpack y
phi    = y(1);
r      = y(2);
pph    = y(3);
prstar = y(4);


% Metric
A  = EOBmet.A;
dA = EOBmet.dA;


% Hamiltoinian
Heff = EOBHam.Heff;
H    = EOBHam.H;
E    = 1;


% Shorhands
u        = 1./r;
u2       = u.^2;


% Angular momentum flux 
rw   = r;
vw   = rw.*Omega;
jhat = pph./(rw.*vw);
x    = vw.^2; 
hatF = EOBFluxDIN(x,Omega,E,Heff,jhat,0, r,prstar,ddotr,EOBopt);


if strcmp(EOBopt.HorizonFlux,'yes')
  
  % Add horizon absorbed contribute
  hatFH = EOBFluxHorizon(x,Heff,jhat,nu, EOBopt);
  hatF  = hatF + hatFH;
    
end

% put together contributes and LO:
Fphi = -32.d0/5.d0 * nu * Omega.^5 .* rw.^4 .* hatF;


% Radial flux
Frstar = 0;

if strcmp(EOBopt.RadialFlux,'yes')

    % Add dissipation of radial momentum. 
    % Consistency argument, Bini & Damour, arXiv:1210.2834
    c1     = (-227/140*nu + 1957/1680);
    c2     = (753/560*nu^2 + 165703/70560*nu - 25672541/5080320);
    Frhat  = 1 + c1.*u + c2.*u2;
    Frstar = -5/3*prstar./pph .*Fphi.*Frhat;

end
