function [Fphi, Frstar] = EOBRRDINddotr_mod(nu, Omega, y, EOBopt, EOBmet, EOBHam, ddotr)

%EOBRRDINddotr Compute radiation reaction with resummed flux for a given ddotr.
%
%   Damour,Iyer & Nagar RR for circular orbits.
%   [Fphi, Frstar] = EOBRRDIN(nu, Omega, y, EOBopt, EOBmet, ddotr)
%
%   Reference(s)
%   Damour,Iyer & Nagar ...
%

    
PNorder = EOBopt.PNorder;
%particle_dynamics = strcmp(EOBopt.Dynamics,'particle');
    
    
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
E    = nu.*H;


% Shorhands
r2       = r.^2;
r3       = r.^3;
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
u4       = u.^4;
u5       = u.^5;
pph2     = pph.^2;
prstar2  = prstar.^2;  


% Angular momentum flux
hatF = 0;
%if particle_dynamics 
% % coded in EOBRRDIN0.m
%  psi  = 1.;
%  rw   = r;
%  vw   = rw.*Omega;
%  jhat = pph./(rw.*vw);
%  x    = vw.^2; 
%  hatF = EOBFluxDIN(x,Omega,E,Heff,jhat,0, r,prstar,ddotr, EOBopt);
%else
sqrtWmo = sqrt(A.*(1 + pph2.*u2))-1;
if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')           
  psi = (1 + 2*nu*sqrtWmo )./(1 - 3*nu.*u2);    
else    
  psi = 2*(1 + 2*nu*sqrtWmo)./(r2.*dA);
end
rw   = r.*psi^(1/3);
vw   = rw.*Omega;
jhat = pph./(rw.*vw);    
x    = vw.^2;
hatF = EOBFluxDIN_mod(x,Omega,E,Heff,jhat,nu, r,prstar,ddotr, EOBopt);
%end

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
