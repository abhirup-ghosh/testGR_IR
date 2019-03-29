function dydt = EOBRhs(t,y, nu, EOBopt)

%EOBRhs Compute r.h.s. of the full EOB equations.
%
%   dydt = EOBRhs(t,y, nu, EOBopt) return the r.h.s. of EOB equations
%   written using prstar, i.e. the conjugate momentum to the generalized
%   tortoise coordinate. 
%
%   Complete reference: Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013)
%


PNorder = EOBopt.PNorder;


% Allocate memory
dydt = zeros(4,1);


% Unpack y
phi    = y(1);
r      = y(2);
pph    = y(3);
prstar = y(4);


% Metric (do not compute 2nd drvts of A)
Metric = EOBMetric(nu,r, EOBopt, 0);
A      = Metric.A;
B      = Metric.B;
dA     = Metric.dA;


% Energy
% H      = \hat{H}       = H/mu
% Heff   = \hat{H}_{eff} = H_{eff}/\mu
Ham = EOBopt.ComputeEOBHam(nu, r,pph,prstar, EOBopt, Metric);

H    = Ham.H;
Heff = Ham.Heff;
E    = nu*H;


% Shorthands 
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
pph2     = pph.^2;
prstar3  = prstar.^3;
sqrtAbyB = sqrt(A./B);


% d\phi/dt 
dydt(1) = A*pph/(r^2*E*Heff );


% dr/dt and (conservative part of) dp_{r*}/dt 
if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')
    
    dydt(2) =  sqrtAbyB.*prstar./(E*Heff);
    dydt(4) = -sqrtAbyB*( pph2.*u2*(dA-2d0*A*u) + dA)/(2.d0*E*Heff);
    
else
    
    dydt(2) =   sqrtAbyB.*(prstar+4.0*nu*(4.0-3.0*nu)*A.*u2.*prstar3)./(E*Heff);
    dydt(4) = - sqrtAbyB*( pph2.*u2*(dA-2d0*A.*u) + dA ...
        + 2.0*nu*(4.0-3.0*nu)*(dA*u2 - 2.0*A.*u3)*prstar^4 )/(2.d0*E*Heff);        
    
end


% Compute radiation reaction force
[Fphi, Frstar] = EOBopt.ComputeRRForce(nu, dydt, y, EOBopt, Metric, Ham);


% dp_{\phi}/dt 
dydt(3) = Fphi;


% dp_{r*}/dt    
dydt(4) = dydt(4) + Frstar;    








