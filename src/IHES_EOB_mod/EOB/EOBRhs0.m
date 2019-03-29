function dydt = EOBRhs0(t,y, nu, EOBopt)

%EOBRhs0 Compute r.h.s. of the full EOB equations in the particle limit (nu<<1). 
%
%   dydt = EOBRhs0(t,y, nu, EOBopt) return the r.h.s. of EOB equations
%   written using prstar, i.e. the conjugate momentum to the generalized
%   tortoise coordinate. 
%


% Allocate memory
dydt = zeros(4,1);


% Unpack y
phi    = y(1);
r      = y(2);
pph    = y(3);
prstar = y(4);


% Shorhands 
r2        = r.^2;
u         = 1./r;
u2        = u.^2;
%u3        = u.^3;
%u4        = u.^4;
%u5        = u.^5;
pph2      = pph.^2;
%prstar2   = prstar.^2;
%prstar3   = prstar.^3;


% Metric
A   = 1 - 2.*u;        
dA  = 2*u2;

Metric.A = A;      
Metric.dA = dA;      


% Energy
Ham = EOBopt.ComputeEOBHam(nu, r,pph,prstar, EOBopt, Metric);
Heff = Ham.Heff;
%Heff = sqrt(A.*(1+pph2./r2)+prstar2);


% d\phi/dt 
dydt(1) = A*pph/( r2*Heff );
    

% dr/dt 
dydt(2) = A*prstar/Heff;


% dp_{r*}/dt (conservative)
dydt(4) = -A*( pph2.*u2.*(dA-2.0*A.*u) + dA )/(2.0*Heff);
    

% Compute radiation reaction force
[Fphi, Frstar] = EOBopt.ComputeRRForce(nu, dydt, y, EOBopt, Metric, Ham);                    


% dp_{\phi}/dt 
dydt(3) = Fphi;


% dp_{r*}/dt    
dydt(4) = dydt(4) + Frstar;   


    
    