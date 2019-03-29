function y0 = EOBIdPostAdiab( nu, r0, EOBopt )

%EOBIdPostAdiab Provides initial data in post-post-circular approximation for
%   EOB dynamics (generic nu).     
%
%   Three steps procedure:
%    1. Compute j                      => circular ID, j!=0, pr =0
%    2. From j, compute pr*            => post circular ID, j!=0, pr!=0  
%    3. From pr* and j, re-compute pph => post-post-circular ID,
%    pph0!=j!=0, pr!=0 
%
%   Y0 = EOBIdPostAdiab( nu, r0, EOBopt ) given the radius, return the
%   initial data vector Y0 = [0,r0,pph,prstar,pr,j] 
%
%   Reference(s)
%


use_DIN_RR = strcmp(EOBopt.RadReac,'din');
PNorder    = EOBopt.PNorder;


% Compute everything in a neighbourhood of r0
% (take a finite diff drvt later, quick way)
N  = 6;
dr = 1e-8;
r  = r0-(N-1)*dr:dr:r0+N*dr;
r2 = r.^2;
r3 = r.^3;
z3 = 2.0*nu*(4.0-3.0*nu);
u  = 1./r;


% Compute metric
[A dA d2A B dB] = EOBMetric( nu, r, EOBopt );


% Angular momentum for circular orbit: circular ID
j2   =  r3.*dA./(2*A-r.*dA);
j    =  sqrt(j2);
j3   =  j.^3;
djdr = -j3./r3.*( 2.0 - 3.0*A./(r.*dA) - A.*d2A./(dA.^2) );


% For circular orbit at r0=r(N)
H0eff    = sqrt(A.*(1.d0 + j2./r2));            % effective Hamiltonian H_0^eff
E0       = sqrt(1.0 + 2.0*nu*(H0eff - 1.0) );  % real Hamiltonian      H_0
H0       = E0/nu;                               % H_0/nu
Omega_j  = A.*j./(nu*r2.*H0.*H0eff);            % Orbital frequency (from Hamilton's equation)
psi      = 2*(1.0 + 2.0*nu*(H0eff - 1.0))./(r2.*dA); % correction factor to the radius
r_omega  = r.*psi.^(1.0/3.0);                         % EOB-corrected radius 
v_phi    = Omega_j.*r_omega;                          % "corrected" azimuthal velocity such that 
                                                      % Kepler's law is satisfied, r_omg^3 Omg_i^2 = 1


if use_DIN_RR
        
    % DIN resummmation of radiation reaction    
    % NOTE: the flux here has pr_star=ddot(r)=0 for simplicity.                
    x    = v_phi.^2;
    jhat = j./(r_omega.*v_phi); % Newton-normalized angular momentum
    hatF = EOBFluxDIN(x,Omega_j,E0,H0eff,jhat,nu,r, 0,0, EOBopt);
    
    if strcmp(lower(EOBopt.HorizonFlux),'yes') 
    
        % Add horizon absorbed contribute
        hatFH = EOBFluxHorizon(x,H0eff,jhat,nu, EOBopt);
        hatF  = hatF + hatFH;

    end    
    
    Fphi = -32.0/5.0*nu * Omega_j.^5 .*  r_omega.^4 .*hatF;

else 
        
    % Taylor form    
    hatF = EOBFluxTaylor(v_phi,nu);    
    Fphi = -32.0/5.0*nu * Omega_j.^(7/3) .* hatF;
    
end


% Radial momentum conjugate to r*: post-circular ID
Ctmp     =  sqrt(B./A)*nu.*H0.*H0eff;
prstar   =  Ctmp.*Fphi./djdr;


% Angular momentum again: post-post-circular ID
rhs       = prstar;;
c1        = (-227/140*nu + 1957/1680);
c2        = (753/560*nu^2 + 165703/70560*nu - 25672541/5080320);
Frhat     = 1 + c1.*u + c2.*u.^2;  
Frstar    = -5/3.*Fphi.*prstar./j .*Frhat;
dprstardt = Fphi./djdr.*FDdrvt(rhs,r,4) - Frstar;  % NOTE: Fr* here


if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')
    pph = j.*sqrt(1 + 2.*Ctmp./dA.*dprstardt);
else
    pph = j.*sqrt(1 + 2.*Ctmp./dA.*dprstardt - z3.*prstar.^4./j2);
end


% Radial momentum conjugate to r
pr = prstar.*sqrt(B./A);


% Compute the variables at r0=r(N)
pr     = pr(N);
prstar = prstar(N);
pph    = pph(N);
j      = j(N);

y0     = [0,r0,pph,prstar,pr,j]; 




