function y0 = EOBIdPostAdiab0( nu, r0, EOBopt )

%EOBIdPostAdiab0 Provides initial data in post-circular (or post-adiabatic) 
%   approximation for point-particle dynamics (nu<<1). 
%
%   For nu<<1, the post-circular approximation is sufficient to
%   yield an initial configuration with negligible 
%   eccentricity. Note that for generic nu one computes initial
%   data in the post-post-circular approximation in a different routine
%   (EOBIdPostAdiab(nu,r0,EOBopt)
%
%   Two steps procedure:
%    1. Compute j                      => circular ID, j!=0, pr =0
%    2. From j, compute pr*            => post circular ID, j!=0, pr!=0  
%
%   Y0 = EOBIdPostAdiab0( nu, r0, EOBopt ) given the radius, return
%   the initial data vector Y0 = [0,r0,pph,prstar,pr,j] 
%
%   Reference(s)
%


use_DIN_RR = strcmp(EOBopt.RadReac,'din');


% Compute everything in a neighbourhood of r0
% (take a finite diff drvt later, quick way)
N  = 6;
dr = 1e-8;
r  = r0-(N-1)*dr:dr:r0+N*dr;
r2 = r.^2;
r3 = r.^3;
u  = 1./r;


% Compute metric
[A dA d2A B] = EOBMetric( nu, r, EOBopt );


% Angular momentum stuff
j2   =  r3.*dA./(2*A-r.*dA);
j    =  sqrt(j2);
j3   =  j2.*j;
djdr = -j3./r3.*( 2.0 - 3.0*A./(r.*dA) - A.*d2A./(dA.^2) );


% Adiabatic energy
H0eff    = sqrt(A.*(1.0 + j2./r2));
E0       = 1.0;
Omega_j  = A.*j./(r2.*H0eff);
r_omega  = r;
v_phi    = Omega_j.*r_omega;


if use_DIN_RR
        
    % DIN resummmation of radiation reaction        
    x    = v_phi.^2;
    jhat = j./(r_omega.*v_phi);
    hatF = EOBFluxDIN(x,Omega_j,E0,H0eff,jhat,nu,r,0,0, EOBopt);
    
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
Ctmp      =  sqrt(B./A).*H0eff;
prstar    =  Ctmp.*Fphi./djdr;
pph       =  j; 


% Radial momentum conjugate to r
pr = prstar.*sqrt(B./A);


% Compute the variables at r(N)
pr     = pr(N);
prstar = prstar(N);
pph    = pph(N);
j      = j(N);


y0 = [0,r0,pph,prstar,pr,j]; 
