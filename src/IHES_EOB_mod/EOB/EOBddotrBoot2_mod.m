function [ddotr dr_dt] = EOBddotrBoot2_mod(nu,y,dydt,EOBopt,EOBmet,EOBHam)

%EOBddotrBoot2.m Computes an approximate value for \ddot{r}.
%
%   The second-order time derivative of the EOB radius r is computed using
%   the Poisson's bracket:
%
%    \ddot{r} = -\left\{H,\dot{r}\right\} + \dfrac{\partial\dot{r}}{\partial t}
%
%   (where H is the EOB (normalized) Hamiltonian) 
%   A two steps bootstrap procedure which calls the flux is employed. 
%
%   [ddotr dr_dt] = EOBddotrBoot2(nu,y,dydt,EOBopt,EOBmet,EOBHam)
%


% Unpack y
[nv,nx] = size(y);
phi    = y(1,:);
r      = y(2,:);
pph    = y(3,:);
prstar = y(4,:);

sz    = size(r);

% r.h.s. of the conservative part of two EoM
dprstar_dt_0 = dydt(4,:);
dr_dt        = dydt(2,:);


% Metric
A  = reshape(EOBmet.A,sz);
dA = reshape(EOBmet.dA,sz);
B  = reshape(EOBmet.B,sz);
dB = reshape(EOBmet.dB,sz);


% Hamiltonian
Heff          = reshape(EOBHam.Heff,sz);
H             = reshape(EOBHam.H,sz);
dHeff_dr      = reshape(EOBHam.dHeff_dr,sz);
dHeff_dprstar = reshape(EOBHam.dHeff_dprstar,sz);
dHeff_dpphi   = reshape(EOBHam.dHeff_dpphi,sz);


% Shorhands
sqrtAbyB  = sqrt(A./B); 
one_A    = 1.d0./A;
one_dA   = 1.d0./dA;
one_B    = 1.d0./B;

E        = nu.*H;
E2       = E.*E;
tmpE     = 1./Heff+nu./E2;
denE     = E.*Heff;
one_denE = 1./denE;

z3       = 2.0*nu*(4.0-3.0*nu);
pph2     = pph.^2;
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
prstar2 = prstar.^2;
prstar3 = prstar.^3;


% Derivatives of \dot{r} wrt to r, prstar and p_phi
ddotr_dr      = sqrtAbyB.*( (prstar + z3.*2.*A.*u2.*prstar3)...
                .*(0.5.*(dA.*one_A-dB.*one_B)-dHeff_dr.*tmpE)...
                + 2.0*z3*(dA.*u2 - 2*A.*u3).*prstar3).*one_denE;
            
ddotr_dprstar = sqrtAbyB.*( 1+z3*6.d0*A.*u2.*prstar2...
                 -(prstar + z3.*2*A.*u2.*prstar3).*dHeff_dprstar.*tmpE).*one_denE;
            
ddotr_dpphi   = -dr_dt.*tmpE.*dHeff_dpphi;            


%  Bootstrap two-step procedure to approximately compute ddotr

% 0. Approximate ddot(r) without Fphi
ddotr = dprstar_dt_0.*ddotr_dprstar + dr_dt.*ddotr_dr;

% Compute RR force with 1.
Omega          = A.*pph.*u2.*one_denE;
[Fphi, Frstar] = EOBRRDINddotr_mod(nu, Omega, y, EOBopt, EOBmet, EOBHam, ddotr);

% 1. Approximate ddot(r) with all three terms
dprstar_dt = dprstar_dt_0 + Frstar;
ddotr      = dprstar_dt.*ddotr_dprstar + dr_dt.*ddotr_dr + ddotr_dpphi.*Fphi;

% Recompute the RR force
[Fphi, Frstar] = EOBRRDINddotr_mod(nu, Omega, y, EOBopt, EOBmet, EOBHam, ddotr);

% 2. Compute \ddot{r} from the three separate contribution
dprstar_dt = dprstar_dt_0 + Frstar; 
ddotr      = dprstar_dt.*ddotr_dprstar + dr_dt.*ddotr_dr + ddotr_dpphi.*Fphi;









