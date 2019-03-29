function [ddotr dr_dt dprstar_dt] = EOBddotrNoFlux(nu,y,dydt,EOBopt,EOBmet,EOBHam)

%EOBddotrBoot2.m Computes an approximate value for \ddot{r}.
%
%   The second-order time derivative of the EOB radius r is computed using
%   the Poisson's bracket:
%
%    \ddot{r} = -\left\{H,\dot{r}\right\} + \dfrac{\partial\dot{r}}{\partial t}
%
%   (where H is the EOB (normalized) Hamiltonian) 
%   A direct calculation without the flux is employed.
%
%   [ddotr dr_dt] = EOBddotrNoFlux(nu,y,dydt,EOBopt,EOBmet,EOBHam)
%


% Unpack y
r      = y(2,:);
pph    = y(3,:);
prstar = y(4,:);

sz     = size(r);

% r.h.s. of the conservative part of two EoM
dprstar_dt = dydt(4,:); 
dr_dt      = dydt(2,:);

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

% Shorthands
sqrtAbyB = sqrt(A./B); 
one_A    = 1.d0./A;
one_B    = 1.d0./B;

E        = nu.*H;
E2       = E.*E;
tmpE     = 1./Heff+nu./E2;
denE     = E.*Heff;
one_denE = 1./denE;

z3       = 2.0*nu*(4.0-3.0*nu);
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
prstar2 = prstar.^2;
prstar3 = prstar.^3;
prstar4 = prstar.^4;

% Computing dpr*_dt, Not necessary: given by dydt. DELME.
%dprstar_dt  = -0.5d0*sqrtAbyB.*( pph.^2.*u2.*(dA-2d0*A.*u) + dA ...
%              + z3*(dA.*u2 - 2.0*A.*u3).*prstar4 ).*one_denE;

% Derivatives of \dot{r} wrt to r, prstar and p_phi
ddotr_dr      = sqrtAbyB.*( (prstar + z3.*2.*A.*u2.*prstar3)...
                .*(0.5.*(dA.*one_A-dB.*one_B)-dHeff_dr.*tmpE)...
                + 2.0*z3*(dA.*u2 - 2*A.*u3).*prstar3).*one_denE;
            
ddotr_dprstar = sqrtAbyB.*( 1+z3*6.d0*A.*u2.*prstar2...
                 -(prstar + z3.*2*A.*u2.*prstar3).*dHeff_dprstar.*tmpE).*one_denE;
            

% Approximate ddot(r) without Flux
ddotr = dprstar_dt.*ddotr_dprstar + dr_dt.*ddotr_dr;


