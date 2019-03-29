function [hatFH FlmH] = EOBFluxHorizon(x,Heff,jhat,nu, EOBopt)

%EOBFluxHorizon Compute horizon-absorbed fluxes.
%
%   [hatF FlmH] = EOBFluxHorizon(x,Heff,jhat,nu, EOBopt)
%
%   Nagar & Akcay, PRD 85, 044025 (2012)
%   Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)


L    = EOBopt.L;
LM2K = EOBopt.LM2K;


% Allocate memory
kmax   = length(L);
rhoHlm = zeros(kmax, length(x));
FlmHLO = rhoHlm;
FlmH   = rhoHlm;


% Shorthands
nu2 = nu.^2;
nu3 = nu.^3;

x2  = x.^2;
x3  = x.^3;
x4  = x.^4;
x5  = x.^5;
x9  = x.^9;
x10 = x.^10;

k22 = LM2K(2,2);
k21 = LM2K(2,1);
k33 = LM2K(3,3);
k32 = LM2K(3,2);
k31 = LM2K(3,1);
k44 = LM2K(4,4);
k43 = LM2K(4,3);
k42 = LM2K(4,2);
k41 = LM2K(4,1);


% The Newtonian asymptotic contribution
FNewt22 = 32/5*x5;
F22_1PN = 1 + (3-17*nu+25*nu2 - 8*nu3)/(1-4*nu+2*nu2)*x;


% Compute leading-order part (nu-dependent)
FlmHLO(k22,:) = 32/5*(1-4*nu+2*nu2)*x9;
FlmHLO(k21,:) = 32/5*(1-4*nu+2*nu2)*x10;


% Compute rho_lm
% Hybrid expression for tilde(rho)22
% The 1PN  term is  the exact analytical one, while the others come
% from the expansion of the (constrained) fit above.
% 4PN fractional accuracy for all rho_lm^H up to l=4.

% Coefficients are pre-computed by EOBFluxHorizonFitCoefs.m
c1 = EOBopt.HorizonFluxCoefs.c{1};
c2 = EOBopt.HorizonFluxCoefs.c{2};
c3 = EOBopt.HorizonFluxCoefs.c{3};
c4 = EOBopt.HorizonFluxCoefs.c{4};

% NOTE: the following is a polynomial evaluation, can be optimized 
rhoH(k22,:) = 1 + c1(k22)*x + c2(k22)*x2 + c3(k22)*x3 + c4(k22)*x4;
rhoH(k21,:) = 1 + c1(k21)*x + c2(k21)*x2 + c3(k21)*x3 + c4(k21)*x4;
%{
rhoH(k33,:) = 1 + c1(k33)*x + c2(k33)*x2 + c3(k33)*x3 + c4(k33)*x4;
rhoH(k32,:) = 1 + c1(k32)*x + c2(k32)*x2 + c3(k32)*x3 + c4(k32)*x4;
rhoH(k31,:) = 1 + c1(k31)*x + c2(k31)*x2 + c3(k31)*x3 + c4(k31)*x4;
rhoH(k44,:) = 1 + c1(k44)*x + c2(k44)*x2 + c3(k44)*x3 + c4(k44)*x4;
rhoH(k43,:) = 1 + c1(k43)*x + c2(k43)*x2 + c3(k43)*x3 + c4(k43)*x4;
rhoH(k42,:) = 1 + c1(k42)*x + c2(k42)*x2 + c3(k42)*x3 + c4(k42)*x4;
rhoH(k41,:) = 1 + c1(k41)*x + c2(k41)*x2 + c3(k41)*x3 + c4(k41)*x4;
%}


% Compute horizon multipolar flux (only l=2)
Heff2 = Heff.^2;
jhat2 = jhat.^2;

FlmH(k22,:) = FlmHLO(k22,:) .* Heff2 .*rhoH(k22,:).^4;
FlmH(k21,:) = FlmHLO(k21,:) .* jhat2 .*rhoH(k21,:).^4;

%{
% PN-expanded flux (check the influence on phasing)
FlmH(k22,:) = FH_LO(k22,:).*F22_1PN;
FlmH(k21,:) = FH_LO(k21,:);
%}


% Sum over multipoles and normalize to the 22 Newtonian multipole
hatFH = sum(FlmH)./FNewt22;


