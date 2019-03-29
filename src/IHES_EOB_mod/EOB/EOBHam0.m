function varargout = EOBHam0(nu, r,pph,prstar, EOBopt, varargin)

%EOBHam0 Computes the EOB real and effective Hamiltonian for particle
%   dynamics (nu=0).
%
%   [H Heff dHeff_dr dHeff_dprstar dHeff_dpphi] = EOBHam0(nu, r,pph,prstar, EOBopt)
%      H             = real EOB Hamiltonian divided by mu=m1m2/(m1+m2)
%      Heff          = effective EOB Hamiltonian (divided by mu)
%      dHeff_dr      = drvtv of Heff wrt r
%      dHeff_dprstar = drvtv of Heff wrt prstar
%      dHeff_dprstar = drvtv of Heff wrt pphi
%
%   [...] = EOBHam0(nu,r,pph,prstar,EOBopt,EOBmet) pass the structure
%   output by EOBMetric.m, and do not call the routine.   
%
%   EOBHAM = EOBHam0( nu, ... ) return a structure with the Hamiltonian
%   variables and derivatives. 


% Manage args in
na = length(varargin);
if (na>1)
    error('too many input args')
elseif na==1
    Metric = varargin{1};    
else
    % Compute the EOB metric
    Metric = EOBMetric( nu, r, EOBopt );
end


% Shorthands & constants
pph2    = pph.^2;
u       = 1./r;
u2      = u.^2;
prstar2 = prstar.^2;


%A       = Metric.A; % Metric not needed
A       = 1-2*u;


% H    
Heff  = sqrt(A.*(1. + pph2.*u2) + prstar2);
H     = ones(size(Heff)); % H/nu

dHeff_dr      = zeros(size(Heff));
dHeff_dprstar = zeros(size(Heff));
dHeff_dpphi   = zeros(size(Heff));
    

% Manage args out
varargout = SetVOutput(nargout, H,Heff,dHeff_dr,dHeff_dprstar,dHeff_dpphi); 

