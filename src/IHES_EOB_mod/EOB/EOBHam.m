function varargout = EOBHam(nu, r,pph,prstar, EOBopt, varargin)

%EOBHam Computes the EOB real and effective Hamiltonian for EOB dynamics
%   (generic nu). 
%
%   [H Heff dHeff_dr dHeff_dprstar dHeff_dpphi] = EOBHam(nu, r,pph,prstar, EOBopt)
%      H             = real EOB Hamiltonian divided by mu=m1m2/(m1+m2)
%      Heff          = effective EOB Hamiltonian (divided by mu)
%      dHeff_dr      = drvtv of Heff wrt r
%      dHeff_dprstar = drvtv of Heff wrt prstar
%      dHeff_dprstar = drvtv of Heff wrt pphi
%
%   [H Heff dHeff_dr dHeff_dprstar] = EOBHam(nu,r,pph,prstar,EOBopt,EOBmet)
%   pass the structure output by EOBMetric.m, and do not call the routine.  
%
%   EOBHAM = EOBHam( nu, ... ) return a structure with the Hamiltonian
%   variables and derivatives. 
%


PNorder = EOBopt.PNorder;


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
z3      = 2.0*nu*(4.0-3.0*nu);
pph2    = pph.^2;
u       = 1./r;
u2      = u.^2;
u3      = u.^3;
prstar2 = prstar.^2;
prstar3 = prstar.^3;
prstar4 = prstar.^4;

A       = Metric.A;
dA      = Metric.dA;


% Ham
if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')
    
    Heff          = sqrt(A.*(1.d0 + pph2.*u2) + prstar2);
    dHeff_dr      = 0.;
    dHeff_dprstar = 0.;
    dHeff_dpphi   = 0.;
    
else
    
    Heff          = sqrt(A.*(1.d0 + pph2.*u2) + prstar2 + z3*A.*u2.*prstar4);
    divHeff       = 1./Heff;     
    dHeff_dr      = 0.5*(dA + (pph2 + z3*prstar4).*(dA.*u2 - 2*A.*u3)).*divHeff;
    dHeff_dprstar = (prstar + z3*2.0*A.*u2.*prstar3).*divHeff;
    dHeff_dpphi   =  A.*pph.*u2.*divHeff;
    
end

H = sqrt( 1d0 + 2d0*nu*(Heff - 1) )/nu;
    

% Manage args out
varargout = SetVOutput( nargout, H,Heff,dHeff_dr,dHeff_dprstar,dHeff_dpphi ); 
