function [Omega E pr Heff] = EOBOmgE(nu, r,pph,prstar, EOBopt, varargin)

%EOBOmgE Computes Energy and orbital frequency along dynamics.
%
%   [Omega E pr Heff] = EOBOmgE(nu, r,pph,prstar, EOBopt) 
%
%   [Omega E pr Heff] = EOBOmgE(nu, r,pph,prstar, EOBopt, EOBmet) pass the
%   structure output by EOBMetric.m, and do not call the routine.  
%


particle_dynamics = strcmp(EOBopt.Dynamics,'particle');
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


% Shorthands
A        = Metric.A;
B        = Metric.B;
pph2     = pph.^2;
u        = 1./r;
r2       = r.^2;
u2       = u.^2;
z3       = 2.0*nu*(4.0-3.0*nu);
prstar2  = prstar.^2;
prstar4  = prstar.^4;


if particle_dynamics
    % Effective Hamiltonian (divided by mu)
    Heff  = sqrt(A.*(1.0 + pph2.*u2) + prstar2);
    % Energy (divided by M)
    E     = ones(size(Heff));                             
    Omega = A.*pph./( r2.*Heff );   
else
    % Effective Hamiltonian (divided by mu)
    if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')    
        Heff = sqrt(A.*(1.0 + pph2.*u2)+ prstar2);    
    else    
        Heff = sqrt(A.*(1.0 + pph2.*u2) + prstar2 + z3.*A.*u2.*prstar4);  
    end
    % Energy (divided by M)
    E     = sqrt( 1.0 + 2.0*nu*(Heff - 1.0) );
    Omega = A.*pph./( r2.*E.*Heff );
end

% radial momentum
pr = sqrt(B./A).*prstar;

