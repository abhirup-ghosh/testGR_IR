function W = EOBPotential(nu, r, j, EOBopt, varargin)

%EOBPotential Compute the adiabatic EOB effective potential.
%
%   W = EOBPotential(nu, r, j, EOBopt) calls EOBMetric.m for the
%   computation of the metric potential A.
%
%   W = EOBPotential(nu, r, j, EOBopt, EOBmet ) pass the structure output
%   by EOBMetric.m, and do not call the routine. 
%

%FIXME documentation, which Eq ?


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

A  = Metric.A;
j2 = j^2;
u2 = (1./r).^2;     

% Effective "EOB" squared potential
W2 = A.*(1+j2.*u2);

% "Real" radial potential
W = sqrt(1+2*nu*(sqrt(W2)-1));
