function f = EOBfLSO(u, nu, EOBopt)

%EOBfLSO Compute function for adiabatic LSO determination.
%
%   F = EOBfLSO( u, nu, EOBopt ) returns f = A'(u) B''(u) - A''(u) B'(u)
%   from EOB metric functions.
%


% Compute the EOB metric
Metric = EOBMetric( nu, 1./u, EOBopt );

A   = Metric.A;
dA  = Metric.dA_u;
d2A = Metric.d2A_u;

%dB  = Metric.dB;
% which should be:
u2  = u.^2;
dB  = u2.*dA + 2*A.*u;

% this has to be computed anyway:
d2B = d2A.*u2 + 4*u.*dA+2*A;

f   = dA.*d2B - d2A.*dB;



