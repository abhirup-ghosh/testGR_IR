function LSO = EOBAdiabLSO(nu, EOBopt)

%EOBAdiabLSO Compute Last-Stable-Orbit (LSO) in the EOB formalism.
%
%   LSO = EOBAdiabLSO( nu, EOBopts ) return a structure with LSO data
%


% LSO
LSO.u = fzero( @(x) EOBfLSO(x, nu, EOBopt), 1/5.7 );
LSO.r = 1./LSO.u;


% Metric potentials
Metric = EOBMetric(nu, LSO.r, EOBopt);
A      = Metric.A;
dA     = Metric.dA;


% Angular momentum
j2    = LSO.r.^3.*dA./(2*A-LSO.r.*dA);
LSO.j = sqrt(j2);

LSO.Omega = EOBOmgE(nu, LSO.r,LSO.j,0, EOBopt, Metric);
LSO.w22   = 2*LSO.Omega;
    

% Newtonian-like connection 
LSO.Newt_w22   = 2*LSO.u^(3/2);
LSO.Newt_Omega = LSO.u.^(3/2);

