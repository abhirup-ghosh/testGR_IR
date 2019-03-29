function dydt = EOBRhsCons(t,y, nu, EOBopt, varargin)

%EOBRhs Compute conservative r.h.s. of the full EOB equations.
%
%   This function is not to be used for the ODE solver.
%   Assume size(y) = nvars x length(t).
%
%   dydt = EOBRhsCons(t,y, nu, EOBopt) return the r.h.s. of EOB equations
%   written using prstar, i.e. the conjugate momentum to the generalized
%   tortoise coordinate. 
%
%   dydt = EOBRhsCons(t,y, nu, EOBopt, EOBmet) pass the structure output by
%   EOBMetric.m, and do not call the routine. 
%
%   dydt = EOBRhsCons(t,y, nu, EOBopt, EOBmet, EOBHam) pass the structure
%   output by EOBHam.m (or EOBHam0.m), and do not call the routine. 
%


% Manage args in
na = length(varargin);
if (na>2)
    error('too many input args')
elseif na==1 && ~isempty(varargin{1})
    Metric = varargin{1};    
    % Compute the EOB Hamiltonian
    Ham = EOBopt.ComputeEOBHam(nu, y(2,:),y(3,:),y(4,:), EOBopt, Metric);    
elseif na==2 && ~isempty(varargin{1}) && ~isempty(varargin{2})
    Metric = varargin{1};    
    Ham    = varargin{2};    
else
    % Compute the EOB metric (do not compute 2nd drvts of A)
    Metric = EOBMetric(nu, y(2,:), EOBopt, 0);
    % Compute the EOB Hamiltonian
    Ham = EOBopt.ComputeEOBHam(nu, y(2,:),y(3,:),y(4,:), EOBopt, Metric);    
end


% Other options
particle_dynamics = strcmp(EOBopt.Dynamics,'particle');
PNorder = EOBopt.PNorder;


% Allocate memory
[nv,nx] = size(y);
dydt = zeros(nv,nx);

% unpack y
phi    = y(1,:);
r      = y(2,:);
pph    = y(3,:);
prstar = y(4,:);


% Metric (reshape)
A      = reshape(Metric.A,1,nx);
B      = reshape(Metric.B,1,nx);
dA     = reshape(Metric.dA,1,nx);


% Energy (reshape)
% H      = \hat{H}       = H/mu
% Heff   = \hat{H}_{eff} = H_{eff}/\mu

H    = reshape(Ham.H,1,nx);
Heff = reshape(Ham.Heff,1,nx);
E    = nu*H;

% Shorthands 
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
pph2     = pph.^2;
prstar3  = prstar.^3;
prstar4  = prstar.^4;
sqrtAbyB = sqrt(A./B);
z3       = 2.0*nu*(4.0-3.0*nu);
denE     = E.*Heff;
one_denE = 1./denE;


if particle_dynamics

    % d\phi/dt 
    dydt(1,:) = A.*pph./( r.^2.*Heff );
    
    % dr/dt 
    dydt(2,:) = A.*prstar./Heff;

    % dp_{r*}/dt (conservative)
    dydt(4,:) = -A.*( pph2.*u2.*(dA-2d0*A.*u) + dA )./(2.d0*Heff);    
    
    return
    
end

    
% d\phi/dt
dydt(1,:) = A.*pph./(r.^2.*E.*Heff );

% dr/dt and (conservative part of) dp_{r*}/dt
if strcmp(PNorder,'1pn') || strcmp(PNorder,'2pn')
    
    dydt(2,:) = sqrtAbyB.*prstar.*one_denE;
    dydt(4,:) = -0.5*sqrtAbyB.*( pph2.*u2.*(dA-2d0*A.*u) + dA).*one_denE;
    
else
            
    dydt(2,:) = sqrtAbyB.*(prstar + z3*2*A.*u2.*prstar3).*one_denE;
    
    dydt(4,:) = - 0.5d0*sqrtAbyB.*( pph2.*u2.*(dA-2*A.*u) + dA ...
              + z3*(dA.*u2 - 2.0*A.*u3).*prstar4 ).*one_denE;
    
end












