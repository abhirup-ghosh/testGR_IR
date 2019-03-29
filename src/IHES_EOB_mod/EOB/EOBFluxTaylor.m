function hatF = EOBFluxTaylor(nu,v)

%EOBFluxTaylor.m Compute the Taylor expanded GW flux at 3.5PN order for
%nu\neq 0. 
%
%
% Nothe that a tunable 4PN term is added by  
% hand. Here we use the value of this term obtained by Mroue, Kidder 
% and Teukolsky, PRD 78 044004 (2008) as input. The value that they 
% obtained by fitting to NR data is
%
% F(8) = -333.75
%
% Our best value (slightly tuned) is
%
% F(8) = -322.5
%
% This is due to a fundamental difference between MKT and our Hamiltonian:
% they use the Hamiltonian with pr, while we are using the one with pr*,
% that is consistently correct at 3PN order.
%
% USAGE: EOB_FluxTaylor(v,nu)
%
% where nu is the symmetric mass ratio and v is some velocity of the system
%   hatF = EOBFluxTaylor(v,nu)
%


nu2     = nu*nu;
nu3     = nu*nu2;
gE      = 0.57721566490153286061;  %Euler gamma
pi2     = pi*pi;

F(1) =  0.d0;
F(2) = -1247.d0/336.d0 - 35.d0/12.d0*nu;
F(3) =  4.d0*pi;
F(4) = -44711.d0/9072.d0 + 9271.d0/504.d0*nu + 65.d0/18.d0*nu2;
F(5) = -(8191.d0/672.d0 + 583.d0/24.d0*nu)*pi;  


% 3.5PN Radiation-Reaction force
fl6   = -1712d0/105d0;
F(6)  =  6643739519d0/69854400d0 + 16d0/3d0*pi2 - 1712/105*gE + ...
          (-134543/7776 + 41/48*pi2)*nu - 94403/3024*nu2 - 775/324*nu3;
F(7) =  (-16285/504 + 214745/1728*nu + 193385/3024*nu2)*pi;
F(8) =  -322.5; % hard-coded after a minimum fine-tuning

hatF = 1 + F(2).*v.^2 + F(3).*v.^3 + F(4).*v.^4 + F(5).*v.^5 + ...
          (F(6) + fl6.*log(4.*v)).*v.^6 + F(7).*v.^7 + F(8).*v.^8;


