function [Mbh abh] = EOBFinalBH(nu)

%EOBFinalBH Compute mass and spin of the final black hole. 
%
%   [Mbh abh] = EOBFinalBH(nu)
%
%   Reference(s)
%   Pan et al. (xxx)
%

Mbh = 1+(sqrt(8/9)-1)*nu - 0.4333*nu^2 - 0.4392*nu^3;
abh = sqrt(12)*nu - 3.871*nu^2 + 4.028*nu^3;


