function [y dy]=EOBQNMTemplate(ck,t, sigmak)

%EOBTemplateQNM Compute QNM waveform.
%
%   [y dy]=EOBQNMTemplate(c,t, sigma) given complex coefficients sigma(k)
%   and c(k), and the time vector t return 
%
%   y = sum_k c_k exp( - sigma_k t )
%   
%   and the derivative
%

[n,m] = size(t);

nt = length(t);
nc = length(sigmak);

t      = reshape(t,1,nt);
sigmak = reshape(sigmak,nc,1);
ck     = reshape(ck,nc,1);

y = ck * ones(1,nt);
y = y .* exp(-sigmak*t);      
y = sum(y);

y = reshape(y,n,m); % give back the same size as t

if nargout==2
    
    dk = sigmak.*ck;    
    dy = - dk * ones(1,nt);
    dy = dy .* exp(-sigmak*t);
    dy = sum(dy);

    dy = reshape(dy,n,m);
    
end
    