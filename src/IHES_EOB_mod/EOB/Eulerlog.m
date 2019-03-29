function y = Eulerlog(x,m)

%EULERLOG Compute Eulerlog function.
%
%   y = Eulerlog(x,m)
%

EulerGamma = 0.57721566490153286061;
y = EulerGamma + log(2.*m*sqrt(x));

