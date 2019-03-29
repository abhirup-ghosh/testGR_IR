function [xm,fm,dfm] = FindMax(x,f,n)

%FindMax Compute extremum of f(x)
%
%   [xm,fm,dfm] = findmax(x,f,n) given a function f(x) tabulated at points
%   x, return the extremum ( f'=0 ) of the fitting polynomial of degree n
%

f = reshape(f,size(x));
[p,s,mu] = polyfit(x,f,n);
dp       = polyder(p);
[xm,dfm] = fzero(@(x) polyval(dp,x,s,mu),(x(end)-x(1))/2);
fm       = polyval(p,xm,s,mu);

