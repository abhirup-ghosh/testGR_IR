function dydt = T4Rhs(t,y, nu)

%T4Rhs Compute r.h.s. of the Taylor T4 approximation
%
%   dydt = T4Rhs(t,y, nu) 
%
%   Reference(s)
%    Boyle et al, XXX YYY Eq. (46).
%


% Allocate memory
x    = y(1,:);
nx   = length(x);
dydt = zeros(2,nx);


% Shorthands
gamma  = 0.5772156649015329;
nu2    = nu^2;
nu3    = nu^3;
pi2    = pi^2;
x2     = x.^2;
x3     = x.^3;
x5     = x.^5;
x32    = x.^(3/2);
x52    = x.^(5/2);
x72    = x.^(7/2);


% d\Phi/dt 
dydt(2,:) = x32;


% dx/dt 
dydt(1,:)  = 64/5*nu*x5.*(1 - (743/336+11/4*nu).*x + 4*pi*x32 ...
        + (34103/18144 + 13661/2016*nu + 59/18*nu2).*x2...
        -(4159/672 + 189/8*nu)*pi.*x52...
        +( 16447322263/139708800 -1712/105*gamma -56198689/217728*nu...
         +541/896*nu2 -5605/2592*nu3 + pi2/48*(256+451*nu)...
         -856/105*log(16*x)).*x3...
        +(-4415/4032 + 358675/6048*nu + 91495/1512*nu2)*pi*x72);
   




