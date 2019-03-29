function rho22 = EOBrho22Pade23(x,nu)

%EOBrho22Pade23 Computes the (2,3) Padé approximant to the rho_22. 
%
%   rho22 = EOBrho22Pade23(x,nu)
%


% Shorthands & constants 
one        =  1.0;
two        =  2.0;
EulerGamma =  0.5772156649015329;
eulerlog22 =  EulerGamma + two*log(two) + 0.5d0*log(x);


% Coefficients of the Taylor expansion (at 5PN) of rho22
% after the square-root resummation. It is:
%
% rho22 = 1 + f1*x + f2*x^2 + f3*x^3 + f4*x^4 + f5*x^5
%
% Note the use of the eulerlog22 function
f1  = -43/42 + 55*nu/84;
f2  = -20555/10584 - 33025/21168*nu + 19583/42336*nu^2;
f3  =  1556919113/122245200 - 428/105*eulerlog22     - 48993925*nu/9779616 ...
      -6292061*nu^2/3259872 + 10620745*nu^3/39118464 + 41*nu*pi^2/192;
f4  = -387216563023/160190110080 + 9202*eulerlog22/2205;
f5  = -16094530514677/533967033600 + 439877/55566*eulerlog22;


% Coefficients of the (2,3) Padé approximant 
c1 =  -f1;
c2 =   f1 - f2/f1;
c3 =  (f1 .* f3-f2.^2)./(f1.*(f1.^2-f2));
c4 =  -f1 .* (f2.^3+f3.^2+f1.^2.*f4-f2.*(two.*f1.*f3+f4))./((f1.^2-f2)*(f1.*f3-f2.^2));
c5 =  -(f1.^2-f2).*(-f3.^3+two.*f2.*f3.*f4-f1.*f4.^2-f2.^2.*f5+f1.*f3.*f5)...
        ./ ( (f1.*f3-f2.^2).*(f2.^3+f3.^2+f1.^2.*f4-f2.*(two*f1.*f3+f4) ) );
    
% Compute the (2,3) Padé approximant
y     = one + c4.*x./(one+c5.*x);
y     = one + c3.*x./y;
y     = one + c2.*x./y;
y     = one + c1.*x./y;

rho22 = 1./y;  


