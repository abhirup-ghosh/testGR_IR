function varargout = EOBMetric( nu, r, EOBopt, varargin )

%EOBMetric Computes the EOB metric potentials and derivatives.
%
%   [A dA d2A B dB dA_u d2A_u] = EOBMetric( nu, r, EOBopt ) 
%   return metric potentials A(r), B(r) together with first and second
%   derivatives. First and second derivatives of A(u) are also returned.
%
%   EOBMETRIC = EOBMetric( nu, r, EOBopt ) 
%   return a structure with the metric potentials and derivatives.
%
%   EOBMETRIC = EOBMetric( nu, r, EOBopt, computed2A ) 
%   if computed2A then compute second derivatives of A (1, yes)
%


% Default options
comp2nder = 1; % yes


% Manage args in
na = length(varargin);
if (na>1)
    error('too many input args')
end
optargs = {comp2nder};
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[comp2nder] = optargs{:};


% Other options
particle_dynamics = strcmp(EOBopt.Dynamics,'particle');

PNorder = EOBopt.PNorder;

a5 = EOBopt.a5;
a6 = EOBopt.a6;


% Shorthands & constants
nu2 = nu^2;
nu3 = nu^3;
nu4 = nu^4;
nu5 = nu^5;
pi2 = pi^2;
pi4 = pi^4;
pi6 = pi^6;

u   = 1./r;
u2  = u.^2;
u3  = u.^3;
u4  = u.^4;
u5  = u.^5;
u6  = u.^6;
u7  = u.^7;
u8  = u.^8;
u9  = u.^9;
u10 = u.^10;
u11 = u.^11;
u12 = u.^12;
u13 = u.^13;
u14 = u.^14;


% Init with zeros and ones
zzz = zeros(size(u));
ooo = ones(size(u));

A     = ooo;
B     = ooo;
D     = ooo;

dA    = zzz;
dB    = zzz;
d2A   = zzz;
dA_u  = zzz;
d2A_u = zzz;


if particle_dynamics
    
    % Set Schwarzschild background 
    % and exit    
    A   = 1. - 2.*u;
    B   = 1./A;
    
    dA  =  2.*u2;
    d2A = -4.*u3;
    
    dA_u = - 2*ooo;

    dB = - dA./A.^2;
    
    varargout = SetVOutput( nargout, A,dA,d2A, B,dB, dA_u,d2A_u );
    return
    
end


% EOB dynamics n~=0

% Compute A function and its derivatives at various PN (resummed) orders
% A'(r) and A''(r) are computed from A'(u) and A''(u)
switch PNorder
    
    case '1pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A    = 1. - 2.*u;        
        dA   = 2.*u2;        
        dA_u = -2.0;        

        if (comp2nder)
            
            %d2A_u = ...
            d2A   = -4.*u3;            
            
        end
                
    case '2pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A     =  1d0 - 2d0.*u + 2.d0*nu.*u3;
        dA    =  2.d0.*u2 - 6d0*nu.*u4;                
        dA_u  =  -2.0 + 6*nu.*u2;        
        
        if (comp2nder)
            
            d2A_u =  12*nu*u;
            d2A   = -4d0.*u3  + 24*nu.*u5;
                        
        end        
        
        
    case '3pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        a4 = (94.0/3.0 - 41.0/32.0*pi2)*nu;
        
        A    = (2*(4.0-nu)+(a4-16.0+8.0*nu)*u)...
            ./( 2*(4-nu)+(a4+4*nu)*u + 2*(a4+4*nu)*u2 + 4*(a4+nu^2)*u3);
        
        dA_u = (-((8 + u*(a4 + 8*(-2 + nu)) - 2*nu).*(a4 + 4*nu + 4*u*(a4 + 4*nu) + 12*u2*(a4 + nu^2)))...
            + (a4 + 8*(-2 + nu))*(8 - 2*nu + u*(a4 + 4*nu) + 2*u2*(a4 + 4*nu) + 4*u3*(a4 + nu^2)))...
            ./(8 - 2*nu + u*(a4 + 4*nu) + 2*u2*(a4 + 4*nu) + 4*u3*(a4 + nu^2)).^2;
        
        dA  = -u2.*dA_u;
        
        %{
        dA   = -u2.*(-((8 + u*(a4 + 8*(-2 + nu)) - 2*nu).*(a4 + 4*nu + 4*u*(a4 + 4*nu) + 12*u2*(a4 + nu^2)))...
        + (a4 + 8*(-2 + nu))*(8 - 2*nu + u*(a4 + 4*nu) + 2*u2*(a4 + 4*nu) + 4*u3*(a4 + nu^2)))...
        ./(8 - 2*nu + u*(a4 + 4*nu) + 2*u2*(a4 + 4*nu) + 4*u3*(a4 + nu^2)).^2;
        %}
        
        if (comp2nder)
                        
            d2A_u = (-48*nu*u.*(18432*nu^4*u3.*(384 + (-3776 + 123*pi2)*u) + 1179648*(-96 + (-3008 + 123*pi2)*u)...
                + 512*nu*(165888 - 35424*(-64 + 3*pi2)*u - 4*(10040320 - 781296*pi2 + 15129*pi2^2)*u2 - 123*(192512 - 16896*pi2 + 369*pi2^2)*u3 + 6*(3008 - 123*pi2)^2*u4)...
                - 16*nu^3*(-110592 - 576*(-4160 + 123*pi2)*u + (-23228416 + 1259520*pi2 - 15129*pi2^2)*u2 + 6*(8482816 - 692736*pi2 + 15129*pi2^2).*u3 ...
                + 24*(11284480 - 834432*pi2 + 15129*pi2^2)*u4) + nu^2*(-21233664 + 84934656*u + 64*(57094144 - 4990848*pi2 + 105903*pi2^2)*u2 ...
                + (-32852148224 + 3713568768*pi2 - 142333632*pi2^2 + 1860867*pi2^3)*u3 + (-66556788736 + 7741513728*pi2 - 296286336*pi4 + 3721734*pi6)*u4)))...
                ./(768 + 384*nu^2*u3 + nu*(-192 + (3392 - 123*pi2)*u + (6784 - 246*pi2)*u2 + (12032 - 492*pi2).*u3)).^3;
            
            d2A = -2./r.*dA + u4.*d2A_u;
            
        end        
        
        
    case '4pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A   = (1 - (4*(768 - 3584*nu - 24*a5*nu + 123*nu*pi2)*u)/(1536 - 3776*nu + 123*nu*pi2))./ ...
            (1 - (2*(-3392*nu - 48*a5*nu + 123*nu*pi2)*u)./(1536 - 3776*nu + 123*nu*pi2) - ...
            (4*(-3392*nu - 48*a5*nu + 123*nu*pi2)*u2)./(1536 - 3776*nu + 123*nu*pi2) - ...
            (2*(-12032*nu - 192*a5*nu - 3776*nu2 + 492*nu*pi2 + 123*nu2*pi2)*u3)/(1536 - 3776*nu + 123*nu*pi2) + ...
            ((73728*a5*nu + 11505664*nu2 - 18432*a5*nu2 - 834432*nu2*pi2 + 15129*nu2*pi4)*u4)/(96.*(1536 - 3776*nu + 123*nu*pi2)));
        
        
        dA_u  = (384*(-113246208 + nu^3*u2.*(1860867*pi6*u.*(-1 + 3*u) - 363096*pi4*(-6 - 424*u + 3*(432 + a5)*u2) + ...
            377856*pi2*(-354 + (-11164 + 9*a5)*u + 3*(11660 + 47*a5)*u2) + 16384*(125316 - 236*(-9892 + 27*a5)*u + (-7550592 - 38466*a5 + 81*a5^2)*u2)) - ...
            147456*nu*(123*pi2*(1 - 2*u).^2.*(1 + 2*u) - 32*(118 - (212 + 3*a5)*u - 2*(176 + 3*a5)*u2 + (752 - 12*a5)*u3 + 36*a5*u4)) + ...
            48*nu2*(15129*pi4*(-1 + 4*u + 4*u2 - 64*u3 + 48*u4) + ...
            7872*pi2*(118 - (448 + 3*a5)*u - 328*u2 + 4*(1648 + 3*a5)*u3 + 24*(-212 + 3*a5)*u4) - ...
            1024*(13924 - 236*(212 + 3*a5)*u + (-21136 - 144*a5 + 9*a5^2)*u2 + 4*(168448 + 348*a5 + 9*a5^2)*u3 + 12*(-44944 + 1416*a5 + 9*a5^2)*u4))))./ ...
            (147456 + nu2*u3.*(15129*pi4*u - 7872*pi2*(3 + 106*u) - 2048*(-354 - 5618*u + 9*a5*u)) - ...
            96*nu*(123*pi2*(-1 + 2*u + 4*u2 + 8*u3) - 32*(-118 + (212 + 3*a5)*u + (424 + 6*a5)*u2 + 4*(188 + 3*a5)*u3 + 24*a5*u4))).^2;
        
        dA  = -u2.*dA_u;
        
        %{
    dA  = -u2.*(384*(-113246208 + nu^3*u2.*(1860867*pi6*u.*(-1 + 3*u) - 363096*pi4*(-6 - 424*u + 3*(432 + a5)*u2) + ...
        377856*pi2*(-354 + (-11164 + 9*a5)*u + 3*(11660 + 47*a5)*u2) + 16384*(125316 - 236*(-9892 + 27*a5)*u + (-7550592 - 38466*a5 + 81*a5^2)*u2)) - ...
        147456*nu*(123*pi2*(1 - 2*u).^2.*(1 + 2*u) - 32*(118 - (212 + 3*a5)*u - 2*(176 + 3*a5)*u2 + (752 - 12*a5)*u3 + 36*a5*u4)) + ...
        48*nu2*(15129*pi4*(-1 + 4*u + 4*u2 - 64*u3 + 48*u4) + ...
        7872*pi2*(118 - (448 + 3*a5)*u - 328*u2 + 4*(1648 + 3*a5)*u3 + 24*(-212 + 3*a5)*u4) - ...
        1024*(13924 - 236*(212 + 3*a5)*u + (-21136 - 144*a5 + 9*a5^2)*u2 + 4*(168448 + 348*a5 + 9*a5^2)*u3 + 12*(-44944 + 1416*a5 + 9*a5^2)*u4))))./ ...
        (147456 + nu2*u3.*(15129*pi4*u - 7872*pi2*(3 + 106*u) - 2048*(-354 - 5618*u + 9*a5*u)) - ...
        96*nu*(123*pi2*(-1 + 2*u + 4*u2 + 8*u3) - 32*(-118 + (212 + 3*a5)*u + (424 + 6*a5)*u2 + 4*(188 + 3*a5)*u3 + 24*a5*u4))).^2;
        %}
        
        %D = ...
        %B = ...
        %dB = ...
               
        if (comp2nder)                        
            
            d2A_u = (384*(3*nu*(1860867*nu2*pi^6*u2.*(-1 + 4*u) + 294912*a5^2*nu*u.*(-1 - 6*u + 6*(-4 + nu)*u2) - ...
                131072*(2544 - 6254*nu + (8448 - 5284*nu - 10443*nu2)*u + (-27072 + 252672*nu - 291814*nu2)*u2 + 89888*nu*(-3 + 14*nu)*u3) - ...
                484128*nu*pi4*(-2 - (4 + 3*nu)*u + (96 - 318*nu)*u2 + 48*(-2 + 27*nu)*u3) + ...
                503808*pi2*(24 - 112*nu + (96 - 164*nu - 177*nu2)*u - 3*(96 - 1648*nu + 2791*nu2)*u2 + 636*nu*(-8 + 55*nu)*u3) - ...
                96*a5*(32*(1536 + nu*(-3776 + 123*pi2)) - 49152*(-4 + nu)*u - 96*(-6144 + 3*nu2*(-3776 + 123*pi2) + nu*(-7424 + 492*pi2))*u2 + ...
                (-2359296 - 3072*nu*(-3776 + 123*pi2) + nu2*(8753152 - 739968*pi2 + 15129*pi4))*u3)).* ...
                (147456 + nu2*u3.*(15129*pi4*u - 7872*pi2*(3 + 106*u) - 2048*(-354 - 5618*u + 9*a5*u)) - ...
                96*nu*(123*pi2*(-1 + 2*u + 4*u2 + 8*u3) - 32*(-118 + (212 + 3*a5)*u + (424 + 6*a5)*u2 + 4*(188 + 3*a5)*u3 + 24*a5*u4))) - ...
                8*nu*(15129*nu*pi4*u3 - 2304*a5*(-1 - 4*u - 12*u2 + 8*(-4 + nu)*u3) - 1968*pi2*(3 + 12*u + 9*(4 + nu)*u2 + 424*nu*u3) + ...
                1024*(159 + 636*u + 9*(188 + 59*nu)*u2 + 11236*nu*u3)).* ...
                (-113246208 + nu3*u2.*(1860867*pi6*u.*(-1 + 3*u) - 363096*pi4*(-6 - 424*u + 3*(432 + a5)*u2) + ...
                377856*pi2*(-354 + (-11164 + 9*a5)*u + 3*(11660 + 47*a5)*u2) + 16384*(125316 - 236*(-9892 + 27*a5)*u + (-7550592 - 38466*a5 + 81*a5^2)*u2)) ...
                - 147456*nu*(123*pi2*(1 - 2*u).^2.*(1 + 2*u) - 32*(118 - (212 + 3*a5)*u - 2*(176 + 3*a5)*u2 + (752 - 12*a5)*u3 + 36*a5*u4)) + ...
                48*nu2*(15129*pi4*(-1 + 4*u + 4*u2 - 64*u3 + 48*u4) + ...
                7872*pi2*(118 - (448 + 3*a5)*u - 328*u2 + 4*(1648 + 3*a5)*u3 + 24*(-212 + 3*a5)*u4) - ...
                1024*(13924 - 236*(212 + 3*a5)*u + (-21136 - 144*a5 + 9*a5^2)*u2 + 4*(168448 + 348*a5 + 9*a5^2)*u3 + 12*(-44944 + 1416*a5 + 9*a5^2)*u4))) ...
                ))./(147456 + nu2*u3.*(15129*pi4*u - 7872*pi2*(3 + 106*u) - 2048*(-354 - 5618*u + 9*a5*u)) - ...
                96*nu*(123*pi2*(-1 + 2*u + 4*u2 + 8*u3) - 32*(-118 + (212 + 3*a5)*u + (424 + 6*a5)*u2 + 4*(188 + 3*a5)*u3 + 24*a5*u4))).^3;
                        
            d2A = -2./r.*dA + u4.*d2A_u;
            
        end
                        
    case '5pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        c0 = (-3*(512 - 3520*nu - 32*a5*nu - 8*a6*nu + 32*nu^2 + 123*nu*pi^2))/(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2);
        c1 = (3392*nu + 48*a5*nu + 24*a6*nu - 96*nu^2 - 123*nu*pi^2)/(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2);
        c2 = (-2*(-3392*nu - 48*a5*nu - 24*a6*nu + 96*nu^2 + 123*nu*pi^2))/(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2);
        c3 = (-2*(-6016*nu - 96*a5*nu - 48*a6*nu - 3392*nu^2 - 24*a5*nu^2 + 246*nu*pi^2 + 123*nu^2*pi^2))/(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2);
        c4 = (36864*a5*nu + 18432*a6*nu + 11431936*nu^2 + 72192*a5*nu^2 - 4608*a6*nu^2 + 18432*nu^3 - 834432*nu^2*pi^2 - 2952*a5*nu^2*pi^2 + 15129*nu^2*pi^4)...
            /(96.*(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2));
        c5 = (36864*a6*nu + 11358208*nu^2 + 325632*a5*nu^2 + 2304*a5^2*nu^2 - 90624*a6*nu^2 + 362496*nu^3 ...
            - 834432*nu^2*pi^2 - 11808*a5*nu^2*pi^2 + 2952*a6*nu^2*pi^2 - 11808*nu^3*pi^2 + 15129*nu^2*pi^4)...
            /(96.*(768 - 3584*nu - 24*a5*nu + 123*nu*pi^2));
        
        NA   = 1 + c0*u;
        dNA  = c0;
        
        DA   = 1 + c1*u + c2*u2 + c3*u3 + c4*u4 + c5*u5;
        dDA  = c1 + 2*c2*u + 3*c3*u2 + 4*c4*u3 + 5*c5*u4;                
        DA2  = DA.^2;
                
        A     = NA./DA;
        dA_u  = dNA./DA - NA.*dDA./DA2;                
        dA    = -u2.*dA_u;
        
        if (comp2nder)
            
            DA3  = DA2.*DA;                
            d2DA = 2*c2 + 6*c3*u + 12*c4*u2 + 20*c5*u3;          
            
            d2A_u = -2*dNA.*dDA./DA2 - NA.*d2DA./DA2 + 2*NA.*dDA.^2./DA3;    
            d2A   = u4.*d2A_u + 2*u3.*dA_u;
                                                
        end

        
        
    case '5pnlog'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        logu = log(u);
        
        % 4PN and 5PN coefficients including all known log terms
        a5tot  = a5  + 64/5*logu;
        a6tot  = a6  + (-7004/105 - 144/5*nu).*logu;
        a5tot2 = a5tot.^2;      
        
        % Coefficients of the Padeed function        
        N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*pi2)))...
             ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
        D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))...
             ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
        D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))...
             ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
        D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*pi2 - 123*nu*pi2))...
             ./(-768 + nu*(3584 + 24*a5tot - 123*pi2));
        D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*pi2)) + nu*(2048*(5582 + 9*nu) - 834432*pi2 + 15129*pi4)))...
             ./(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
        D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-3008 - 96*nu + 123*pi2))))...
            ./(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
        
        % First derivatives        
        dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
            ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*u);
        dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
            ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*u);
        dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
            ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*u);
        dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
            ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*u);
        dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))...
            ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*u);
        dD5 = (nu*(-8400*nu*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.).*(1536 + nu*(-3776 + 123*pi2)) ...
            + nu*(-2304*power(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.).*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) ...
            - (1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2))).*(4128768*logu*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*pi2) - 756*nu*(1536 + nu*(-3776 + 123*pi2))))))...
            ./(2625.*power(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),2).*u);
        
        % Numerator and denominato of the Pade
        Num = 1 + N1.*u;
        Den = 1 + D1.*u + D2.*u2 + D3.*u3 + D4.*u4 + D5.*u5;
        A   = Num./Den;
                
        % First derivative
        dNum  = dN1.*u + N1;
        dDen  = D1 + u.*(dD1 + 2*D2) + u2.*(dD2 + 3*D3) + u3.*(dD3 + 4*D4) + u4.*(dD4 + 5*D5) + dD5.*u5;
                
        % Derivative of A function with respect to u
        prefactor = A./(Num.*Den);
        dA_u     = prefactor.*(dNum.*Den - dDen.*Num);
        
        % Derivative of A with respect to r
        dA    = -u2.*dA_u;
        
        if (comp2nder)
            
            % Second derivatives of Pade coefficients            
            d2N1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu)...
                + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))...
                ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3).*u2);
            d2D1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) ...
                + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))...
                ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3).*u2);
            d2D2 = (320*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) ...
                + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))...
                ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3).*u2);
            d2D3 = (640*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) ...
                + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))...
                ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3).*u2);
            d2D4 = (320*(-4 + nu)*nu*(-828672 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) ...
                + nu*(5006848 + 42024*a5 + 8064*a6 - 32256*nu - 174045*pi2)).*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2)))...
                ./(7.*power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3).*u2);
            d2D5 = (nu*(power(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2).*(4128768*logu*nu - 7680*(1751 + 756*nu) + nu*(64*(808193 ...
                + 5040*a5 + 223020*nu) - 615*(3095 + 756*nu)*pi2)) + 3072*nu*(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))...
                .*(4128768*logu*nu - 7680*(1751 + 756*nu) + 5*nu*(64*(174541 + 1008*a5 + 44604*nu) - 123*(3095 + 756*nu)*pi2)) ...
                + 25804800*nu2*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.).*(1536 + nu*(-3776 + 123*pi2)) ...
                + nu*(-2304*power(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.).*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) ...
                + 42000*nu*(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2)).*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) ...
                + nu*(-2304*power(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.).*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2)))))...
                ./(13125.*power(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),3).*u2);
            
            % Second derivative of numerator and denominator
            d2Num = 2*dN1 + d2N1.*u;
            d2Den = 2*(D2 + dD1) + u.*(6*D3 + 4*dD2 + d2D1) + u2.*(12*D4 + 6*dD3 + d2D2) ...
                + u3.*(20*D5 + 8*dD4 + d2D3) + u4.*(10.*dD5 + d2D4) + u5.*d2D5;
            
            % Second derivative with respect of u
            d2A_u    = prefactor.*(2.*dDen.^2.*A - 2.*dNum.*dDen + Den.*d2Num - d2Den.*Num);
            
            % second derivative with respect of r
            d2A = u4.*d2A_u + 2*u3.*dA_u;
            
        end
        
    otherwise
        error('unknown option %s for PNorder',PNorder);
        
end


% Compute B and D functions
switch PNorder
        
    case '1pn'

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %D  = ooo;
        B   = 1./A;
        
        %dB = ...
        
    case '2pn'
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %D  = ooo;
        B   =  1./A.*(1.d0-6.d0*nu.*u2);
        
        %dB  =  ...
        
    otherwise
        
        if strcmp(EOBopt.resumD,'pade03')
            
            % Resummation of the D function. One uses the Pade P03[Taylor[D]] to
            % avoid the catastrophic zero of the standard Taylor-expanded
            Dp  = 1.0 + 6*nu*u2 - 2*(3.0*nu-26.0)*nu*u3;
            D   = 1./Dp;
            dD  = 6*u2.*(2*nu*u-(3*nu-26)*nu*u2).*D.^2;
            
        else
            
            % Taylor-expanded D (4PN)
            Dp = 1.0 + 6*nu*u2 - 2*nu*(3.0*nu-26.0)*u3 + nu*(226 + 592/15*log(u)).*u4;
            D  = 1./Dp;
            dD = 2*nu/15*u3.*(90 - 45*u*(3*nu-26) + u2.*(7076+1184*log(u))).*D.^2;
            
            %D   = 1.0 - 6*nu*u2 + 2*(3.0*nu-26.0)*nu*u3;
            %dD  = 12*nu*u3 - 6*(3*nu-26)*nu*u4;
            
        end
        
        B   = D./A;
        dB  = (dD.*A - D.*dA)./A.^2;
        
end


% Finalize
varargout = SetVOutput( nargout, A,dA,d2A, B,dB, dA_u,d2A_u );