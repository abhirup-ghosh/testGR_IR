function deltalm = EOBdeltalm(Hreal,Omega,nu, EOBopt)


%EOBdeltalm Residual phase corrections delta_{lm} up to l=m=5.
%
%   EOBdeltalm(Hreal,Omega,nu, EOBopt);
%
% 
%   nu=0, the delta_lmm are written in Taylor-expanded form and all terms
%   up to 4.5PN accuracy are included. 
%   nu~=0, 
%         (i) the 4.5PN, test-mass terms are excluded because of the too
%         large PN-gap with the nu-dependent terms 
%         (ii) the deltalmm are replaced by suitable Pade' approximants for
%         multipoles (2,2), (2,1), (3,3), (3,1)
%   The l=m=2 residual phase includes the 3.5PN (nu-dependent) correction
%   obtained in Faye et al. 
%
%   Reference(s)
%    Damour, Iyer & Nagar, PRD 79, 064004 (2008) 
%    Fujita & Iyer, PRD 82 044051 (2010)
%    Faye et al., Class. Q. Grav. 29 175004 (2012)
%    Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013)
%


% TODO: this routine requires optimization
% - precompute coefficients c(nu)
% - evaluate efficiently polynomials

L         = EOBopt.L;
LM2K      = EOBopt.LM2K;


% Shorthands
pi2    = pi.^2;
nu2    = nu^2;

y      = (Hreal.*Omega).^(2/3);
sqrt_y = sqrt(y);
y3     = y.^3;
y32    = Hreal.*Omega;
y92    = y.^(9/2);



% Compute deltalm
kmax  = length(L);
deltalm = zeros(kmax,length(y));

if nu==0
        
    % Residual phases in Taylor-expanded form    
    
    % l=2 ---------------------------------------------------------------
    
    deltalm(LM2K(2,2),:) =   7d0/3d0*y32 + 428d0/105d0*pi * y3  + (-2203/81 + 1712/315*pi2) * y92;
    deltalm(LM2K(2,1),:) =   2d0/3d0*y32 + 107d0/105d0*pi * y3  + (-272/81  +  214/315*pi2) * y92;
    
    % l=3 ---------------------------------------------------------------
    
    deltalm(LM2K(3,3),:) = 13d0/10d0 * y32  + 39/7 *pi * y3 + (-227827/3000  + 78/7  *pi2) * y92;
    deltalm(LM2K(3,2),:) = 2d0/2     * y32  + 52/21*pi * y3 + (-9112/405     + 208/63*pi2) * y92;
    deltalm(LM2K(3,1),:) = 13d0/30d0 * y32  + 13/21*pi * y3 + (-227827/81000 + 26/63 *pi2) * y92;
        
    % l=4 ---------------------------------------------------------------

    deltalm(LM2K(4,4),:) = 112d0/120d0  * y32 + 25136/3465*pi*y3 + (-55144/375 + 201088/10395*pi2)*y92;
    deltalm(LM2K(4,3),:) = 486d0/810d0  * y32 + 1571/385*pi*y3;
    deltalm(LM2K(4,2),:) =   7d0/15d0   * y32 + 6284/3465*pi*y3  + ( -6893/375 +  25136/10395*pi2)*y92;
    deltalm(LM2K(4,1),:) =   2d0/10d0   * y32 + 1571/3465*pi*y3;
    
    % l=5 ---------------------------------------------------------------
    
    deltalm(LM2K(5,5),:) = 96875/131250 * y32;
    
    return
    
end


% Residual phases in Pade-resummed form when possible

% Leading order contributions
delta22LO = 7d0/3d0   * y32;
delta21LO = 2d0/3d0   * y32;
delta33LO = 13d0/10d0 * y32;
delta31LO = 13d0/30d0 * y32;

% l=2 ------------------------------------------------------------------

% Pade(2,2) approximant
num        = (808920*nu*pi*sqrt(y) + 137388*pi2*y + 35*nu2*(136080 + (154975 - 1359276*nu)*y));
den        = (808920*nu*pi*sqrt(y) + 137388*pi2*y + 35*nu2*(136080 + (154975 + 40404*nu)*y));
deltalm(LM2K(2,2),:) = delta22LO.*num./den;

% Pade(1,2) approximant
num        = 69020*nu + 5992*pi*sqrt_y;
den        = 5992*pi*sqrt_y + 2456*nu*(28+493*nu*y);
deltalm(LM2K(2,1),:) = delta21LO.*num./den;

% l=3 ------------------------------------------------------------------

% Pade(1,2) approximant
num        = 1   + 94770*pi/(566279*nu)*sqrt_y;
den        = num + 80897*nu/3159*y;
deltalm(LM2K(3,3),:) = delta33LO.*num./den;

% Taylor-expanded form
deltalm(LM2K(3,2),:) = (10d0+33d0*nu)/(15*(1-3*nu))*y32    +  52/21*pi*y3;

% Pade(1,2) approximant
num        = 4641*nu + 1690*pi*sqrt_y;
den        = num + 18207*nu2*y;
deltalm(LM2K(3,1),:) = delta31LO.*num./den;

% l=4 ------------------------------------------------------------------

deltalm(LM2K(4,4),:) =  (112d0+219*nu)/(120d0*(1-3d0*nu)) * y32 + 25136/3465*pi*y3;
deltalm(LM2K(4,3),:) = (486d0+4961*nu)/(810*(1d0-2d0*nu)) * y32 +   1571/385*pi*y3;
deltalm(LM2K(4,2),:) =  7d0*(1+6d0*nu)/(15d0*(1-3d0*nu))  * y32 +  6284/3465*pi*y3;
deltalm(LM2K(4,1),:) =   (2d0+507*nu)/(10d0*(1-2d0*nu))   * y32 +  1571/3465*pi*y3;

% l=5 ------------------------------------------------------------------

deltalm(LM2K(5,5),:) = (96875 + 857528*nu)/(131250*(1-2d0*nu)) * y32;

