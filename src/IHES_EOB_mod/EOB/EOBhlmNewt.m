function hlmNewt = EOBhlmNewt(r,Omega,phi, nu, EOBopt)

%EOB_hlmNewt Computes the leading-order (Newtonian) prefactor of the
%multipolar resummed waveform.  
%
%   hlmNewt = EOB_hlmNewt(r,Omega,phi, nu)
%
%   Reference(s)
%    Damour, Iyer & Nagar, PRD 79, 064004 (2009)
%

% TODO: this routine requires optimization
% - precompute coefficients c(nu)
% - evaluate efficiently polynomials

L         = EOBopt.L;
M         = EOBopt.M;
LM2K      = EOBopt.LM2K;


% Shorthands
nu2   = nu.^2;
nu3   = nu.^3;

vphi  = r.*Omega;
vphi2 = vphi.^2;
vphi3 = vphi.^3;
vphi4 = vphi.^4;
vphi5 = vphi.^5;
vphi6 = vphi.^6;
vphi7 = vphi.^7;
vphi8 = vphi.^8;


% Poly in nu
p1 = 1;
p2 = sqrt(1-4*nu);
p3 = (3*nu-1);
p4 = (2*nu-1)*sqrt(1-4*nu);
p5 = 1-5*nu+5*nu2;
p6 = (1-4*nu+3*nu2)*sqrt(1-4*nu);
p7 = 7*nu3 - 14*nu2 + 7*nu -1;


% Allocate memory
kmax  = length(L);
nx    = length(r);
hlmNewt = zeros(kmax,length(r));
%hlmNewt = sparse(zeros(kmax,length(r)));


% Phase factor
M   = reshape(M,kmax,1);
phi = reshape(phi,1,nx);

expimphi = M*phi;
expimphi = exp(- 1i* expimphi);


% Compute hlmNewt (without phase factor)

% l=2 ------------------------------------------------------------------

hlmNewt(LM2K(2,2),:) = -8.*sqrt(pi/5)         * p1 * vphi2; ... .*exp(-2*1i.*phi);
hlmNewt(LM2K(2,1),:) = -8/3*1i*sqrt(pi/5)     * p2 * vphi3; ... .*exp(-1i.*phi);

% l=3 ------------------------------------------------------------------

hlmNewt(LM2K(3,3),:) = 3*1i*sqrt(6*pi/7)     * p2 * vphi3; ... .*exp(-3*1i.*phi);
hlmNewt(LM2K(3,2),:) = 8/3*sqrt(pi/7)        * p3 * vphi4; ... .*exp(-2*1i.*phi);
hlmNewt(LM2K(3,1),:) = -1/3*1i*sqrt(2*pi/35) * p2 * vphi3; ... .*exp(-1i.*phi);

% l=4 ------------------------------------------------------------------

hlmNewt(LM2K(4,4),:) = -64/9*sqrt(pi/7)       * p3 * vphi4; ... .*exp(-4*1i.*phi);
hlmNewt(LM2K(4,3),:) = -9/5*1i*sqrt(2*pi/7)   * p4 * vphi5; ... .*exp(-3*1i.*phi);
hlmNewt(LM2K(4,2),:) = 8/63*sqrt(pi)          * p3 * vphi4; ... .*exp(-2*1i.*phi);
hlmNewt(LM2K(4,1),:) = 1/105*1i*sqrt(2*pi)    * p4 * vphi5; ... .*exp(-1i.*phi);

% l=5 ------------------------------------------------------------------

hlmNewt(LM2K(5,5),:) = 125/12*1i*sqrt(5*pi/66)* p4 * vphi5; ... .*exp(-5*1i.*phi);
hlmNewt(LM2K(5,4),:) = 256/45*sqrt(pi/33)     * p5 * vphi6; ... .*exp(-4*1i.*phi);
hlmNewt(LM2K(5,3),:) = -9*1i/20*sqrt(3*pi/22) * p4 * vphi5; ... .*exp(-3*1i.*phi);
hlmNewt(LM2K(5,2),:) = -16/135*sqrt(pi/11)    * p5 * vphi6; ... .*exp(-2*1i.*phi);
hlmNewt(LM2K(5,1),:) = 1i/180*sqrt(pi/77)     * p4 * vphi5; ... .*exp(-1i.*phi);

% l=6 ------------------------------------------------------------------

hlmNewt(LM2K(6,6),:) = -432/5*sqrt(pi/715)      * p5 * vphi6; ... .*exp(-6.*1i.*phi);
hlmNewt(LM2K(6,5),:) = -1i*625/63*sqrt(5*pi/429)* p6 * vphi7; ... .*exp(-5.*1i.*phi);
hlmNewt(LM2K(6,4),:) = 1024/495*sqrt(2*pi/195)  * p5 * vphi6; ... .*exp(-4.*1i.*phi);
hlmNewt(LM2K(6,3),:) = 1i*81/385*sqrt(pi/13)    * p6 * vphi7; ... .*exp(-3.*1i.*phi);
hlmNewt(LM2K(6,2),:) = -16/1485*sqrt(pi/13)     * p5 * vphi6; ... .*exp(-2.*1i.*phi);
hlmNewt(LM2K(6,1),:) = -1i/2079*sqrt(2*pi/65)   * p6 * vphi7; ... .*exp(-1i.*phi);

% l=7 ------------------------------------------------------------------

hlmNewt(LM2K(7,7),:) = 16807/180*1i*sqrt(7*pi/4290)  * p6 * vphi7; ... .*exp(-7.*1i.*phi);
hlmNewt(LM2K(7,6),:) = 648/35*sqrt(3*pi/715)         * p7 * vphi8; ... .*exp(-6.*1i.*phi);
hlmNewt(LM2K(7,5),:) = -1i*3125/3276*sqrt(5*pi/66)   * p6 * vphi7; ... .*exp(-5.*1i.*phi);
hlmNewt(LM2K(7,4),:) = -1024/1365*sqrt(2*pi/165)     * p7 * vphi8; ... .*exp(-4.*1i.*phi);
hlmNewt(LM2K(7,3),:) = 1i*243/20020*sqrt(3*pi/10)    * p6 * vphi7; ... .*exp(-3.*1i.*phi);
hlmNewt(LM2K(7,2),:) = 8/3003*sqrt(pi/15)            * p7 * vphi8; ... .*exp(-2.*1i.*phi);
hlmNewt(LM2K(7,1),:) = -1i/108108*sqrt(pi/10)        * p6 * vphi7; ... .*exp(-1i.*phi);

% l=8 ------------------------------------------------------------------

hlmNewt(LM2K(8,8),:) = -131072/315*sqrt(2*pi/17017)  * p7 * vphi8; ... .*exp(-8.*1i.*phi);


% Multiply by phase factor
hlmNewt = hlmNewt .* expimphi;

