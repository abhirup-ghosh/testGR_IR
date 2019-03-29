function [Fphi, Frstar] = EOBRRNewt(nu, dydt, y, EOBopt, EOBmet, EOBHam)

%EOBRRno Compute radiation reaction for generic orbits at Newtonian order.
%
%   [Fphi, Frstar] = EOBRRNewt(nu, dydt, y, EOBopt, EOBmet, EOBHam)
%
%

%FIXME ref in doc :   Iyer-Will 2005
%FIXME Frstar 

% take these from EOBopt
ab = [-1 0]; % Harmonic gauge
%ab = [4 5]; % ADM gauge


alpha = ab(1);
beta  = ab(2);



% Unpack y
phi    = y(1);
r      = y(2);
pph    = y(3);
prstar = y(4);


% Shorhands
r2       = r.^2;
r3       = r.^3;
u        = 1./r;
u2       = u.^2;
u3       = u.^3;
u4       = u.^4;
u5       = u.^5;
pph2     = pph.^2;
prstar2  = prstar.^2;

Fphi = -nu*8/5*pph*( -u3.*prstar2.*(1+2*alpha) + ...
    (2-alpha)*u4 + pph2*(2+alpha)*u5 );

Fr   = nu*8/15*( u3.*prstar3.*( 2 + 6*(alpha-beta)) ...
    + prstar .*( u4.*(17+9*(alpha-beta)) + 3*u5.*pph2.*(1+3*(3*beta-alpha)) ) );

warning(' Frstar need to be fixed in this routine')
Frstar   = Fr; %FIXME!