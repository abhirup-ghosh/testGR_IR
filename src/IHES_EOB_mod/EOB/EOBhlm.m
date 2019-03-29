function varargout = EOBhlm(nu, t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, EOBmet)

%EOBhlm Compute the multipolar resummed waveform.
%
%   [hlm,phi,omglm,domglm,d2omglm, psilm,Alm,dpsilm] = ...
%   EOBhlm(nu,t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, EOBmet) 
%   return (complex) multipolar wave, phase, frequency, and derivatives,
%   and the RWZ normalized (complex) wave, its amplitude and (complex)
%   derivative.
%
%   WAVE = EOBhlm( .... ) return a structure with the wave
%


particle_dynamics = strcmp(EOBopt.Dynamics,'particle');
rholmarg          = EOBopt.rholmarg;
newton_waves      = strcmp(EOBopt.NewtonWaves,'yes');

L    = EOBopt.L;
M    = EOBopt.M;
kmax = length(L);

nx   = length(Omega);


% Shorthands
pph2 = pph.^2;
r2   = r.^2;
u    = 1./r;
u2   = u.^2;

A    = EOBmet.A;
dA   = EOBmet.dA;


% Compute argument
if particle_dynamics
    rw        = r;
    vphi      = rw.*Omega;
    nu_in_wav = 0;        
else
    W         = A.*(1 + pph2.*u2);
    psi       = 2*(1 + 2*nu*(sqrt(W) - 1.))./(r2.*dA);
    rw        = r.*psi.^(1./3.);
    vphi      = rw.*Omega;
    nu_in_wav = nu;
end

switch rholmarg
    case 'v_omg'
        x = Omega.^(2/3);
    case 'v_phi'
        x = (Omega.*rw).^2;
    case '1_r'
        x = u;
    otherwise
        error('unknown option %s for rholmarg',rholmarg)
end



% Newtonian waveform
hlmNewt = EOBhlmNewt(rw,Omega,phi, nu_in_wav, EOBopt);

if newton_waves
  hlm   = hlmNewt;
  nell  = (sqrt((L+2).*(L+1).*L.*(L-1)))';
  nell  = nell * ones(1,nx);
  psilm = hlm./nell;
  [Alm,philm, omglm,domglm,d2omglm, dpsilm] = EOBExtractAPO(psilm,t);
  varargout = SetVOutput( nargout, hlm,philm,omglm,domglm,d2omglm, psilm,Alm,dpsilm );
  return
end


% Compute corrections
jhat   = pph./(rw.*vphi);

% flm = (rho_lm)^ell
if particle_dynamics    
    flm = EOBflm0(x,EOBopt);    
else    
    flm = EOBflm(x,nu_in_wav,EOBopt);    
    % Further resummation of rho_22 or f_22 (if needed)
    if strcmp(EOBopt.resumf22,'pade23')
        flm(LM2K(2,2),:) = EOBrho22Pade23(xarg,nu).^2;
    elseif strcmp(EOBopt.resumf22,'oldpade')
        flm(LM2K(2,2),:) = EOBf22Pade(xarg,nu);
    end    
end

% Tail
r0  = 1.213061319425267e+00;   % 2/sqrt(e);
tlm = EOBhhatlmTail( L,M, Omega,E, r0);

% Prefactor
Heff = reshape(Heff,1,nx);
jhat = reshape(jhat,1,nx);
lmeven  = (~mod(L+M,2))';
lmodd   = ( mod(L+M,2))';
prefact = lmeven * Heff + lmodd * jhat;

% Compute \hat{h}_lm
hhatlm = prefact .* flm.*tlm;

% Residual phase corrections delta_{lm}
deltalm    = EOBdeltalm(E,Omega, nu_in_wav, EOBopt);


% Complete h_{lm}
hlm   = hlmNewt.* hhatlm .* exp(1i*deltalm);


% Compute phase, amplitude and frequency and RWZ normalized waveform
nell  = (sqrt((L+2).*(L+1).*L.*(L-1)))';
nell  = nell * ones(1,nx);
psilm = hlm./nell;

[Alm,philm, omglm,domglm,d2omglm, dpsilm] = EOBExtractAPO(psilm,t);


% Finalize
varargout = SetVOutput( nargout, hlm,philm,omglm,domglm,d2omglm, psilm,Alm,dpsilm);

