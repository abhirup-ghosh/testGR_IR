function varargout = T4Wav(nu,t,x,phi, varargin)

%T4Wav Taylor-T4 waveform at 2.5 or 3.0 PN accuracy.
%
%   [h22, psi22] = T4Wav(nu,t,x,phi)
%
%   [h22, psi22] = T4Wav(nu,t,x,phi, pnorder) specify PN order ('2.5PN' or
%   '3PN') 
%
%   [h22, psi22] = T4Wav(nu,t,x,phi, [],addelta) re-define phase to absorb
%   the log(x/x0) (if addelta='yes') 
%  
%   wav = T4Wav(nu,t,x,phi) return a structure with the wave
%
%   Reference(s)
%    Boyle et al, XXX YYY Eq. (30),(46).
%    Kidder, PRD 77, 044016 (2008), Eq.(79).
%


% Manage args in
order    = '3pn';
addelta  = 'no';

optargs = {order addelta};
na = length(varargin);
if (na>2)
    error('too many input args')
end
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[order addelta] = optargs{:};


% Shorthands
gamma =  0.5772156649015329;  % Euler constants
pi2   =  pi^2;                  
N22   = -8*nu*sqrt(pi/5.0);      % The normalization
b     =  2.0;
nu2   = nu^2;
nu3   = nu^3;

x2    = x.^2;
x3    = x.^2;
x32   = x.^(3/2);
x52   = x.^(5/2);
logx  = log(x);


% Wave at given order
% TODO: the following can be optimized
switch lower(order)
    
    case '2.5pn'
    
        htilde22 = N22* x.* ( 1-x.*(107/42-55/42*nu) + x32.*2*pi...
            -(2173/1512 + 1069/216*nu - 2047/1512*nu2).*x2...
            -x52.*( (107/21 - 34/21*nu)*pi+24*i*nu)   );    
    
    case '3pn'
    
        htilde22 = N22* x.* ( 1-x.*(107/42-55/42*nu) + x32.*2*pi...
            - (2173/1512 + 1069/216*nu - 2047/1512*nu2).*x2...
            - x52.*( (107/21 - 34/21*nu)*pi+24*i*nu) ...
            +  x3.*( 27027409/646800 - 856/105*gamma + 2/3*pi2 ...
            - 1712/105*log(2) -428/105*logx...
            - (278185/33264-41/96*pi2)*nu -20261/2772*nu2 ...
            + 114635/99792*nu3 + 428/105*i*pi) );
    
    otherwise
        
        error(' only 2.5PN or 3PN allowed');
        
end


% Re-define phase to absorb the log(x/x0)  (if required)
% Eq.(30) of Boyle et al.
if strcmp(lower(addelta),'yes')
    logx0 = 11/18 -2.0/3.0*gamma + 2.0/3.0*log(1.0/(4*b));
    delta  = -3*(1-x/8.0).*x32.*(logx-logx0);
    phi = phi+delta;            
end
   

% Wave
emitp  = exp(-2*i*phi);
h22    = htilde22 .* emitp;
psi22 = h22/sqrt(24.0);


% Output
varargout = SetVOutput( nargout, psi22,h22 );


