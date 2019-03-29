function varargout = T4Run(q,r0,Phi0,dt,tmax, varargin)

%T4Run Taylor-T4 approximant.
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax)
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, ODEoptions) specify ODExx optional
%   structure (see odeset.m)
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],pnorder) specify PN order
%   ('2.5PN' or '3PN')  
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],[], addelta) re-define waveform
%   phase to absorb the log(x/x0) (if addelta='yes') 
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],[],[],computeAPO) compute
%   modulus, phase, frequency and curvature waveform (if computeAPO='yes')
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],[],[],[], verbose) screen info
%   level  
%             0 = minimum text (default)
%             1 = text
%             2 = text + figures 
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],[],[],[],[], saveme) save data, 
%   0 = no 1 = yes (default) 
%
%   [dyn,wav] = T4Run(q,r0,Phi0,dt,tmax, [],[],[],[],[],[], basedir ) specify
%   output directory, default = 'mydata'.
% 
%   DATA = T4Run(q,r0,Phi0,dt,tmax,..) return structure output
%
%   Reference(s)
%    Boyle et al, XXX YYY (ZZZ)
%    Kidder, PRD 77, 044016 (2008).
%


fprintf('===> T4 run\n');
tottime = tic;


% Manage args in
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'Refine',1);
pnorder = '3pn';
addelta = 'no';
compAPO = 'no';
verbose = 0;   % do not print info
saveme  = 1;   % do save results
basedir = 'mydata/';

optargs = {options pnorder addelta compAPO verbose saveme basedir};
na = length(varargin);
if (na>7)
    error('too many input args')
end
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[options, pnorder, addelta, compAPO, verbose, saveme, basedir] = optargs{:};


% Init run
fprintf(' ===> Init run\n');
tstart = tic;

if verbose<0 | verbose >2
    error('unknown value for ''verbose'' par, {0,1,2} => no, screen info, screen info and fig');
end
fhIdx = 0; % figures index


% Create output dirs (if necessary)
if saveme
    
    basedir = strcat(basedir,'/');
    [ok,mess,messid] = mkdir(basedir);       if ~ok, error(mess); end;
    
    fid = fopen(strcat(basedir,'log.txt'),'w');
    fprintf(fid,'%s\n Notes:\n',datestr(now, 'mm/dd/yyyy HH:MM:SS'));
    fclose(fid);
    
    if verbose==2
        [ok,mess,messid] = mkdir(basedir,'Fig'); if ~ok, error(mess); end;
    end
    
end


% Initial data
x0  = 1/r0; 
nu  = q/(1+q)^2; 
tspan = [0:dt:tmax];
y0 = [x0; Phi0];

fprintf(' Done %4.3f sec\n',toc(tstart));

if verbose
          
    fprintf(' nu           = %g\n',nu);
    fprintf(' q            = %g\n',q);
    fprintf(' PNorder      = %s\n',pnorder);
    
    fprintf(' dt           = %+.12e\n',dt);
    fprintf(' tmax         = %+.12e\n',tmax);    
    fprintf(' ODESolRelTol = %+.12e\n',options.RelTol);
    fprintf(' ODESolAbsTol = %+.12e\n',options.AbsTol);
    fprintf(' ODESolRefine = %d\n',options.Refine);
      
    fprintf(' r            = %+.12e\n',y0(1));    
    
end


% Dynamics
fprintf(' ===> Dynamics ...\n');
tstart = tic;
[t,y]  = ode45(@(t,y) T4Rhs(t,y,nu),tspan,y0,options);
fprintf(' Done %4.3f sec\n',toc(tstart));

y = y.';
t = t.';


% Compute Omega
fprintf(' ===> Compute orbital frequency ...\n');
tstart = tic;
dydt  = T4Rhs(t,y,nu);
Omega = dydt(2,:); 
fprintf(' Done %4.3f sec\n',toc(tstart));

% Pack dynamics
dyn.t     = t;
dyn.y     = y;
dyn.dydt  = dydt;
dyn.Omega = Omega;


% Wave
fprintf(' ===> Waveform ...\n');
tstart = tic;
wav = T4Wav(nu,t,y(1,:),y(2,:), pnorder, addelta);
fprintf(' Done %4.3f sec\n',toc(tstart));

wav.t = t;

if verbose
    
    % todo
    fprintf(' # of orbits     = %g\n',y(2,end)/(2*pi));
end


% Compute modulus, phase, frequency and curvature wave (if required)
if strcmp(lower(compAPO),'yes')
            
    coeff      =  12/5*2^(1/3)*nu;
            
    % Metric waveform
    wav.phi22  = - unwrap(angle(wav.psi22));
    wav.omg22  = FDdrvt(wav.phi22,t,4); 
    wav.domg22 = FDdrvt(wav.omg22,t,4);
    wav.A_omg  = wav.domg22./( coeff*wav.omg22.^(11/3) );
    wav.Q_omg  = wav.omg22.^2./( wav.domg22 );
    
    % curvature waveform
    [tmp,Psi4_22] = FDdrvt(wav.h22,t,4);    
    ModPsi4_22    = abs(Psi4_22);
            
    phi    = -unwrap(angle(Psi4_22));
    omg22  = FDdrvt(phi,t,4);
    domg22 = FDdrvt(omg22,t,4);        
    
    curv.Psi4_22 = Psi4_22;
    curv.ModPsi4_22 = ModPsi4_22;
    curv.omega22 = omg22;
    curv.A_omega = domg22./(coeff*omg22.^(11/3));    
    curv.Q_omega = omg22.^2./domg22;

    wav.curv = curv;
    
end


if saveme
    
    fprintf(' ===> Save *.mat files ...\n');
    tstart = tic;        
        
    filename  = strcat(basedir,'Dynamics');
    save(filename,'dyn');
    [err mess] = system(sprintf('du -h %s.mat ',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);    

    filename  = strcat(basedir,'Waveform');
    save(filename,'wav');
    [err mess] = system(sprintf('du -h %s.mat ',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);    
    
    fprintf(' Done %4.3f sec\n',toc(tstart));
    
end


% Output
varargout = SetVOutput( nargout, dyn,wav );
fprintf('End %4.3f sec\n',toc(tottime));

