function varargout = EOBRun_DIN_mod_find_Mbh_abh( nu, ...
    ainput, ... % a1,a2,a3,a4,a5,a6 (rem: a1...a4 multipolar!)
    r0, tmax, dt, nqnm,lmaxqnm, ...
    varargin)

%EOBRun Run EOB simulation.
%
%   DATA = EOBRun( nu, a, r0, tmax,dt, nqnm,lmaxqnm ) perform EOB
%   run with default options (see EOBSet.m) and parameters
%
%    nu      : symmetric mass ratio (scalar)
%    a       : a{1} ,..., a{4}, a{5}, a{6} NQC parameters for the
%    multipolar waveform and deformation paramters (1x6 cell array
%    of vectors)  
%    r0      : initial value of EOB radius
%    tmax    : final time for dynamics
%    dt      : time spacing output
%    nqnm    : number of positive and negative QNM (vector size=2) 
%    lmaxqnm : max l-value to be used in QNM match (scalar)
%
%   Return a data structure with dynamics and waveform.
%
%   DATA = EOBRun( ... , lmaxqnm, EOBOPTIONS ) specify options from
%   structure EOBOPTIONS (see EOBSet.m)
% 
%   DATA = EOBRun( ... , lmaxqnm, [], modus ) specify modus
%   operandi 
%             0 = dynamics
%             1 = dynamics + inspl wave
%             2 = dynamics + inspl-merg wave 
%             3 = dynamics + inspl-merg-ringdown wave (default)
%
%   DATA = EOBRun( ... , lmaxqnm, [],[], verbose ) screen info
%   level 
%             0 = minimum text (default)
%             1 = text
%             2 = text + figures 
%
%   DATA = EOBRun( ... , lmaxqnm, [],[],[], saveme ) save data 0 =
%   no 1 = yes (default)
%
%   DATA = EOBRun( ... , lmaxqnm, [],[],[],[], basedir ) specify
%   output directory, default = 'mydata'.
%
%   [dyn,wav] = EOBRun( ... ) Return structures with dynamics and
%   waveform seprately. 
%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings

fprintf('===> EOB run\n');
tottime = tic;

% Manage args in
EOBopt  = [];                      % struct of options (EOBSet.m)
modus   = 3;                       %
verbose = 0;                       % do not print info
saveme  = 1;                       % do save results
basedir = 'mydata/';

optargs = {EOBopt modus verbose saveme basedir modus};
na = length(varargin);
if (na>6)
    error('too many input args')
end
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[EOBopt, modus, verbose, saveme, basedir] = optargs{:};

% Set default EOB options
defaultopt = EOBSet_mod(...
    'Dynamics','eob',...
    'RadReac','din_mod',...
    'ddotrMethodDyn','noflux',...
    'ddotrMethodWav','noflux',...
    'HorizonFlux','yes',...
    'RadialFlux','yes',...
    'FIterms','yes',...
    'rholmarg','v_phi',...
    'PNorder','5pnlog',...
    'resumD','pade03',...
    'resumf22','no',...
    'Tidal','no',...
    'PNTidal','nnlo',...
    'NewtonWaves','no',...
    'DetermineNQCab','no',...   
    'fac',1.0,...
    'a1',1.0,...
    'a2',1.0,...
    'NQCFileNR','../DataNR/NRpts4NQC.dat',...
    'QNMDataDir','../DataQNMs/',...
    'textend',500,...
    'ODEEventFunRmin',1.0,...
    'ODESolverRelTol',1e-9,...
    'ODESolverAbsTol',1e-12,...
    'ODESolverRefine',1);

% Overwrite EOB options
EOBopt = EOBSet_mod(defaultopt, EOBopt );

% Lower the string entries in EOBoptions
Name   = fieldnames(EOBopt);
idx    = strcmp(cellfun(@(x) class(EOBopt.(x)),Name,'UniformOutput',0),'char');
EOBopt = cellfun(@(x) setfield(EOBopt,x,lower(EOBopt.(x))),Name(idx));
EOBopt = EOBopt(1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize run

fprintf(' ===> Init run\n');
tstart = tic;

% Compute multipolar indexes
lmax = 8;
[LM2K, L,M] = Multipidx(lmax);

EOBopt.LM2K = LM2K;
EOBopt.L    = L;
EOBopt.M    = M;

kmax = length(L);

% Store dynamics pars
EOBopt.a5  = ainput{5};
EOBopt.a6  = ainput{6};

% Set function handles
switch EOBopt.Dynamics
    case 'eob'
        EOBopt.ComputeEOBHam = @EOBHam;
        EOBopt.ComputeRhs    = @EOBRhs;
        EOBopt.ComputeID     = @EOBIdPostAdiab_mod;                
    case 'particle'
        EOBopt.ComputeEOBHam = @EOBHam0;
        EOBopt.ComputeRhs    = @EOBRhs0;
        EOBopt.ComputeID     = @EOBIdPostAdiab0;        
    otherwise
        error('unknown option %s for Dynamics.',EOBopt.Dynamics);
end

switch EOBopt.RadReac
    case 'din'
        if strcmp(EOBopt.Dynamics,'particle')
            EOBopt.ComputeRRForce = @EOBRRDIN0;
        else
            EOBopt.ComputeRRForce = @EOBRRDIN;
        end        
    case 'din_mod'
            EOBopt.ComputeRRForce = @EOBRRDIN_mod;
    case 'pn'
            EOBopt.ComputeRRForce = @EOBRRPN;
    case 'newt'    
        EOBopt.ComputeRRForce = @EOBRRNewt;
    case 'geo'
        EOBopt.ComputeRRForce = @EOBRRGeo;
    otherwise
        error('unknown option %s for RadReac.',EOBopt.RadReac);
end

switch EOBopt.ddotrMethodDyn
    case 'boot2'
        if strcmp(EOBopt.Dynamics,'particle')
            EOBopt.Computeddotr = @EOBddotrBoot20;            
        else
            EOBopt.Computeddotr = @EOBddotrBoot2;            
        end
    case 'noflux'
        EOBopt.Computeddotr = @EOBddotrNoFlux;
    otherwise
        error('unknown option %s for ddotrMethodDyn.',EOBopt.ddotrMethodDyn);
end

EOBopt.ODESolverEventF = @(t,y) EOBStop(t,y, EOBopt.ODEEventFunRmin);

% Set options for ODE solver and time vector
ODEopt = odeset( 'RelTol',EOBopt.ODESolverRelTol,...
    'AbsTol',EOBopt.ODESolverAbsTol,...
    'Refine',EOBopt.ODESolverRefine,...
    'Events',EOBopt.ODESolverEventF );

% Build time vector
% NOTE: dt given as input is the width of the matching grid. The spacing between each
%       point of the "matching grid" is given by dt/(#QNMs - 1).
%       FIXME: this can be probably done just locally after the integration.

dt    = dt/(sum(nqnm)-1); 
tspan = 0:dt:tmax;

% Init NQC
if strcmp(EOBopt.DetermineNQCab,'yes') 
    % Allocate memory
    a{1} = zeros(1,kmax);
    a{2} = a{1};
    a{3} = a{1};
    b{1} = a{1};
    b{2} = a{1};
    b{3} = a{1};
    % Set some for dynamics
    % a' and b's will be overwritten by the call EOBNQCabFind.m
    a{1}(1:length(ainput{1})) = ainput{1};
    a{2}(1:length(ainput{2})) = ainput{2};
    a{3}(1:length(ainput{3})) = ainput{3};
    a{4}(1:length(ainput{4})) = ainput{4};
    if exist(EOBopt.NQCFileNR,'file')
        NRPTS = importdata(EOBopt.NQCFileNR,' ',3);
        NRPTS = NRPTS.data;
    else
        error('file %s does not exist',EOBopt.NQCFileNR);
    end
else
  % do not determine a's and b's
  % use pre-computed
  if strcmp(EOBopt.Dynamics,'particle')    
    [a,b] = EOBNQCab0(EOBopt);
  else
    [a,b] = EOBNQCabFit(nu, EOBopt);
  end
end

EOBopt.NQC.a = a;
EOBopt.NQC.b = b;

if (strcmp(lower(EOBopt.RadReac),'din') || strcmp(lower(EOBopt.RadReac),'din_mod')) ...
        && strcmp(lower(EOBopt.HorizonFlux),'yes')
    % Init coefficients used in horizon flux
    c = EOBFluxHorizonFitCoefs(nu, EOBopt);
    EOBopt.HorizonFluxCoefs.c = c;
end


if modus<0 || modus >3
    error('unknown value for ''modus'' par, {0,1,2,3}');
end

if verbose<0 || verbose >2    
    error('unknown value for ''verbose'' par, {0,1,2} => no, screen info, screen info and fig');
end
fhIdx = 0; % figures index

if nqnm(2)>nqnm(1) || nqnm(1)<0
    error('nQNM must have two entries, nQNM(1)>=0 and nQNM(2)<=nQNM(1)')
end

% Create output dirs (if necessary)
if saveme
    
    basedir = strcat(basedir,'/');
    [ok,mess,messid] = mkdir(basedir);       
    if ~ok, error(mess); end;
    
    fid = fopen(strcat(basedir,'log.txt'),'w');
    fprintf(fid,'%s\n Notes:\n',datestr(now, 'mm/dd/yyyy HH:MM:SS'));
    fclose(fid);
    
    if verbose==2
        [ok,mess,messid] = mkdir(basedir,'Fig'); if ~ok, error(mess); end;
    end
    
end

fprintf(' Done %4.3f sec\n',toc(tstart));

if verbose
  PrintSFields(EOBopt,1);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial data

fprintf(' ===> Compute initial data for evolution\n');
tstart = tic;
y  = EOBopt.ComputeID(nu, r0, EOBopt );
fprintf(' Done %4.3f sec\n',toc(tstart));

y0      = y(1:4); % this goes into theDyn     ODE call
r0      = y0(2);
pph0    = y0(3);
prstar0 = y0(4);
j0      = y(6);

if verbose
    fprintf(' r      = %+.12e\n',y0(2));
    fprintf(' pph    = %+.12e\n',y0(3));
    fprintf(' prstar = %+.12e\n',y0(4));
    fprintf(' pr     = %+.12e\n',y(5));
    fprintf(' j      = %+.12e\n',y(6)); 
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adiabatic LSO, EOB potential and Hamiltonian (if required)

if verbose
    
    fprintf(' ===> Compute Adiabatic LSO\n');
    tstart = tic;
    LSO = EOBAdiabLSO(nu, EOBopt);
    fprintf(' Done %4.3f sec\n',toc(tstart));    
    
    fprintf(' Omega        [LSO Adiab] = %+.12e\n',LSO.Omega);
    fprintf(' omega22      [LSO Adiab] = %+.12e\n',LSO.w22);
    fprintf(' Omega   Newt [LSO Adiab] = %+.12e\n',LSO.Newt_Omega);
    fprintf(' omega22 Newt [LSO Adiab] = %+.12e\n',LSO.Newt_w22);
    fprintf(' r            [LSO Adiab] = %+.12e\n',LSO.r);
    fprintf(' u            [LSO Adiab] = %+.12e\n',LSO.u);
        
    fprintf(' ===> Compute EOB metric, potential and Hamiltonian\n');
    tstart = tic;
    R          = 1.1:0.01:100;
    Metric     = EOBMetric( nu, R, EOBopt );
    Wadiab     = EOBPotential(nu, R, j0  , EOBopt, Metric );
    W          = EOBPotential(nu, R, pph0, EOBopt, Metric );
    [H0 Heff0] = EOBopt.ComputeEOBHam(nu, r0,pph0,prstar0, EOBopt);
    fprintf(' Done %4.3f sec\n',toc(tstart));    
    
    fprintf(' nu H0 = %+.12e\n',nu*H0);
    fprintf(' Heff0 = %+.12e\n',Heff0);
        
    if verbose==2
        
        fhIdx  = fhIdx+1; fhName{fhIdx} = 'Potential'; fh(fhIdx) = figure;
        plot(R,W,'b-',R,Wadiab,'r--');
        hold on        
        if strcmp(EOBopt.Dynamics,'particle')
            line([min(R) max(R)],[H0 H0],'Color','k','LineStyle','--');
        else
            line([min(R) max(R)],[nu*H0 nu*H0],'Color','k','LineStyle','--');
        end
        line([r0 r0],[min(W) max(W)],'Color','k','LineStyle','--');
        xlabel('r');
        legend('W[p_{\phi}(0)]','W[j(0)]',1);
        title('EOB potential')
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolve

fprintf(' ===> Compute dynamics\n');
tstart = tic;
[t,y] = ode113(@(t,y) EOBopt.ComputeRhs(t,y, nu,EOBopt), tspan, y0, ODEopt);
fprintf(' Done %4.3f sec\n',toc(tstart));

if verbose
            
    if verbose==2
        
        fhIdx  = fhIdx+1; fhName{fhIdx} = 'Dyn_rt'; fh(fhIdx) = figure;
        plot(t,y(:,2))
        xlabel('t/M'); ylabel('r/M');
        title('EOB dynamics r(t)')
        
        fhIdx  = fhIdx+1; fhName{fhIdx} = 'Dyn_rphi'; fh(fhIdx) = figure;
        polar(y(:,1),y(:,2))
        title('EOB dynamics r(\phi)')
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack and Extend variables (if necessary)

nt     = length(t);

phi    = y(:,1);
r      = y(:,2);
pph    = y(:,3);
prstar = y(:,4);

text = EOBopt.textend;
tdyn = t(end);

text = [(tdyn+dt):dt:(tdyn+text)]';
t    = [t; text];
nt   = length(t);

phi     = phi(end)   *ones(nt,1);
r       = r(end)     *ones(nt,1);
pph     = pph(end)   *ones(nt,1);
prstar  = prstar(end)*ones(nt,1);

Idx = (t<=tdyn);
phi(Idx)    = y(:,1);
r(Idx)      = y(:,2);
pph(Idx)    = y(:,3);
prstar(Idx) = y(:,4);

y = [phi r pph prstar]';



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy, angular momentum and frequency along the EOB dynamics

fprintf(' ===> Compute energy, angular momentum and orbital frequency\n');
tstart = tic;
Metric            = EOBMetric( nu, r, EOBopt );
[Omega E pr Heff] = EOBOmgE( nu, r,pph,prstar, EOBopt, Metric );
fprintf(' Done %4.3f sec\n',toc(tstart));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of orbital frequency and ddotr

fprintf(' ===> Compute derivatives of orbital frequency ...\n');
tstart = tic;
dotOmega = FDdrvt(Omega,t,4).';
A_Omg    = dotOmega./(96*nu/5*Omega.^(11/3));
Q_Omg    = Omega.^2./dotOmega;
fprintf(' Done %4.3f sec\n',toc(tstart));

fprintf(' ===> Compute ddot r ...\n');
switch EOBopt.ddotrMethodWav
    case 'boot2'
        Ham   = EOBopt.ComputeEOBHam(nu,r,pph,prstar,EOBopt,Metric);
        dydt  = EOBRhsCons(t,y, nu, EOBopt, Metric, Ham);
        ddotr = EOBddotrBoot2(nu, y,dydt, EOBopt, Metric, Ham);
    case 'noflux'
        Ham   = EOBopt.ComputeEOBHam(nu,r,pph,prstar,EOBopt,Metric);
        dydt  = EOBRhsCons(t,y, nu, EOBopt, Metric, Ham);
        ddotr = EOBddotrNoFlux(nu, y,dydt, EOBopt, Metric, Ham);
    case 'fd'
        if strcmp(EOBopt.Dynamics,'particle')
            dotr = Metric.A.*prstar./Heff;
        elseif strcmp(EOBopt.PNorder,'1pn') || strcmp(EOBopt.PNorder,'2pn')
            dotr = sqrt(Metric.A ./ Metric.B).*prstar./(E.*Heff);
        else
            dotr = sqrt(Metric.A./Metric.B).*(prstar + ...
                4*nu*(4-3*nu)*Metric.A.*(1./r).^2.*prstar.^3)./(E.*Heff);
        end
        ddotr = FDdrvt(dotr,t,4);
    otherwise
        error('unknown option %s for ddotrMethodWav.',EOBopt.ddotrMethodWav);
end
ddotr = reshape(ddotr,size(r));
fprintf(' Done %4.3f sec\n',toc(tstart));

% Pack dynamics
dyn.t      = t;
dyn.r      = r;
dyn.phi    = phi;
dyn.pph    = pph;
dyn.E      = E;
dyn.e      = (E-1)/nu;
dyn.prstar = prstar;
dyn.pr     = pr;
dyn.ddotr  = ddotr;
dyn.Omg    = Omega;
dyn.Heff   = Heff;
dyn.dotOmg = dotOmega;
dyn.A_Omg  = A_Omg;
dyn.Q_Omg  = Q_Omg;
      

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSO and LR crossing (if required)

if verbose
    
    fprintf(' ===> Compute LSO and LR crossing\n');
    tstart = tic;
    
    % Locate the crossing of adiabatic LSO position and write into LSO struct
    % NOTE: upper case vars are the dynamical ones !
    idxLSO     = find(r <= LSO.r,1,'last');
    LSO.idxLSO = idxLSO;
    LSO.T      = t(idxLSO);
    LSO.OMEGA  = Omega(idxLSO);
    LSO.R      = r(idxLSO);
    LSO.U      = 1/r(idxLSO);
        
    % Locate the maximum of orbital frequency (EOB "light-ring crossing")
    [MaxOmega, idxLR] = max(Omega);
    
    jIdx = idxLR-5:idxLR+5;    
    [tmax,MaxOmega] = FindMax( t(jIdx), Omega(jIdx), 2 );
    LR.idxLR   = idxLR;
    LR.OMEGA   = MaxOmega;
    LR.omega22 = 2 * MaxOmega;
    LR.T   = tmax;
    LR.R   = LagInt1d( 4, t, r, tmax );
    LR.U   = 1/LR.R;
    LR.E   = LagInt1d( 4, t, dyn.E  , tmax );
    LR.pph = LagInt1d( 4, t, dyn.pph, tmax );
    
    % Pack data
    dyn.LSO = LSO;
    dyn.LR = LR;
    
    fprintf(' Done %4.3f sec\n',toc(tstart));
    
    fprintf(' t            [LSO Dyn] = %+.12e\n',LSO.T);
    fprintf(' Omega        [LSO Dyn] = %+.12e\n',LSO.OMEGA);
    fprintf(' r            [LSO Dyn] = %+.12e\n',LSO.R);
    fprintf(' u            [LSO Dyn] = %+.12e\n',LSO.U);
    
    fprintf(' t            [LR]      = %+.12e\n',LR.T);
    fprintf(' Omega        [LR]      = %+.12e\n',LR.OMEGA);
    fprintf(' omega22      [LR]      = %+.12e\n',LR.omega22);
    fprintf(' r            [LR]      = %+.12e\n',LR.R);
    fprintf(' u            [LR]      = %+.12e\n',LR.U);
    fprintf(' E            [LR]      = %+.12e\n',LR.E);
    fprintf(' pph          [LR]      = %+.12e\n',LR.pph);
    
    if verbose==2
        
        fhIdx  = fhIdx+1; fhName{fhIdx} = 'EOBdyn_epph'; fh(fhIdx) = figure;
        plot(pph,(E-1)/nu,'k')
        title('Binding energy vs angular momentum');
        xlabel('p_\phi');
        ylabel('e = (E-1)/\nu');
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save (if required) ... End here ?

if saveme
    
    % *.mat file(s)
    fprintf(' ===> Save *.mat files ...\n');
    tstart = tic;        
    
    filename  = strcat(basedir,'Settings');
    save(filename,'EOBopt');
    [err mess] = system(sprintf('du -h %s.mat ',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);    
    
    filename  = strcat(basedir,'Dynamics');
    save(filename,'dyn');
    [err mess] = system(sprintf('du -h %s.mat ',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);    
    
    fprintf(' Done %4.3f sec\n',toc(tstart));
    
    if verbose==2
        % *.fig(s)        
        fprintf(' ===> Save *.fig files ...\n');
        for k=1:fhIdx
            saveas(fh(k),strcat(basedir,'/Fig/',fhName{fhIdx}),'fig');
        end        
        fprintf(' Done %4.3f sec\n',toc(tstart));
        %close all;  
        fhIdx = 0;      
    end
    
end

if modus==0
    varargout = SetVOutput( nargout, EOBopt,dyn );
    fprintf('End %4.3f sec\n',toc(tottime));
    return
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insplunge waveform

fprintf(' ===> Insplunge waveform ...\n');
tstart = tic;
wav = EOBhlm(nu, t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, Metric);
fprintf(' Done %4.3f sec\n',toc(tstart));

wav.t = t;

if verbose
  
    fprintf(' # of orbits     = %g\n',phi(end)/(2*pi));

    if verbose==2
        
        fhIdx  = fhIdx+1; fhName{fhIdx} = 'Wav_inspl'; 
        fh(fhIdx) = figure('Name','Insplunge waveform');
        plot(t,abs(wav.psilm(LM2K(2,2),:)),'b-',...
            t,wav.omglm(LM2K(2,2),:),'r--',...
            t,Omega,'g-.')        
        title('Insplunge waveform');
        xlabel('t/M');
        ylim([0 0.56]);
        legend('|\Psi_{22}|/\nu','M\omega_{22}','M\Omega');        
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save (if required) ... End here ?

if saveme
    
    fprintf(' ===> Save *.mat files ...\n');
    tstart = tic;        
    
    filename  = strcat(basedir,'WavInspl');
    save(filename,'wav');
    [err,mess] = system(sprintf('du -h %s.mat ',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);
    
    fprintf(' Done %4.3f sec\n',toc(tstart));
    
    if verbose==2
        % *.fig(s)        
        fprintf(' ===> Save *.fig files ...\n');
        for k=1:fhIdx
            saveas(fh(k),strcat(basedir,'/Fig/',fhName{fhIdx}),'fig');
        end        
        fprintf(' Done %4.3f sec\n',toc(tstart));
        %close all;    
        fhIdx = 0; 
    end
    
end

if modus==1
    varargout = SetVOutput( nargout, EOBopt,dyn,wav );    
    fprintf('End %4.3f sec\n',toc(tottime));
    return
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insplunge+merger waveform

% Determine a's and b's (if necessary)
if strcmp(EOBopt.DetermineNQCab,'yes')    
    fprintf(' ===> Determine QNM from NR data ...\n');
    tstart = tic;    
    n = struct2cell(EOBNQCn(r,prstar, Omega,ddotr));    
    [a,b] = EOBNQCabFind( t,Omega, n, ...
        wav.Alm, wav.omglm, wav.domglm, wav.d2omglm,  NRPTS );    
    fprintf(' Done %4.3f sec\n',toc(tstart));
end

% NQC corrections
fprintf(' ===> Insplunge+NQC correction waveform ...\n');
tstart = tic;
psilmnqc = EOBhlmNQC( nu, r,prstar, Omega,ddotr, a,b, EOBopt );
fprintf(' Done %4.3f sec\n',toc(tstart));

% Update  waveform
fprintf(' ===> Update waveform ...\n');
if strcmp(EOBopt.Dynamics,'eob')
    cnu = 12/5*2.0^(1/3)*nu;
else
    cnu = 1; % FIXME ?
end
nell = sqrt((L+2).*(L+1).*L.*(L-1))' * ones(1,nt);

t         = wav.t;
psilm     = wav.psilm .* psilmnqc;
hlm       = nell.*psilm;

wav          = EOBExtractAPO(psilm,t);
wav.t        = t;
wav.psilm    = psilm;
wav.hlm      = hlm;
wav.psilmnqc = psilmnqc;

wav.A_omg  = wav.domglm./( cnu*wav.omglm.^(11/3) ) ;
wav.Q_omg  = wav.omglm.^2./wav.domglm;

wav.NQCa = a; 
wav.NQCb = b; 

fprintf(' Done %4.3f sec\n',toc(tstart));


if verbose
    
    if verbose==2
        
        fhIdx     = fhIdx+1; fhName{fhIdx} = 'Wav_inspl'; 
        fh(fhIdx) = figure('Name','Insplumerg waveform');
        plot(t,abs(wav.psilm(LM2K(2,2),:)),'b-',...
            t,wav.omglm(LM2K(2,2),:),'r--',...
            t,Omega,'g-.')
        ylim([-0.32 0.55]);        
        title('Insplunge+NQC waveform');
        xlabel('t/M');
        legend('|\psi_{22}|/\nu','\omega_{22}','\Omega')
        
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save (if required) ... End here ?

if saveme
    
    fprintf(' ===> Save *.mat files ...\n');
    tstart = tic;    
        
    filename  = strcat(basedir,'WavInsplMerg');
    save(filename,'wav');
    [err mess] = system(sprintf('du -h %s.mat',filename));
    if err, error(mess); end;
    fprintf(' %s',mess);
    
    if strcmp(EOBopt.DetermineNQCab,'yes')
        filename  = strcat(basedir,'NQCcorr');
        save(filename,'a','b');
        [err mess] = system(sprintf('du -h %s.mat',filename));
        if err, error(mess); end;
        fprintf(' %s',mess);        
    end
    
    fprintf(' Done %4.3f sec\n',toc(tstart));
     
    if verbose==2
        % *.fig(s)        
        fprintf(' ===> Save *.fig files ...\n');
        for k=1:fhIdx
            saveas(fh(k),strcat(basedir,'/Fig/',fhName{fhIdx}),'fig');
        end        
        fprintf(' Done %4.3f sec\n',toc(tstart));
        %close all;   
        fhIdx = 0; 
    end   

end

if modus==2
    varargout = SetVOutput( nargout, EOBopt,dyn,wav );       
    fprintf('End %4.3f sec\n',toc(tottime));
    return
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insplunge + merger + ringdown waveform

ferrs = @(s)(Compute_Mbh_abh_errs(nqnm,lmaxqnm, EOBopt,s(1)*LR.E,s(2)*LR.pph*nu/(s(1)*LR.E)^2,Omega,t, dt, wav,nu,H0,pph0));

%ferrs = @(s)(Compute_Mbh_abh_errs(nqnm,lmaxqnm, EOBopt,(1+atan(s(1))/(pi/2))*LR.E/2,(1+atan(s(2))/(pi/2))*LR.pph*nu/2/((1+atan(s(1)))*LR.E/2)^2,Omega,t, dt, wav,nu,H0,pph0));

%s0 = [1 - 0.02/LR.E, (1 - 0.06/(LR.pph*nu))/(1 - 0.02/LR.E)^2];

[Mbh_fit, abh_fit] = EOBFinalBH(nu);

%Mbh_fit = 0.866; %0.878; %0.896; %0.928; % Final value for q = 1, a2 = 300 (a2 = 200, 100, 20 commented out), as an initial guess for higher a2
%abh_fit = 0.351; %0.410; %0.495; %0.621; % Ibid.

%Mbh_fit = 0.892; %0.901; %0.916; %0.942; % Final value for q = 2, a2 = 300, (a2 = 200, 100, 20 commented out), as an initial guess for higher a2
%abh_fit = 0.311; %0.363; %0.440; %0.548; % Ibid.

%Mbh_fit = 0.883; %0.913; % Final value for q = 1, iQNMfac = 0.5, with 0.6 commented out, used as an initial guess for smaller values of iQNMFac
%abh_fit = 0.486; %0.585; % Ibid.

%Mbh_fit = 0.86; % Initial guess values for q = 1, iQNMfac = 0.45 
%abh_fit = 0.3; % Ibid.

%Mbh_fit = 0.915; % Final value for q = 2, iQNMfac = 0.6, used as an initial guess for smaller values of iQNMFac
%abh_fit = 0.486; % Ibid.

%Mbh_fit = 0.87; % Initial guess values for q = 2, iQNMfac = 0.5--didn't converge
%abh_fit = 0.4; % Ibid.

s0 = [Mbh_fit/LR.E, abh_fit/(LR.pph*nu/LR.E^2)]; % Using the fit to give an initial guess

%s0 = [0.98, ]; % Using the 

%s0t = (pi/2).*tan(2.*s0 - 1); 

%optionss = optimset('Display','iter','MaxIter',200,'PlotFcns',{@optimplotx,@optimplotfval});
optionss = optimset('Display','iter','MaxIter',200);

[sout, ferrsout, exitflagsout, outputsout] = fminsearch(ferrs, s0,optionss);

%soutp = (1 + atan(sout)./(pi/2))./2;

Mbhabhout = [sout(1)*LR.E,sout(2)*LR.pph*nu/(sout(1)*LR.E)^2];

fprintf('Mbh (LR) = %+.12e\n',LR.E);
fprintf('abh (LR) = %+.12e\n',LR.pph*nu/LR.E^2);
fprintf('Mbh (consistent) = %+.12e\n',sout(1)*LR.E);
fprintf('abh (consistent) = %+.12e\n',sout(2)*LR.pph*nu/(sout(1)*LR.E)^2);

varargout = SetVOutput( nargout, EOBopt,dyn,wav,Mbhabhout,ferrsout,exitflagsout,outputsout);
fprintf('End %4.3f sec\n',toc(tottime));
