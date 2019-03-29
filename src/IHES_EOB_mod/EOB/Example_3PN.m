% Example of EOB run: specify to order 3PN 


% 1. Choose input parameters (these are needed!)

% mass ratio
q = 1;

% Sym mass ratio
nu  = q/((1+q)^2);
nu2 = nu^2;

% EOB free pars
% a_i's for radiation reaction (from fit)
a1 = 0; 
a2 = 0; 
a3 = 0; 
a4 =  0;
% 4PN and 5PN effective corrections
a5 =  23.5;
a6 = -122+147*(1-4*nu);
% pack them in ainput (see explanation in example.m)
ainput{1} = a1;
ainput{2} = a2;
ainput{3} = a3;
ainput{4} = a4;
ainput{5} = a5; 
ainput{6} = a6; 

r0   = 16;   % Initial radius
tmax = 5000; % Final time
dt   = 0.7; % Width of the matching comb

nQNM = [5 0]; % no QNM [ positive negative ]
lmaxQNM = 2; % lmax for QNM matching


% 2. Choose optional input

options = EOBSet('PNorder','3pn');

% use default options
%options = []; % Set options in this struct (see EOBSet.m)

% Reminder of all default options
%{
options =EOBSet('Dynamics','eob',...
    'RadReac','din',...
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
    'NQCFileNR','./DataNR/NRpts4NQC.dat',...
    'QNMDataDir','./dataQNMs/',....
    'textend',100,...
    'ODESolverRelTol',1e-9,...
    'ODESolverAbsTol',1e-12,...
    'ODESolverRefine',1);
%}
modus   = 2; % (1=> inspl) (2=> inspl+mrg) (3=> inspl+mrg+rng)
verbose = 2; % 0 = minimum, 1 = text ,2 = text+figs
saveme  = 1; % save data ?
outputdir = './tmp';


% 3. Get down to work 

out = EOBRun(nu, ainput, r0,tmax,dt, nQNM,lmaxQNM, ...
    ... the following pars are optional 
    options,... 
    modus,...
    verbose, ...
    saveme, ...
    'tmp')
 

