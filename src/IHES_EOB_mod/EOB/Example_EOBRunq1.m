% Main example of EOB run


% 1. Choose input parameters (these are needed!)

% mass ratio
q = 1;

% Sym mass ratio
nu  = q/((1+q)^2);
nu2 = nu^2;

% EOB free pars
% Note: these values are used only with the option 
%         EOBopt.DetermineNQCab = 'yes'
% as best guesses for the dynamics, then overwritten by determination
% procedure. Otherwise the values used for both dynamics and waves are
% those taken from the fit, ie as output of EOBNQCabFit.m
%  
% a_i's for radiation reaction (from fit)
a1 =  2.160054e+00*nu2 - 1.093667e+00*nu + 7.926107e-02;
a2 = -1.080736e+01*nu2 + 7.141956e+00*nu + 7.034821e-01;
a3 = -2.766573e+00*nu2 - 1.768795e-01*nu + 1.012170e-01;
a4 =  0;
% 4PN and 5PN effective corrections
a5 =  23.5;
a6 = (-110.5 + 129*(1-4*nu)).*sqrt(1 - 1.5e-5/((nu-0.26).^2));
% pack them in ainput
% Remark: each a_i i = 1...4 can contain all the multipoles so, here below,
% the order matters ! A quick way to get the correct index without thinking
% is to use the Multipidx.m function. As an example in the following we
% pass the 22 mode only:
LM2K = Multipidx(2);
ainput{1} = zeros(1,length(LM2K)); 
ainput{2} = zeros(1,length(LM2K)); 
ainput{3} = zeros(1,length(LM2K)); 
ainput{4} = zeros(1,length(LM2K)); 
ainput{1}(LM2K(2,2)) = a1;
ainput{2}(LM2K(2,2)) = a2;
ainput{3}(LM2K(2,2)) = a3;
ainput{4}(LM2K(2,2)) = a4;

ainput{5} = a5; 
ainput{6} = a6; 

r0   = 16;   % Initial radius
tmax = 10000; % Final time

dt   = 0.7;  % Width of the matching comb
nQNM = [5 0]; % no QNM [ positive negative ]
lmaxQNM = 4;  % lmax for QNM matching


% 2. Choose optional input

options = []; % Set options in this struct (see EOBSet.m)

% Reminder of all default options

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
    'NQCFileNR','/Users/nkjm/PN/EOB/EOB_1202_5PNlogNospin/DataNR/NRpts4NQC.dat',...
    'QNMDataDir','/Users/nkjm/PN/EOB/EOB_1202_5PNlogNospin/DataQNM/',....
    'textend',100,...
    'ODEEventFunRmin',1.0,...
    'ODESolverRelTol',1e-9,...
    'ODESolverAbsTol',1e-12,...
    'ODESolverRefine',1);

modus   = 3; % (1=> inspl) (2=> inspl+mrg) (3=> inspl+mrg+rng)
verbose = 2; % 0 = minimum, 1 = text ,2 = text+figs
saveme  = 1; % save data ?
outputdir = './tmp';

% 3. Get down to work 

out = EOBRun(nu, ainput, r0,tmax,dt, nQNM,lmaxQNM, options, modus, verbose,  saveme, 'tmp')
 

