% Example of EOB run: run the Matlab profiler


% 1. Choose input parameters (these are needed!)

% Sym mass ratio
nu  = 1/4;
nu2 = nu.^2;

% EOB free pars
a1 =  2.160054e+00*nu2 - 1.093667e+00*nu + 7.926107e-02;
a2 = -1.080736e+01*nu2 + 7.141956e+00*nu + 7.034821e-01;
a3 = -2.766573e+00*nu2 - 1.768795e-01*nu + 1.012170e-01;
a4 = 0;
a5 = +23.5;
a6 = (-110.5 + 129*(1-4*nu)).*sqrt(1 - 1.5e-5/((nu-0.26).^2));
ainput{1} = a1;
ainput{2} = a2;
ainput{3} = a3;
ainput{4} = a4;
ainput{5} = a5; 
ainput{6} = a6; 

r0 = 16; % Initial radius
tmax = 5000; % Final time
dt = 0.5; % spacing time vector

nQNM = [5 0]; % no QNM [ positive negative ]
lmaxQNM = 2; % lmax for QNM matching


% 2. Choose optional input

options = []; % Set options in this struct (see EOBSet.m)
modus   = 3; 
verbose = 0; % 0 = minimum, 1 = text ,2 = text+figs
saveme  = 0; % save data ?
outputdir = './tmp';


% 3. Get down to work 
profile on

[options,dyn,wav] = EOBRun(nu, ainput, r0,tmax,dt, nQNM,lmaxQNM, ...
    ... the following pars are optional 
    options,... 
    modus,...
    verbose, ...
    saveme, ...
    'tmp')

profile off


% 4. Check what to improve
profile viewer
%profsave(profile('info'),strcat(datestr(now, 'yyyymmdd'),'_myprofile_results'))
