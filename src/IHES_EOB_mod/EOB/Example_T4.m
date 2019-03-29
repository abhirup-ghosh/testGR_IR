% Example of T4 run


% 1. Choose input parameters (these are needed!)
q = 1;

r0 = 16; % Initial radius
Phi0 = 0;

tmax = 5000; % Final time
dt = 0.5; % spacing time vector


% 2. Choose optional input
options = []; % Set ODE options 
pnorder = '3pn';
addelta = 'no';
computePAO = 'yes';
verbose = 0; % 0 = minimum, 1 = text ,2 = text+figs
saveme  = 0; % save data ?
outputdir = './tmp';


% 3. Get down to work 
[dyn,wav] = T4Run(q,r0,Phi0,dt,tmax,...
    ... the following pars are optional 
    options,...   
    pnorder,...
    addelta,...
    computePAO,...
    verbose, ...
    saveme, ...
    'tmp')

figure
plot(wav.t,wav.h22/0.25)