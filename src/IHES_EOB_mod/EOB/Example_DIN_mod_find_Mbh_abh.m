% Main example of EOB run


% 1. Choose input parameters (these are needed!)

% mass ratio
q = 1.775;

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

r0   = 16; %80; %32; % 16;  % Initial radius
%tmax = 100000; %2100000; %100000; % 10000 % Final time
tmax = (5/256)*r0^4*(1+q)^2/q+10000; % "Newtonian" coalescence time, plus an added margin for safety

dt   = 0.7;  % Width of the matching comb
nQNM = [5 0]; % no QNM [ positive negative ]
lmaxQNM = 7; %4; % lmax for QNM matching


% 2. Choose optional input

options = []; % Set options in this struct (see EOBSet.m)

% Reminder of all default options

options =EOBSet_mod('Dynamics','eob',...
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
    'a2',20.0,...
    'iQNMfac',1.0,...
    'NQCFileNR','/home/nathan/testGR_IR/src/IHES_EOB_mod/DataNR/NRpts4NQC.dat',...
    'QNMDataDir','/home/nathan/testGR_IR/src/IHES_EOB_mod/DataQNM/',....
    'textend',100,...
    'ODEEventFunRmin',1.75,...
    'ODESolverRelTol',1e-9,...
    'ODESolverAbsTol',1e-12,...
    'ODESolverRefine',1);

    %'Cp',0.0,...
    %'pow',1.5,...
    %'Ce',0.0,...
    %'g',0.25,...

modus   = 3; % (1=> inspl) (2=> inspl+mrg) (3=> inspl+mrg+rng)
verbose = 1; % 0 = minimum, 1 = text ,2 = text+figs
saveme  = 0; % save data (in Matlab format)?
savedat = 0; % save data (in text format)?

% 3. Get down to work 

out = EOBRun_DIN_mod_find_Mbh_abh(nu, ainput, r0,tmax,dt, nQNM,lmaxQNM, options, modus, verbose,  saveme, 'tmp_DIN_mod_q2s_fac2')
 
% Output t, Re h22, Im h22

if savedat
    
    modfac = nu*out.EOBopt.fac^0.5; % Factor to multiply the waveforms to take the modification into account (currently just includes fac); also includes the factor of nu
    
    outfile = fopen('/Users/nkjm/testing_GR/testGR_IR/src/IHES_EOB_mod/waveforms/q1_fac1_GR_final_M_and_a.dat','w');
    fprintf(outfile, '%-10s %-13s %-13s\n', 't','Re h22', 'Im h22');
    fprintf(outfile, '%.4f %.6e %.6e\n',vertcat(out.wav.t.',modfac.*real(out.wav.hlm(LM2K(2,2),:)),modfac.*imag(out.wav.hlm(LM2K(2,2),:))));
    fclose(outfile);
    
    outfile2 = fopen('/Users/nkjm/testing_GR/testGR_IR/src/IHES_EOB_mod/waveforms/q1_fac0p5_E_vs_j_test.dat','w');
    fprintf(outfile2, '%-10s %-13s\n', 'j', 'E');
    fprintf(outfile2, '%.6e %.6e\n',vertcat(out.dyn.pph.',(out.dyn.E.' - 1)./nu));
    fclose(outfile2);
    
end

