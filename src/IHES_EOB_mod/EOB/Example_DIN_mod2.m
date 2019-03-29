% Main example of EOB run

% Read in masses for injections; this expects a two-column file of masses
% separated by commas. Note that the current setting of r0 = 34 would need
% to be changed if the minimum mass was below (for flow = 10 Hz).

warning off;

datfile = fopen('/home/abhirup/Documents/Work/testGR_IR/scripts/mass_list_popsynth.txt','r');
injfile = fopen('/home/abhirup/Documents/Work/testGR_IR/scripts/inj_thresh8_matlab.txt', 'r');

m1m2 = fscanf(datfile, '%f %f', [2 Inf]);
inj_list = fscanf(injfile, '%f', [2 Inf]);

fclose(datfile)
fclose(injfile)

% Iterate over masses

for n = inj_list(126:150)%1:length(m1m2)

    fprintf('-----------------------\n');
    
    fprintf('Injection: %i\n',n);
    
    % 1. Choose input parameters (these are needed!)
    
    % mass ratio
    
    q = max(m1m2(2,n)/m1m2(1,n),m1m2(1,n)/m1m2(2,n));

    mtot = m1m2(1,n) + m1m2(2,n);
    
    fprintf('Mass ratio: %+.12e\n',q);
    
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
    LM2K = Multipidx(7); %Multipidx(2);
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
    
    r0   = 34; %32; %80; %32; % 16;  % Initial radius
    %tmax = 240000; %110000; %500000; %2100000; %100000; % 10000 % Final time
    tmax = (5/256)*r0^4/nu+10000; % "Newtonian" coalescence time, plus an added margin for safety
    
    dt   = 0.7; %0.005; %0.7;  % Width of the matching comb
    nQNM = [5 0]; % no QNM [ positive negative ]
    lmaxQNM = 7;  % lmax for QNM matching
    
    
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
        'a2',8000./(mtot*mtot),...
        'iQNMfac',1.0,...
        'NQCFileNR','/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/DataNR/NRpts4NQC.dat',...
        'QNMDataDir','/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/DataQNM',....
        'textend',300,...
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
    savedat = 1; % save data (in text format)?
    
    % 3. Get down to work
    
    out = EOBRun_DIN_mod2(nu, ainput, r0,tmax,dt, nQNM,lmaxQNM, options, modus, verbose,  saveme, 'tmp_DIN_mod_q2s_fac2')
   
    Mbh = out.Mbhabhout(1);
    abh = out.Mbhabhout(2);
    Mbhrad= out.Mbhabhradout(1);
    abhrad = out.Mbhabhradout(2);

    % Output t, Re h22, Im h22 (and higher modes)
    
    if savedat
        
        modfac = nu*out.EOBopt.fac^0.5; % Factor to multiply the waveforms to take the modification into account (currently just includes fac); also includes the factor of nu
        
        %outfile = fopen('/home/nathan/testGR_IR/src/IHES_EOB_mod/waveforms/q2_a2_400_r0_32.dat','w');
        %fprintf(outfile, '%-10s %-13s %-13s\n', 't','Re h22', 'Im h22');
        %fprintf(outfile, '%.4f %.6e %.6e\n',vertcat(out.wav.t.',modfac.*real(out.wav.hlm(LM2K(2,2),:)),modfac.*imag(out.wav.hlm(LM2K(2,2),:))));
        %fclose(outfile);
        
        %outfile = fopen(strcat('/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/waveforms/popsynth_modGR/injection_',int2str(n),'_q_',num2str(q,3),'_a2_',num2str(out.EOBopt.a2,3),'_.dat'),'w');
        outfile = fopen(strcat('/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/waveforms/popsynth_modGR_mtotalmod/injection_',int2str(n),'_q_',num2str(q,3),'_.dat'),'w');
        %fprintf(outfile, '%-10s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s\n', 't','Re h22', 'Im h22','Re h21', 'Im h21','Re h33', 'Im h33','Re h32', 'Im h32','Re h31', 'Im h31','Re h44', 'Im h44','Re h43', 'Im h43','Re h42', 'Im h42','Re h41', 'Im h41','Re h55', 'Im h55','Re h54', 'Im h54','Re h53', 'Im h53','Re h52', 'Im h52','Re h51', 'Im h51');
        fprintf(outfile, '%-10s %-13s %-13s\n', 't','Re h22', 'Im h22');
        %fprintf(outfile, '%.4f %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n',vertcat(out.wav.t.',modfac.*real(out.wav.hlm(LM2K(2,2),:)),modfac.*imag(out.wav.hlm(LM2K(2,2),:)),modfac.*real(out.wav.hlm(LM2K(2,1),:)),modfac.*imag(out.wav.hlm(LM2K(2,1),:)),modfac.*real(out.wav.hlm(LM2K(3,3),:)),modfac.*imag(out.wav.hlm(LM2K(3,3),:)),modfac.*real(out.wav.hlm(LM2K(3,2),:)),modfac.*imag(out.wav.hlm(LM2K(3,2),:)),modfac.*real(out.wav.hlm(LM2K(3,1),:)),modfac.*imag(out.wav.hlm(LM2K(3,1),:)),modfac.*real(out.wav.hlm(LM2K(4,4),:)),modfac.*imag(out.wav.hlm(LM2K(4,4),:)),modfac.*real(out.wav.hlm(LM2K(4,3),:)),modfac.*imag(out.wav.hlm(LM2K(4,3),:)),modfac.*real(out.wav.hlm(LM2K(4,2),:)),modfac.*imag(out.wav.hlm(LM2K(4,2),:)),modfac.*real(out.wav.hlm(LM2K(4,1),:)),modfac.*imag(out.wav.hlm(LM2K(4,1),:)),modfac.*real(out.wav.hlm(LM2K(5,5),:)),modfac.*imag(out.wav.hlm(LM2K(5,5),:)),modfac.*real(out.wav.hlm(LM2K(5,4),:)),modfac.*imag(out.wav.hlm(LM2K(5,4),:)),modfac.*real(out.wav.hlm(LM2K(5,3),:)),modfac.*imag(out.wav.hlm(LM2K(5,3),:)),modfac.*real(out.wav.hlm(LM2K(5,2),:)),modfac.*imag(out.wav.hlm(LM2K(5,2),:)),modfac.*real(out.wav.hlm(LM2K(5,1),:)),modfac.*imag(out.wav.hlm(LM2K(5,1),:))));
        fprintf(outfile, '%.4f %.6e %.6e\n',vertcat(out.wav.t.',modfac.*real(out.wav.hlm(LM2K(2,2),:)),modfac.*imag(out.wav.hlm(LM2K(2,2),:))));
        fclose(outfile);
        
	%outfile2 = fopen(strcat('/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/waveforms/popsynth_modGR/injection_',int2str(n),'_q_',num2str(q,3),'_a2_',num2str(out.EOBopt.a2,3),'_Mbhabh.dat'),'w');
	outfile2 = fopen(strcat('/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/waveforms/popsynth_modGR_mtotalmod/injection_',int2str(n),'_q_',num2str(q,3),'_Mbhabh_.dat'),'w');
	fprintf(outfile2, ' Mbh      [LSO Dyn] = %+.12e\n',Mbh);
	fprintf(outfile2, ' abh      [LSO Dyn] = %+.12e\n',abh);

	fprintf(outfile2,' Mbh (rad) = %+.12e\n', Mbhrad);
    	fprintf(outfile2,' abh (rad) = %+.12e\n', abhrad);

    	fprintf(outfile2,'Mbh fractional error = %+.2e\n',Mbh/Mbhrad-1);
    	fprintf(outfile2,'abh fractional error = %+.2e\n',abh/abhrad-1);

	fclose(outfile2);

        %outfile2 = fopen('/Users/nkjm/testing_GR/testGR_IR/src/IHES_EOB_mod/waveforms/q1_fac0p5_E_vs_j_test.dat','w');
        %fprintf(outfile2, '%-10s %-13s\n', 'j', 'E');
        %fprintf(outfile2, '%.6e %.6e\n',vertcat(out.dyn.pph.',(out.dyn.E.' - 1)./nu));
        %fclose(outfile2);
        
    end

end

