function [hatF Flm FNewtlm flm] = ...
    EOBFluxDIN_mod(x,Omega,E,Heff,jhat,nu,r,prstar,ddotr, EOBopt)

%EOBFluxDIN Computes the resummed Newton-normalized GW energy flux.
%
%   [hatF Flm FNewtlm flm] = ...
%   EOBFluxDIN(x,Omega,E,Heff,jhat,nu,r,prstar,ddotr, EOBopt) 
%
%   Reference(s)
%    Damour, Iyer & Nagar, PRD 79, 064004 (2008) 
%


L                 = EOBopt.L;
M                 = EOBopt.M;
LM2K              = EOBopt.LM2K;

rholmarg          = EOBopt.rholmarg;
f22_resum         = EOBopt.resumf22;
NQC               = EOBopt.NQC;

particle_dynamics = strcmp(EOBopt.Dynamics,'particle');


% Compute amplitude of resummed waveform

% Argument for the modulus of hhatlm
% Note this has to be used only in the hhatlm modulus calculations, 
% everywhere else use x !
switch rholmarg    
    case 'v_omg'
        xarg = Omega.^(2/3);
    case 'v_phi'
        xarg = x;
    case '1_r'
        xarg = 1./r;
    otherwise
        error('unknown option %s for rholmarg',rholmarg)
end

if ~particle_dynamics
        
    % Compute f_lm = (rho_lm)^ell    
    flm = EOBflm(xarg,nu,EOBopt);
        
    % Further resummation of rho_22 or f_22 (if needed)
    if strcmp(f22_resum,'pade23')        
        flm(LM2K(2,2),:) = EOBrho22Pade23(xarg,nu).^2;   
    elseif strcmp(f22_resum,'oldpade')   
        flm(LM2K(2,2),:) = EOBf22Pade(xarg,nu);    
    end    
    
    % Compute NQC correction to the modulus of the (l,m) waveform        
    [n1,n2,n3]= EOBNQCn(r,prstar, Omega,ddotr);
    %warning(' check routine EOBNQCn and use it here ');
    %{ 
    n1 = (prstar./(r.*Omega)).^2;
    n2 = ddotr./(r.*Omega.^2);
    n3 = n1.*prstar.^2;
    %}
    
    hnqc_22 = (1+NQC.a{1}(LM2K(2,2))*n1 + NQC.a{2}(LM2K(2,2))*n2 + NQC.a{3}(LM2K(2,2))*n3);
    hnqc_21 = (1+NQC.a{1}(LM2K(2,1))*n1 + NQC.a{2}(LM2K(2,1))*n2 + NQC.a{3}(LM2K(2,1))*n3);
    hnqc_33 = (1+NQC.a{1}(LM2K(3,3))*n1 + NQC.a{2}(LM2K(3,3))*n2 + NQC.a{3}(LM2K(3,3))*n3);

    % FIXME: to include NQC corrections for all multipoles one should just
    % use something like
    %
    %hnqc = (1+NQC.a1 * n1 + NQC.a2 * n2 + NQC.a3 * n3);
    %
    % where hnqc is the matrix of size = (no multipoles )x( length(Omega) )
    % Code below should be fixed accordingly.

else
        
    % Compute flm = (rho_lm)^ell
    flm  = EOBflm0(xarg,EOBopt);    
    
    hnqc_22 = 1;
    hnqc_21 = 1;
    hnqc_33 = 1;        
    
end


% Compute The PN correction to the waveform amplitude
% Note that jhat must be put to 1 when 1 has factored
% the full angular momentum J Omega in the Newtonian part 
% Beware this is the absolute value of the \hat{h}_lm
%MTlm  = EOBModTail( L, M, E.*Omega );
% MEX-file version ( speed - O(100) )
MTlm  = EOBModTailm( L, M, E.*Omega );


% Calculate prefactor (extend vectors to matrices) 
kmax   = length(L);
nx     = length(xarg);
lmeven = (~mod(L+M,2))';
lmodd  = ( mod(L+M,2))';
Heff   = reshape(Heff,1,nx); % do not resize back, if not needed below
jhat   = reshape(jhat,1,nx);              

prefact = lmeven*Heff + lmodd*jhat; 


% Compute modulus of hhat_lm 
Modhhatlm = prefact .* MTlm .* flm;


% Include NQC corrections
% (only for some multipoles)
Modhhatlm(LM2K(2,2),:) = Modhhatlm(LM2K(2,2),:) .* hnqc_22;
Modhhatlm(LM2K(2,1),:) = Modhhatlm(LM2K(2,1),:) .* hnqc_21;
Modhhatlm(LM2K(3,3),:) = Modhhatlm(LM2K(3,3),:) .* hnqc_33;

% FIXME: for all multipoles this becomes
%
%Mod_hhatlm = Mod_hhatlm .* hnqclm;
%


% Compute Newtonian flux
FNewtlm = EOBflmNewt(x,nu,EOBopt);
FNewt22 = FNewtlm( LM2K(2,2) , : );


% Total flux multipoles
Flm = Modhhatlm.^2 .* FNewtlm; 

% Define vOmg, the v we obtain from Omega, and dotvOmg, the time derivative
% of vOmg
%vOmg = Omega.^(1/3);

% Multiply the modes that first enter at 1PN by a constant factor

Flm(LM2K(2,1),:) = Flm(LM2K(2,1),:) .* EOBopt.a1;
Flm(LM2K(3,3),:) = Flm(LM2K(3,3),:) .* EOBopt.a1;
Flm(LM2K(3,1),:) = Flm(LM2K(3,1),:) .* EOBopt.a1;

% Also multiply the modes that first enter at 2PN by a constant factor

Flm(LM2K(3,2),:) = Flm(LM2K(3,2),:) .* EOBopt.a2;
Flm(LM2K(4,4),:) = Flm(LM2K(4,4),:) .* EOBopt.a2;
Flm(LM2K(4,2),:) = Flm(LM2K(4,2),:) .* EOBopt.a2;

% Sum over multipoles and normalize to the 22 Newtonian multipole
hatF = EOBopt.fac .* sum(Flm)./FNewt22;
%hatF = EOBopt.fac .* (1 + EOBopt.Cp .* vOmg.^EOBopt.pow + EOBopt.Ce .* exp(-EOBopt.g ./vOmg)).^2 .* sum(Flm)./FNewt22;

