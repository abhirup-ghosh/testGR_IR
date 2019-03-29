function varargout = EOBhlmQNMMatch(Omega,t, sigma, nqnm, dt, wav )

%EOBhlmQNMMatch Perform the matching to QNM multipolar wave.
%   Uses a nonlinear fit
%
%   [psilm,dpsilm,omglm] = EOBhlmQNMMatch(Omega,t, sigma, nqnm, dt, wav )
%
%   WAVE = EOBhlm( .... ) return a structure with the wave
%


Omegamax_interpolate = 1; % hardcoded option FIXME

% Unpack wav
psilm  = wav.psilm;
dpsilm = wav.dlm;
omglm  = wav.omglm;

psiqnm  = 0*psilm;
dpsiqnm = psiqnm;
omgqnm  = psiqnm;


% Compute max(Omega)
[Omgmax,jmax]  = max(Omega);
tmax = t(jmax);
if Omegamax_interpolate    
    jIdx = jmax-3:jmax+3;
    [tmax,Omgmax] = FindMax( t(jIdx), Omega(jIdx), 2 );
end
tmatch = tmax;


% Grid around t=0 (tmatch)
%{ 
% matching grid with #pts > #QNMs
n = 2*sum(nqnm);      % no of real QNM coefficients (incl. pos and neg QNM)
n = n + 2 + mod(n,2); % add some degree of freedom, n is even
tgrid = [-n/2:1:n/2]*dt;
%}
% #pts = #QNMs
n=sum(nqnm);
if mod(n,2)==0
    % even number of QNMs
    tgrid = (-n:1:n-1)*dt;
else
    % odd number of QNMs
    tgrid = (-(round(n/2)-1):1:round(n/2)-1)*dt;
end


% Interpolate psi on new grid, shift tmatch to t=0
if n>6
    order=6;
elseif n>4
    order=4;
    warning('reducing intepolation order to %d',order)
else
    order=2;
    warning('reducing intepolation order to %d',order)
end


    

% Fit multipoles
[kmax nnn] = size(sigma);
if nnn~=sum(nqnm)
    error('dimension problem')
end
knnz = 1:kmax;
knnz = knnz(all(sigma~=0,2)); % index of nnz multipoles

for k=knnz
    
    % Interpolate    
    tmpr = LagInt1d( order, t-tmatch, real(psilm(k,:)), tgrid );
    tmpi = LagInt1d( order, t-tmatch, imag(psilm(k,:)), tgrid );        
    
    psigrid = tmpr + 1i*tmpi;          
        
    % Fit with complex coefficients    
    Ck = ones(nnn,1) + 1i*ones(nnn,1);               
    Ck = EOBnlinfit(tgrid,psigrid, @(b,x) EOBQNMTemplate(b,x, sigma(k,:)), Ck); 
    
    % Calculate the QNM wave on t-tmatch
    [psi dpsi] = EOBQNMTemplate(Ck,t-tmatch, sigma(k,:));    
    
    psiqnm(k,:)  = psi;
    dpsiqnm(k,:) = dpsi;
    omgqnm(k,:)  = -imag(dpsi./psi);          

end


% Do the matching at t=tmatch
tIdx = (t<=tmatch); 
psiqnm(:,tIdx)  = 0;
dpsiqnm(:,tIdx) = 0;
omgqnm(:,tIdx)  = 0;            
         
tIdx = (t>tmatch); 
psilm(:,tIdx)  = 0;
dpsilm(:,tIdx) = 0;
omglm(:,tIdx)  = 0;

psilm  = psilm  + psiqnm;
dpsilm = dpsilm + dpsiqnm;
omglm  = omglm  + omgqnm;


% Finalize
varargout = SetVOutput( nargout, psilm,dpsilm,omglm );







