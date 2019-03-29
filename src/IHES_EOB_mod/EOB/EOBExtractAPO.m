function varargout = EOBExtractAPO(hlm,t)

%EOBExtractAPO Compute amplitude, phase and frequency of waveform.
%
%   [Alm,philm, omglm,domglm,d2omglm, dlm,d2lm] = EOBExtractAPO(hlm,t)
%   return amplitude, phase, frequency and derivatives and first and second
%   complex derivative of the complex multipolar waveform. 
%   Assume size(hlm) = (no multipoles) x (time).
%
%   wav = EOBExtractAPO(hlm,t) return a structure with fields as above
%


% Assume size = (no multipoles) x (time)
[kmax nx] = size(hlm); 


% Amplitude
Alm    = abs(hlm);

% Phase
philm   = - unwrap(angle(hlm),[],2);
omglm   = 0 *hlm;
domglm  = omglm;
d2omglm = omglm;
dlm     = omglm;
d2lm    = omglm;


% Frequency and derivatives, loop on nonzeros multipoles only 
% Drvts calculation requires optimization 
knnz = 1:kmax;
knnz = knnz(all(hlm~=0,2)); % index of nnz multipoles
%tic
for k=knnz

  omglm(k,:)   = FDdrvt( philm(k,:) , t, 4);    
  
  %{
  domglm(k,:)  = FDdrvt( omglm(k,:) , t, 4);
  d2omglm(k,:) = FDdrvt( domglm(k,:), t, 4);
  dlm(k,:)     = FDdrvt( hlm(k,:)   , t, 4);
  d2lm(k,:)    = FDdrvt( dlm(k,:)   , t, 4);
  %}
  % save computation:
  [domglm(k,:),d2omglm(k,:)] = FDdrvt( omglm(k,:) , t, 4);   
  [dlm(k,:),d2lm(k,:)]       = FDdrvt( hlm(k,:)   , t, 4);  
  
end
%toc


% Finalize
varargout = SetVOutput( nargout, Alm,philm, omglm,domglm,d2omglm, dlm,d2lm );