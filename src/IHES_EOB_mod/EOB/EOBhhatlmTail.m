function tlm = EOBhhatlmTail( L,M, w,e, bphys)

%EOBhhatlmTail Computes the tail contribution to the resummed wave.
%
%   tlm = EOBTail(L,M, Omega,E, bphys) 
%
%   Damour, Iyer & Nagar, PRD 79, 064004 (2009)

kmax  = length(L);
nx    = length(w);

L    = reshape(L,kmax,1);
M    = reshape(M,kmax,1);
w    = reshape(w,1,nx);
e    = reshape(e,1,nx);

k    = M * w;
hatk = M * (w.*e);

L    = (ones(nx,1)*L')'; 

ratio = GammaComplex(L + 1-2*1i*hatk)./GammaComplex(L+1);
% MEX-file version (slower, Matlab optimized for complex vars)
%ratio = GammaComplexm(L + 1-2*1i*hatk)./GammaComplexm(L+1);


%tlm = ratio .* exp(pi.*hatk) .*exp(2.*1i*hatk.*log(2.d0*k.*bphys)) ;
tlm   = ratio .* exp( pi.*hatk + 2.*1i*hatk.*log(2.*k.*bphys) );




