function mTlm = EOBModTail( L,M, w )

%EOBModTail Compute the modulus of the complex tail factor |T_lm|.
%
%   mTlm = EOBModTail( L, M, Omega ) return |T_lm(mOmega)| given the vectors
%   of multipolar linear indexes L(k) and M(k) and value(s) of Omega.
%

% FIXME: this routine should be optimized
% - vectorize the k-loop
% - vectorize the i-loop (during dynamics w is scalar, so this should
% be less important)


kmax  = length(L);
nx    = length(w);
mTlm  = zeros(kmax,nx);


L    = reshape(L,kmax,1);
M    = reshape(M,kmax,1);
w    = reshape(w,1,nx);
hatk = M * w;
x2   = (-2.*hatk).^2;
y    = 4*pi*hatk;

% v2 (some optimization wrt v1)

oofact2 = ( 1./(factorial(L) ) ).^2;

for i=1:nx
    
    yk  = y(:,i);    

    for k=1:kmax
        s2  = [1:L(k)].^2;
        tmp(k) = prod(s2+x2(k,i));                    
    end
    
    mTlm(:,i) = oofact2 .* yk./(1 - exp(-yk)).*tmp.';
    
end

mTlm = sqrt(mTlm);


% here below old versions and tests

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
% v1 (rm inner s-loop)
for i=1:nx  
    for k=1:kmax
            
        x2ki = x2(k,i);
        yki  = y(k,i);
        
        s2  = [1:L(k)].^2;
        tmp = prod(s2+x2ki);
    
        mTlm2 = 1./(factorial(L(k)).^2) .* yki./(1 - exp(-yki)).*tmp;            
        mTlm(k,i)  = sqrt(mTlm2);    
        
    end                
end
%}

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 

% v0 (unvectorized)
%for k=1:kmax
for i=1:nx         
  for k=1:kmax
            
            hatk = M(k)*w(i);
            x    = -2.*hatk;
            x2   = x.^2;
            y    = 4*pi*hatk;
            prod = 1;

            for s=1:L(k)
                s2   = s.^2;
                prod = prod.*(s2+x2);
                %fprintf(' %d  -> %g %g -> %g\n',s,s2,x2, prod);                
            end
            
            %fprintf('%d (%d %d) %d -> %g -> %g\n',k,M(k),L(k),i, hatk,prod)
            %fprintf('%d -> y = %g\n',k,y./(1 - exp(-y)))
            
            mTlm2 = 1./(factorial(L(k)).^2) .* y./(1 - exp(-y)).*prod;            
            mTlm(k,i)  = sqrt(mTlm2);    
            
    end    
end

%}

