function [LM2K, L,M] = Multipidx(lmax)

%MULTIPIDX Compute maps between multipolar indexes.
%
%   [LM2K, L,M] = Multipidx(lmax) return sparse matrix LM2K mapping
%   (l,m)->k, and vectors L,M mapping k->l,m
%
%   NOTE: Consider only multipoles m>0.
%


k = 0;
for l=2:lmax
    for m=1:l
        k = k + 1;
        
        % map (l,m)->k
        LM2K(l,m) = k;
        
        % map k->l,m
        L(k) = l;
        M(k) = m;
        
    end
end

LM2K = sparse(LM2K);





