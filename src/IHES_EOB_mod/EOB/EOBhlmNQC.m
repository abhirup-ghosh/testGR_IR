function psilmnqc = EOBhlmNQC( nu, r,prstar, Omega,ddotr, a,b, EOBopt )

%EOBhNQC Calculate NQC corrections to the RWZ multipolar waveform.
%
%   psilmnqc = EOBhNQC( nu, r,prstar, Omega,ddotr, a,b, EOBopt )
%


% Compute n
NQC = EOBNQCn(r,prstar, Omega,ddotr);
n = struct2cell(NQC);


% Reshape things
kmax = length(a{1});
nx   = length(n{1});
n    = cellfun(@(x) reshape(x,1,nx),n,'UniformOutput',false);
a    = cellfun(@(x) reshape(x,kmax,1),a,'UniformOutput',false);
b    = cellfun(@(x) reshape(x,kmax,1),b,'UniformOutput',false);


% NQC multipolar correction factor
psilmnqc = (1 + a{1}*n{1} + a{2}*n{2} + a{3}*n{3})...
    .*exp( 1i*(b{1}*n{4} + b{2}*n{5} + b{3}*n{6}) );




