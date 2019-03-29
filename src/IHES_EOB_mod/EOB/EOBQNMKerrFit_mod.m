function sigma = EOBQNMKerrFit_mod(nqnm,lmax, L,M,abh, datadir, iQNMfac)

%EOBQNMKerrFit Fit to the perturbative QNM sigma's.
%
%   sigma = EOBQNMKerrFit_mod(nqnm, L,M, abh, datadir, iQNMfac) return complex
%   frequencies of Kerr QNM fitted to perturbative data, for a given values
%   of angular momentum abh, with the imaginary part of the QNM multiplied by
%   iQNMfac. It considers only positive frequency QNMs.
%
%   Reference(s)
%    Berti et al Class. Quant. Grav. 26, 163001 (2009)
%


if nqnm > 8
    error('no data available for n>8, use a smaller number of QNM (%d)',nqnm);
end
if lmax > 7
    warning('no data available for l>7, sigma''s are set to 0');
end


% Allocate memory
kmax = length(L);
sigma = zeros(kmax,nqnm);


% Convert indexes
idx = (L<=lmax);
M = M(idx);
L = L(idx);
L = strread(num2str(L),'%s');
M = strread(num2str(M),'%s');


% Interpolate to abh
order = 6 ;
for n=1:nqnm
    
    filename   = strcat(datadir,'/l',L,'/n',num2str(n),'l',L,'m',M,'.dat');       
    sigmak     = cellfun(@(x) InterpData_mod(x, abh, order, iQNMfac), filename);      
    sigma(idx,n) = sigmak;

end



function sigma = InterpData_mod(filename, a, order, iQNMfac)

%InterpData Load file and interpolate to a, multiplying the imaginary part of the QNM by iQNMfac.
%
%   sigma = InterpData_mod(filename, a, order, iQNMfac)
%

y = load(filename);

omega = LagInt1d( order, y(:,1), y(:,2), a );
beta  = LagInt1d( order, y(:,1), y(:,3), a );

sigma =  (1i*omega - iQNMfac*beta);