function varargout = EOBNQCn(r,prstar, Omega,ddotr)

%EOBNQCn Compute n's quantities for NQC.
%
%   [n1,n2,n3,n4,n5,n6] = EOBNQCn(r,prstar, Omega,ddotr);
%
%   NQC = EOBNQCn(r,prstar, Omega,ddotr);
%


% NQC corrections to the modulus
n1 = (prstar./(r.*Omega)).^2;  % [pr*/(r Omg)]^2
n2 = ddotr./(r.*Omega.^2);     % ddot{r}./(r Omg^2)
n3 = n1.*prstar.^2;            % [pr*/(r Omg)]^2 *(pr*)^2


% NQC corrections to the phase
n4 = prstar./(r.*Omega);       %  pr*/(r Omg)
n5 = n4.*(r.*Omega).^2;        %  (pr*)*(r Omg)
n6 = n5.*prstar.^2;            %  (pr*^3)*(r Omg)


varargout = SetVOutput( nargout, n1,n2,n3,n4,n5,n6 );
