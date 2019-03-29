function [a,b] = EOBNQCabFind( t,Omega,n, ...
    Alm,omglm,domglm,d2omglm, NRPTS )

%EOBNQCabFind Compute a's and b's
%
%   [a,b] = EOBNQCabFind( t,Omega,n,  Alm,omglm,domglm,d2omglm, NRPTS )
%
%



warning('this routine must be tested');


Omegamax_interpolate = 0; % hardcoded option FIXME


% Allocate memory
[kmax nt] = size(Alm);
a{1} = zeros(kmax,1);
a{2} = a{1};
a{3} = a{1};
b{1} = a{1};
b{2} = a{1};
b{3} = a{1};


% Numerical derivatives of the phase n_i's
% Drvts calculation requires optimization (TODO)
% TODO: there're 3 drvts: use 6th order ?

dn  = n; % dn's are defined of size(n) for simplicity
d2n = n; % so to use indexes dn{4}
d3n = n; % however dn's{1...3} not needed

%tic
[dn{4}, d2n{4}] = FDdrvt(n{4},t,4);
[dn{5}, d2n{5}] = FDdrvt(n{5},t,4);
[dn{6}, d2n{6}] = FDdrvt(n{6},t,4);
d3n{4}          = FDdrvt(d2n{4},t,4);
d3n{5}          = FDdrvt(d2n{5},t,4);
d3n{6}          = FDdrvt(d2n{6},t,4);
%toc


% Determine the maximum of the orbital frequency
% and compute dn's at the maximum (overwrite)
% do not overwrite n : needed below !
[Omegamax,jmax]  = max(Omega);

if Omegamax_interpolate
    
    jIdx = jmax-3:jmax+3;
    [tmax,Omegamax] = FindMax( t(jIdx), Omega(jIdx), 2 );
    
    %maxn   = cellfun(@(x) LagInt1d( 6, t, x, tmax ),n);
    dn  = cellfun(@(x) LagInt1d( 6, t, x, tmax ),dn);
    d2n = cellfun(@(x) LagInt1d( 6, t, x, tmax ),d2n);
    d3n = cellfun(@(x) LagInt1d( 6, t, x, tmax ),d3n);
    
else
    
    tmax = t(jmax);
    
    %maxn   = cellfun(@(x) x(jmax),n);
    dn  = cellfun(@(x) x(jmax),dn);
    d2n = cellfun(@(x) x(jmax),d2n);
    d3n = cellfun(@(x) x(jmax),d3n);
    
end


% Matrix N is independent on multipolar index => pre-compute it
NN = [dn(4) dn(5) dn(6); d2n(4) d2n(5) d2n(6); d3n(4) d3n(5) d3n(6)];


% Compute a's and b's for each multipole
ktmp = 1:kmax;
knnz = intersect( ktmp(all(Alm~=0,2)) , ktmp(all(NRPTS~=0,2)) );
for k=knnz
    
    % NR pts for kth multipole
    NRAk     = NRPTS(k,3);
    NRdAk    = NRPTS(k,4);
    NRd2Ak   = NRPTS(k,5);    
    NRomgk   = NRPTS(k,6);
    NRdomgk  = NRPTS(k,7);
    NRd2omgk = NRPTS(k,8);    
        
    % Solve linear systems    

    % Matrix M
    m11 = n{1}.'.*Alm(k,:);
    m12 = n{2}.'.*Alm(k,:);
    m13 = n{3}.'.*Alm(k,:);    
    
    m21 = FDdrvt(m11,t,4);
    m22 = FDdrvt(m12,t,4);
    m23 = FDdrvt(m13,t,4);    
    
    m31 = FDdrvt(m21,t,4);
    m32 = FDdrvt(m22,t,4);
    m33 = FDdrvt(m23,t,4);    
    
    % Vectors
    p1tmp = Alm(k,:);
    p2tmp = FDdrvt(p1tmp,t,4);
    p3tmp = FDdrvt(p2tmp,t,4);    
    
    q1tmp = omglm(k,:);
    q2tmp = domglm(k,:);
    q3tmp = d2omglm(k,:);
        
    % Compute at max(Omega)    
    % overwrite tmp vars
    if Omegamax_interpolate
        
        m11 = LagInt1d( 6, t, m11, tmax );
        m12 = LagInt1d( 6, t, m12, tmax );
        m13 = LagInt1d( 6, t, m13, tmax );        
        
        m21 = LagInt1d( 6, t, m21, tmax );
        m22 = LagInt1d( 6, t, m22, tmax );
        m23 = LagInt1d( 6, t, m23, tmax );        
        
        m31 = LagInt1d( 6, t, m31, tmax );
        m32 = LagInt1d( 6, t, m32, tmax );
        m33 = LagInt1d( 6, t, m33, tmax );        
        
        p1tmp = LagInt1d( 6, t, p1tmp, tmax );
        p2tmp = LagInt1d( 6, t, p2tmp, tmax );
        p3tmp = LagInt1d( 6, t, p3tmp, tmax );
        
        q1tmp = LagInt1d( 6, t, q1tmp, tmax );
        q2tmp = LagInt1d( 6, t, q2tmp, tmax );
        q3tmp = LagInt1d( 6, t, q3tmp, tmax );
        
    else
        
        m11 = m11(jmax);
        m12 = m12(jmax);
        m13 = m13(jmax);        
        
        m21 = m21(jmax);
        m22 = m22(jmax);
        m23 = m23(jmax);        

        m31 = m31(jmax);
        m32 = m32(jmax);
        m33 = m33(jmax);        
        
        p1tmp = p1tmp(jmax);
        p2tmp = p2tmp(jmax);
        p3tmp = p3tmp(jmax);
        
        q1tmp = q1tmp(jmax);
        q2tmp = q2tmp(jmax);
        q3tmp = q3tmp(jmax);
        
    end
        
    % Computation of ai's
    P   = [NRAk - p1tmp; NRdAk - p2tmp; NRd2Ak - p3tmp];
    MM  = [m11 m12 m13; m21 m22 m23;  m31 m32 m33];   
    ak  = MM\P;
        
    a{1}(k) = ak(1);
    a{2}(k) = ak(2);
    a{3}(k) = ak(3);
    
    % Computation of bi's
    Q  = [q1tmp - NRomgk; q2tmp  - NRdomgk; q3tmp - NRd2omgk];
    bk = NN\Q;            
        
    b{1}(k) = bk(1);
    b{2}(k) = bk(2);
    b{3}(k) = bk(3);
    
end


a{4} = 0*a{1};
b{4} =   a{4};







    
    