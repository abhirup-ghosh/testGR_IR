function [flm, rholm] = EOBflm0(x, EOBopt)

%EOBflm Compute the resummed amplitudes in the nu=0 case.
%
%
%   [flm, rholm] = EOBflm0(x, EOBopt)
%
%        x   :: PN parameter
%        eps :: eps=1 switches on Iyer-Fujita terms; eps=0 switches them off.
%

% TODO: this routine requires optimization
% - precompute coefficients c(nu)
% - evaluate efficiently polynomials


LM2K = EOBopt.LM2K;
L    = EOBopt.L;
eps  = strcmp(EOBopt.FIterms,'yes');


% Shorthands
x2  = x.^2;
x3  = x.^3;
x4  = x.^4;
x5  = x.^5;


% Compute EulerLogs
el1 = Eulerlog(x,1);
el2 = Eulerlog(x,2);
el3 = Eulerlog(x,3);
el4 = Eulerlog(x,4);
el5 = Eulerlog(x,5);
el6 = Eulerlog(x,6);
el7 = Eulerlog(x,7);


% Compute rho_{lm}
kmax  = length(L);
rholm = zeros(kmax,length(x));
flm   = rholm;


% l=2 ------------------------------------------------------------------

rholm(LM2K(2,2),:) = 1. - (43*x)/42. - (20555*x2)/10584. + (x3.*(12455352904 - 3986357760*el2))/9.779616e8 ...
    + x4.*(-2.4172313935587004 + (9202*el2)/2205.) + x5.*(-30.14143102836864 + (439877*el2)/55566.);

rholm(LM2K(2,1),:)= 1. - (59*x)/56. - (47009*x2)/56448. + x3.*(2.9192806270460925 - (107*el1)/105.) ...
    + eps.*x5.*(-3.8466571723355227 + (5029963*el1)/5.92704e6) + x4.*(-1.28235780892213 + (6313*el1)/5880.);

% l=3 ------------------------------------------------------------------

rholm(LM2K(3,3),:) = 1. - (7*x)/6. - (6719*x2)/3960. + x3.*(14.10891386831863 - (26*el3)/7.) ...
    + x4.*(-6.723375314944128 + (13*el3)/3.) ...
    + eps.*x5.*(-29.568699895427518 + (87347*el3)/13860.);

rholm(LM2K(3,2),:) = 1. - (164*x)/135. - (180566*x2)/200475. + x3.*(6.220997955214429 - (104*el2)/63.) ...
    + eps*x4.*(-3.4527288879001268 + (17056*el2)/8505.);

rholm(LM2K(3,1),:) = 1. - (13*x)/18. + (101*x2)/7128. ...
    +  x3.*(1.9098284139598072 - (26*el1)/63.) + x4.*(0.5368150316615179 + (169*el1)/567.)...
    + eps*x5.*(1.4497991763035063 - (1313*el1)/224532.);

% l=4 ------------------------------------------------------------------

rholm(LM2K(4,4),:) = 1. - (269*x)/220. - (14210377*x2)/8.8088e6 + x3.*(15.108111214795123 - (12568*el4)/3465.) ...
    + eps*x4.*(-8.857121657199649 + (845198*el4)/190575.);


rholm(LM2K(4,3),:) = 1. - (111*x)/88. - (6894273*x2)/7.04704e6 ...
    + eps*( x3.*(8.519456157072423 - (1571*el3)/770.)...
    +       x4.*(-5.353216984886716 + (174381*el3)/67760.));

rholm(LM2K(4,2),:) = 1. - (191*x)/220. - (3190529*x2)/8.8088e6 + (x3.*(848238724511 - 199276197120*el2))/2.197619424e11 ...
    + eps*x4.*(-0.6621921297263365 + (300061*el2)/381150.);

rholm(LM2K(4,1),:) = 1. - (301*x)/264. - (7775491*x2)/2.114112e7 ...
    + eps*(x3.*(0.6981550175535535 - (1571*el1)/6930.) ...
    + x4.*(-0.7931524512893319 + (67553*el1)/261360.));

% l=5 ------------------------------------------------------------------

rholm(LM2K(5,5),:) = 1. - (487*x)/390. - (3353747*x2)/2.1294e6 ...
    + eps*(x3.*(15.939827047208668 - (1546*el5)/429.) ...
    + x4.*(-10.272578060123237 + (376451*el5)/83655.));

rholm(LM2K(5,4),:) = 1. - (2908*x)/2275. ...
    + eps*(- (16213384*x2)/1.5526875e7 + x3.*(10.252052781721588 - (24736*el4)/10725.));

rholm(LM2K(5,3),:) = 1. - (25*x)/26. - (410833*x2)/709800. ...
    + eps*( x3.*(5.733973288504755 - (4638*el3)/3575.) ...
    + x4.*(-1.9573287625526001 + (2319*el3)/1859.));

rholm(LM2K(5,2),:) = 1. - (2638*x)/2275. ...
    + eps*(- (7187914*x2)/1.5526875e7 + x3.*(2.354458371550237 - (6184*el2)/10725.));

rholm(LM2K(5,1),:) = 1. - (319*x)/390. - (31877*x2)/304200. ...
    + eps*( x3.*(0.642701885362399 - (1546*el1)/10725.) ...
    + x4.*(-0.07651588046467575 + (22417*el1)/190125.));

% l=6 ------------------------------------------------------------------

rholm(LM2K(6,6),:) = 1. - (53*x)/42.   - eps*(1025435*x2)/659736.       + eps*x3.*(16.645950799433503 - (3604*el6)/1001.);

rholm(LM2K(6,5),:) = 1. - (185*x)/144. - eps*(59574065*x2)/5.4286848e7  + eps*x3.*(11.623366217471297 - (22525*el5)/9009.);

rholm(LM2K(6,4),:) = 1. - (43*x)/42.   - eps*(476887*x2)/659736.        + eps*x3.*(7.359388663371044 - (14416*el4)/9009.);

rholm(LM2K(6,3),:) = 1. - (169*x)/144. - eps*(152153941*x2)/2.7143424e8 + eps*x3.*(4.002558222882566 - (901*el3)/1001.);

rholm(LM2K(6,2),:) = 1. - (37*x)/42.   - eps*(817991*x2)/3.29868e6      + eps*x3.*(1.7942694138754138 - (3604*el2)/9009.);

rholm(LM2K(6,1),:) = 1. - (161*x)/144. - eps*(79192261*x2)/2.7143424e8  + eps*x3.*(0.21653486654395454 - (901*el1)/9009.);

% l=7 ------------------------------------------------------------------

rholm(LM2K(7,7),:) = 1. - (151*x)/119.   - eps*(32358125*x2)/2.0986602e7 + eps*x3.*(17.255875091408523 - (11948*el7)/3315.);

rholm(LM2K(7,6),:) = 1. - (1072*x)/833.  - eps*(195441224*x2)/1.71390583e8;

rholm(LM2K(7,5),:) = 1. - (127*x)/119.   - eps*(17354227*x2)/2.0986602e7 + eps*x3.*(8.750589067052443 - (59740*el5)/32487.);

rholm(LM2K(7,4),:) = 1. - (8878*x)/7497. - eps*(2995755988*x2)/4.627545741e9;

rholm(LM2K(7,3),:) = 1. - (111*x)/119.   - eps*(7804375*x2)/2.0986602e7 + eps*x3.*(3.0835293524055283 - (35844*el3)/54145.);

rholm(LM2K(7,2),:) = 1. - (8416*x)/7497. - eps*(1625746984*x2)/4.627545741e9;

rholm(LM2K(7,1),:) = 1. - (103*x)/119.   - eps*(1055091*x2)/6.995534e6 + eps*x3.*(0.2581280702019663 - (11948*el1)/162435.);

% l=8 ------------------------------------------------------------------

rholm(LM2K(8,8),:) = 1. - (1741*x)/1368. - eps*(9567401*power(x,2))/6.23808e6;

rholm(LM2K(8,7),:) = 1. - (3913*x)/3040. - eps*(195527087*power(x,2))/1.663488e8;

rholm(LM2K(8,6),:) = 1. - (167*x)/152.   - eps*(376847*power(x,2))/415872.;

rholm(LM2K(8,5),:) = 1. - (725*x)/608.   - eps*(4804679*power(x,2))/6.653952e6;

rholm(LM2K(8,4),:) = 1. - (1333*x)/1368. - eps*(1387201*power(x,2))/2.911104e6;

rholm(LM2K(8,3),:) = 1. - (3433*x)/3040. - eps*(7756983*power(x,2))/1.84832e7;

rholm(LM2K(8,2),:) = 1. - (1231*x)/1368. - eps*(9876487*power(x,2))/4.366656e7;

rholm(LM2K(8,1),:) = 1. - (3337*x)/3040. - eps*(44651567*power(x,2))/1.663488e8;


% Compute the f_lm = rho_lm^l
%pL  = repmat(L',1,length(x));
%flm = rholm.^pL;
for k=1:kmax
    flm(k,:) = rholm(k,:).^L(k);
end


