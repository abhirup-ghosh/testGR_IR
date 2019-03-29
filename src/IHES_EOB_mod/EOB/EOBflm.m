    function [flm, rholm] = EOBflm(x,nu, EOBopt)

%EOBflm Compute the resummed amplitudes in the general nu-dependent case.
%
%
%   [flm, rholm] = EOBflm(x,nu, EOBopt)
%
%        x   :: PN parameter
%        nu  :: symmetric mass ratio. EMRL case is nu=0
%        eps :: eps=1 switches on Iyer-Fujita terms; eps=0 switches them off.
%
%   Reference(s)
%   Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
%   Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
%   Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]

% TODO: this routine requires optimization
% - precompute coefficients c(nu)
% - evaluate efficiently polynomials
% - How does Matlab deal with 1 + small ? (add 1 later?)


LM2K = EOBopt.LM2K;
L    = EOBopt.L;
eps  = strcmp(EOBopt.FIterms,'yes');


% Shorthands
x2  = x.^2;
x3  = x.^3;
x4  = x.^4;
x5  = x.^5;
nu2 = nu^2;
nu3 = nu^3;
nu4 = nu^4;


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

rholm(LM2K(2,2),:) = 1. + (-1.0238095238095237 + 0.6547619047619048*nu)*x + (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2)*x2...
    + x3.*(12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2)...
    + x4.*(-2.4172313935587004 + 4.173242630385488*el2)...
    + x5.*(-30.14143102836864 + 7.916297736025627*el2);
                     
rholm(LM2K(2,1),:) = 1. + (-1.0535714285714286 + 0.27380952380952384*nu)*x ...
    + (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2)*x2 ...
    + x3.*(2.9192806270460925 - 1.019047619047619*el1) ...
    + x4.*(-1.28235780892213 + 1.073639455782313*el1)...
    + eps*x5.*(-3.8466571723355227 + 0.8486467106683944*el1);


% l=3 ------------------------------------------------------------------

rholm(LM2K(3,3),:) = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 ...
    + x3.*(14.10891386831863 - 3.7142857142857144*el3) + x4.*(-6.723375314944128 + 4.333333333333333*el3) ...
    + eps*x5.*(-29.568699895427518 + 6.302092352092352*el3);

rholm(LM2K(3,2),:) = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)./(-1. + 3.*nu)...
    + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)./power(-1. + 3.*nu,2) ...
    + x3.*(6.220997955214429 - 1.6507936507936507*el2) ...
    + eps*x4.*(-3.4527288879001268 + 2.005408583186361*el2);

rholm(LM2K(3,1),:) = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 ...
    + x3.*(1.9098284139598072 - 0.4126984126984127*el1)...
    + x4.*(0.5368150316615179 + 0.2980599647266314*el1)...
    + eps*x5.*(1.4497991763035063 - 0.0058477188106817735*el1);


% l=4 ------------------------------------------------------------------

rholm(LM2K(4,4),:) = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)./(-1. + 3.*nu) ...
    + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/power(-1. + 3.*nu,2)...
    + x3.*(15.108111214795123 - 3.627128427128427*el4)...
    + eps*x4.*(-8.857121657199649 + 4.434988849534304*el4);

rholm(LM2K(4,3),:) = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)./(-1. + 2.*nu) - 0.9783218202252293*x2 ...
    + eps*(x3.*(8.519456157072423 - 2.0402597402597404*el3) ...
    +      x4.*(-5.353216984886716 + 2.5735094451003544*el3));

rholm(LM2K(4,2),:) = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)./(-1. + 3.*nu) ...
    - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/power(-1. + 3.*nu,2)...
    + 4.550378418934105e-12*x3.*(8.48238724511e11 - 1.9927619712e11*el2) ...
    + eps*x4.*(-0.6621921297263365 + 0.787251738160829*el2);

rholm(LM2K(4,1),:) = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)./(-1. + 2.*nu) - 0.36778992787515513*x2 ...
    + x3.*(0.6981550175535535 - 0.2266955266955267*el1) ...
    + eps*x4.*(-0.7931524512893319 + 0.2584672482399755*el1);


% l=5 ------------------------------------------------------------------

rholm(LM2K(5,5),:) = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)./(-1. + 2.*nu) - 1.5749727622804546*x2 ...
    + eps*(x3.*(15.939827047208668 - 3.6037296037296036*el5) ...
    +      x4.*(-10.272578060123237 + 4.500041838503377*el5));

rholm(LM2K(5,4),:) = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)./(1. - 5.*nu + 5.*nu2) ...
    + eps*(- 1.0442142414362194*x2 ...
    +   x3.*(10.252052781721588 - 2.3063869463869464*el4));

rholm(LM2K(5,3),:) = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)./(-1. + 2.*nu) - 0.5788010707241477*x2 ...
    + eps*(x3.*(5.733973288504755 - 1.2973426573426574*el3) ...
    +      x4.*(-1.9573287625526001 + 1.2474448628294783*el3));

rholm(LM2K(5,2),:) = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)./(1. - 5.*nu + 5.*nu2) ...
    + eps*(- 0.4629337197600934*x2 ...
    +      x3.*(2.354458371550237 - 0.5765967365967366*el2));

rholm(LM2K(5,1),:) = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)./(-1. + 2.*nu) - 0.1047896120973044*x2 ...
    + eps*(x3.*(0.642701885362399 - 0.14414918414918415*el1) ...
    +      x4.*(-0.07651588046467575 + 0.11790664036817883*el1));


% l=6 ------------------------------------------------------------------

rholm(LM2K(6,6),:) = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)./(1. - 5.*nu + 5.*nu2) ...
    + eps*(- 1.5543111183867486*x2 + x3.*(16.645950799433503 - 3.6003996003996006*el6));

rholm(LM2K(6,5),:) = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 1.0973940686333457*x2 + x3.*(11.623366217471297 - 2.5002775002775004*el5));

rholm(LM2K(6,4),:) = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)./(1. - 5.*nu + 5.*nu2) ...
    + eps*(- 0.7228451986855349*x2 + x3.*(7.359388663371044 - 1.6001776001776002*el4));

rholm(LM2K(6,3),:) = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 0.5605554442947213*x2 + x3.*(4.002558222882566 - 0.9000999000999002*el3));

rholm(LM2K(6,2),:) = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)./(1. - 5.*nu + 5.*nu2)...
    + eps*( - 0.24797525070634313*x2 + x3.*(1.7942694138754138 - 0.40004440004440006*el2));

rholm(LM2K(6,1),:) = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 0.29175486850885135*x2 + x3.*(0.21653486654395454 - 0.10001110001110002*el1));


% l=7 ------------------------------------------------------------------

rholm(LM2K(7,7),:) = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 1.5418467934923434*x2 + x3.*(17.255875091408523 - 3.6042232277526396*el7));

rholm(LM2K(7,6),:) = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3) ...
    - 1.1403265020692532*eps*x2;

rholm(LM2K(7,5),:) = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 0.8269193364414116*x2 + x3.*(8.750589067052443 - 1.838889401914612*el5));

rholm(LM2K(7,4),:) = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 0.6473746896670599*eps*x2;

rholm(LM2K(7,3),:) = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*(- 0.37187416047628863*x2 + x3.*(3.0835293524055283 - 0.6620001846892604*el3));

rholm(LM2K(7,2),:) = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 0.351319484450667*eps*x2;

rholm(LM2K(7,1),:) = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)./(1. - 4.*nu + 3.*nu2) ...
    + eps*( - 0.1508235111143767*x2 + x3.*(0.2581280702019663 - 0.07355557607658449*el1));


% l=8 ------------------------------------------------------------------

rholm(LM2K(8,8),:) = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 1.5337092502821381*eps*x2;

rholm(LM2K(8,7),:) = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)./(-1. + 6.*nu - 10.*nu2 + 4.*nu3)...
    - 1.175404252991305*eps*x2;

rholm(LM2K(8,6),:) = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 0.9061610303170207*eps*x2;

rholm(LM2K(8,5),:) = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)./(-1. + 6.*nu - 10.*nu2 + 4.*nu3)...
    - 0.7220789990670207*eps*x2;

rholm(LM2K(8,4),:) = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 0.47652059150068155*eps*x2;

rholm(LM2K(8,3),:) = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)./(-1. + 6.*nu - 10.*nu2 + 4.*nu3)...
    - 0.4196774909106648*eps*x2;

rholm(LM2K(8,2),:) = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)./(-1. + 7.*nu - 14.*nu2 + 7.*nu3)...
    - 0.2261796441029474*eps*x2;

rholm(LM2K(8,1),:) = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)./(-1. + 6.*nu - 10.*nu2 + 4.*nu3)...
    - 0.26842133517043704*eps*x2;


% Compute the f_lm = rho_lm^l
%pL  = repmat(L',1,length(x));
%flm = rholm.^pL;
for k=1:kmax
    flm(k,:) = rholm(k,:).^L(k);
end

% TEST: make parallel these loops?
%matlabpool
%parfor k=1:kmax
%    flm(k,:) = rholm(k,:).^L(k);
%end








