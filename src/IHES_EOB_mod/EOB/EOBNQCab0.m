function [a,b]=EOBNQCab0(EOBopt)

%EOBNQCabFit Coefficients for NQC in the test-mass limit
%
%   [a,b] = EOBInitNQCFit(EOBopt )
%
%   multipoles up to l=7 are included.

L    = EOBopt.L;
LM2K = EOBopt.LM2K;


% Shorthands
k22 = LM2K(2,2);
k21 = LM2K(2,1);
k33 = LM2K(3,3);
k32 = LM2K(3,2);
k31 = LM2K(3,1);
k44 = LM2K(4,4);
k43 = LM2K(4,3);
k42 = LM2K(4,2);
k41 = LM2K(4,1);
k55 = LM2K(5,5);
k54 = LM2K(5,4);
k53 = LM2K(5,3);
k52 = LM2K(5,2);
k51 = LM2K(5,1);
k66 = LM2K(6,6);
k65 = LM2K(6,5);
k64 = LM2K(6,4);
k63 = LM2K(6,3);
k62 = LM2K(6,2);
k61 = LM2K(6,1);
k77 = LM2K(7,7);
k76 = LM2K(7,6);
k75 = LM2K(7,5);
k74 = LM2K(7,4);
k73 = LM2K(7,3);
k72 = LM2K(7,2);
k71 = LM2K(7,1);

% Allocate memory
kmax = length(L);
a{1} = zeros(kmax,1);
a{2} = a{1};
a{3} = a{1};
b{1} = a{1};
b{2} = a{1};
b{3} = a{1};


% l=2 -------------------------------------------------------------------

a{1}(k22) = 0.114911894068123; 
a{2}(k22) = 1.055397427331828; 
a{3}(k22) = 0.004783628924899;

b{1}(k22) = 0.130101234415829;  
b{2}(k22) = 0.686482054008241;  
b{3}(k22) = -3.71565356628882; 


a{1}(k21) =  0.130980389482817;
a{2}(k21) =  0.116834175645851; 
a{3}(k21) =  0.192309022550489; 

b{1}(k21) = 0.247438431691121; 
b{2}(k21) = 1.349181281053549; 
b{3}(k21) = -5.338284611946910e-04; 


% l=3 -------------------------------------------------------------------

a{1}(k33) =  0.184346606177275;
a{2}(k33) =  0.985254785576045; 
a{3}(k33) =  0.430859956818454; 

b{1}(k33) =  0.237663568616349; 
b{2}(k33) =  0.613227566892548; 
b{3}(k33) = -3.857440858170929;


a{1}(k32) = 0.123584920956022;
a{2}(k32) = 0.718239011763655; 
a{3}(k32) = 0.435676621761100; 

b{1}(k32) = 0.327542710052521; 
b{2}(k32) = 1.545493532207134; 
b{3}(k32) = 0.026886292653762;


a{1}(k31) =  3.120025486398302;	  
a{2}(k31) = -6.364866203100245;
a{3}(k31) =  0.822824316192125;

b{1}(k31) = 0.432624202212259; 
b{2}(k31) = 3.581267677294043; 
b{3}(k31) = 0.142892516331895;


% l=4 ------------------------------------------------------------------

a{1}(k44) =  0.208068168798727;
a{2}(k44) =  1.109711167298422; 
a{3}(k44) =  1.104100835589372;

b{1}(k44) =  0.324601353901579; 
b{2}(k44) =  0.629835495958264; 
b{3}(k44) = -4.591953989712212;


a{1}(k43) = 0.130073218124257;
a{2}(k43) = 0.894269912607387;
a{3}(k43) = 0.973650412443810;

b{1}(k43) = 0.407457220010151;
b{2}(k43) = 1.073851041535936;
b{3}(k43) = 0.479553412802911;


a{1}(k42) =  1.769972705392476;
a{2}(k42) = -4.236446654920220;
a{3}(k42) =  3.958797395614128;

b{1}(k42) = 0.513959717911504;
b{2}(k42) = 2.995785971549151;
b{3}(k42) = 1.379786531806007;


a{1}(k41) =   4.293804798214098;
a{2}(k41) = -19.887254640783745;
a{3}(k41) =   8.528746919511963; 

b{1}(k41) = 0.534371036659466;
b{2}(k41) = 5.508805361152453;
b{3}(k41) = 4.970542502941439;


% l=5 ------------------------------------------------------------------

a{1}(k55) =  0.132200777268344;
a{2}(k55) =  0.863959711652405; 
a{3}(k55) =  2.452369193415340;

b{1}(k55) =  0.371519386654171;
b{2}(k55) = -2.931907164419524;
b{3}(k55) =  3.416763316233754;


a{1}(k54) = 0.094930427236047;
a{2}(k54) = 0.931954988025597;
a{3}(k54) = 1.919135924406371;

b{1}(k54) =  0.480718530745968;
b{2}(k54) = -0.950874234320440;
b{3}(k54) =  4.323526525478892;


a{1}(k53) =  1.098761342838276;
a{2}(k53) = -4.198817007149198;
a{3}(k53) =  7.710773597932726;

b{1}(k53) = 0.596095143384297;
b{2}(k53) = 1.616068028965650;
b{3}(k53) = 3.959506825236499;


a{1}(k52) =  1.760534186892328;
a{2}(k52) = -8.485273555600328;
a{3}(k52) =  7.589878371399100;

b{1}(k52) = 0.630805250280266;
b{2}(k52) = 3.934919256773048;
b{3}(k52) = 5.556481890964092;


a{1}(k51) = 23.605092474870169;
a{2}(k51) = -87.647612083213232;
a{3}(k51) = 1.498260127945150e+02;

b{1}(k51) = 0.870401381376948;
b{2}(k51) = 6.012078310787134;
b{3}(k51) = 2.159540510773954;


% l=6 ------------------------------------------------------------------

a{1}(k66) = -0.054709879421598;
a{2}(k66) =  0.407148031989311;
a{3}(k66) =  4.531439405040409;

b{1}(k66) =  0.443359356208811;
b{2}(k66) = -3.841714473752655;
b{3}(k66) =  4.605140566125312;


a{1}(k65) = -0.032835441650184; 
a{2}(k65) =  0.591972168017137; 
a{3}(k65) =  3.557530507529137; 

b{1}(k65) =  0.553925226968460; 
b{2}(k65) = -1.626151369204670; 
b{3}(k65) =  4.961600468755351;


a{1}(k64) =  0.356077132193026;
a{2}(k64) = -5.308875788470417; 
a{3}(k64) = 12.983788587706304; 

b{1}(k64) = 0.664966745233224; 
b{2}(k64) = 0.944067590854806; 
b{3}(k64) = 4.965395482449235; 


a{1}(k63) =  0.771704971898766;
a{2}(k63) = -6.730129617199153; 
a{3}(k63) =  10.648855998537769;

b{1}(k63) = 0.756568491277079; 
b{2}(k63) = 3.308816540717364; 
b{3}(k63) = 5.126120203788691; 


a{1}(k62) =  2.283527339342328;
a{2}(k62) = -43.932937814815901; 
a{3}(k62) =  91.026341209451161; 

b{1}(k62) = 0.859439021445532; 
b{2}(k62) = 5.926023551623834; 
b{3}(k62) = 4.222778113868663; 


a{1}(k61) =  42.158896041606852;
a{2}(k61) = -2.744995541393576e+02; 
a{3}(k61) =  5.140989295851671e+02;

b{1}(k61) =  2.258106462660969;
b{2}(k61) =  14.228698864841554;
b{3}(k61) = -36.541201521600250;


% l=7 ------------------------------------------------------------------

a{1}(k77) = -0.398350811783813;
a{2}(k77) = -0.578800426315505; 
a{3}(k77) =  7.667807359334019; 

b{1}(k77) =  0.510644040798544; 
b{2}(k77) = -4.525099915864052;
b{3}(k77) =  5.294705255493486;


a{1}(k76) = -0.494732510643102;
a{2}(k76) = -4.769705227136186;
a{3}(k76) =  8.327199181852539;

b{1}(k76) = 0.610103684022684;
b{2}(k76) = -2.34616792314864; 
b{3}(k76) = 6.016084938993719;


a{1}(k75) = -0.644041179825725;
a{2}(k75) = -7.563087510368737;
a{3}(k75) = 20.324616997912592;

b{1}(k75) = 0.745729161887539;
b{2}(k75) = 0.286739769811991;
b{3}(k75) = 5.483410575876488;


a{1}(k74) = -0.131154014410066;
a{2}(k74) = -10.970353971944270;
a{3}(k74) = 17.005199100223695;

b{1}(k74) = 0.207024146182337;
b{2}(k74) = 1.379590243136047;
b{3}(k74) = 24.034043337723823;


a{1}(k73) = -3.003406666981470;
a{2}(k73) = -40.051470493989385;
a{3}(k73) = 91.221185937417260;

b{1}(k73) = 0.922692987364494;
b{2}(k73) = 5.092973429201423;
b{3}(k73) = 5.657341198364528;


a{1}(k72) =  2.441808050468793;
a{2}(k72) = -84.603761901808824;
a{3}(k72) =  1.551777722664123e+02;

b{1}(k72) =  14.342014627725572;
b{2}(k72) =  35.566197719519877;
b{3}(k72) = -3.815066617815733e+02;


a{1}(k71) = -1.115011802967045e+03;
a{2}(k71) = -5.142353017288521e+03;
a{3}(k71) =  1.420327315943565e+04;

b{1}(k71) =   2.707949700044824;
b{2}(k71) =  17.133414667584038;
b{3}(k71) = -41.205549884638295;

end





