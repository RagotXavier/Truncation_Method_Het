var
r w Z K GDP u I C
rt wt Zt Kt GDPt ut It Ct  
a1 at1 c1 a2 at2 c2 a3 at3 c3 a4 at4 c4 a5 at5 c5 a6 at6 c6 a7 at7 c7 a8 at8 c8 a9 at9 c9 a10 at10 c10 
a11 at11 c11 a12 at12 c12 a13 at13 c13 a14 at14 c14 a15 at15 c15 a16 at16 c16 a17 at17 c17 a18 at18 c18 a19 at19 c19 a20 at20 c20 
a21 at21 c21 a22 at22 c22 a23 at23 c23 a24 at24 c24 a25 at25 c25 ; 

varexo eps; 

parameters
beta alpha abar delta gamma rho_z;

alpha   = 0.36;
beta    = 0.98;
abar    = 0;
delta   = 0.025;
gamma   = 1;
rho_z   = 0.95;


model;

c1 = 0.181000000000*w + (1+r)*at1-a1; % equation 1
at1 = 10^-10 + 0.980700000000*a1(-1)+ 0.019300000000*a2(-1); % equation 2
a1 =  steady_state(a1); % equation 3
c2 = 0.181000000000*w + (1+r)*at2-a2; % equation 4
at2 = 10^-10 + 0.004800000000*a6(-1)+ 0.980800000000*a7(-1)+ 0.014400000000*a8(-1); % equation 5
1.144749942929*(c2^-gamma) = beta*(1+r(+1))*(+ 0.980700000000*(c1(+1)^-gamma)+ 0.019300000000*(c6(+1)^-gamma)); % equation 6
c3 = 0.181000000000*w + (1+r)*at3-a3; % equation 7
at3 = 10^-10 ; % equation 8
a3 =  steady_state(a3); % equation 9
c4 = 0.181000000000*w + (1+r)*at4-a4; % equation 10
at4 = 10^-10 ; % equation 11
a4 =  steady_state(a4); % equation 12
c5 = 0.181000000000*w + (1+r)*at5-a5; % equation 13
at5 = 10^-10 ; % equation 14
a5 =  steady_state(a5); % equation 15
c6 = 0.374000000000*w + (1+r)*at6-a6; % equation 16
at6 = 10^-10 + 0.980700000000*a1(-1)+ 0.019300000000*a2(-1); % equation 17
a6 =  steady_state(a6); % equation 18
c7 = 0.374000000000*w + (1+r)*at7-a7; % equation 19
at7 = 10^-10 + 0.004800000000*a6(-1)+ 0.980800000000*a7(-1)+ 0.014400000000*a8(-1); % equation 20
0.643312271893*(c7^-gamma) = beta*(1+r(+1))*(+ 0.005494799726*(c2(+1)^-gamma)+ 0.630960676273*(c7(+1)^-gamma)+ 0.009943940773*(c12(+1)^-gamma)); % equation 21
c8 = 0.374000000000*w + (1+r)*at8-a8; % equation 22
at8 = 10^-10 + 0.009600000000*a12(-1)+ 0.980800000000*a13(-1)+ 0.009600000000*a14(-1); % equation 23
0.863143214338*(c8^-gamma) = beta*(1+r(+1))*(+ 0.005494799726*(c2(+1)^-gamma)+ 0.630960676273*(c7(+1)^-gamma)+ 0.009943940773*(c12(+1)^-gamma)); % equation 24
c9 = 0.374000000000*w + (1+r)*at9-a9; % equation 25
at9 = 10^-10 ; % equation 26
a9 =  steady_state(a9); % equation 27
c10 = 0.374000000000*w + (1+r)*at10-a10; % equation 28
at10 = 10^-10 ; % equation 29
a10 =  steady_state(a10); % equation 30
c11 = 0.773000000000*w + (1+r)*at11-a11; % equation 31
at11 = 10^-10 ; % equation 32
a11 =  steady_state(a11); % equation 33
c12 = 0.773000000000*w + (1+r)*at12-a12; % equation 34
at12 = 10^-10 + 0.004800000000*a6(-1)+ 0.980800000000*a7(-1)+ 0.014400000000*a8(-1); % equation 35
0.690551442596*(c12^-gamma) = beta*(1+r(+1))*(+ 0.008286174858*(c8(+1)^-gamma)+ 0.794064307420*(c13(+1)^-gamma)+ 0.009359594633*(c18(+1)^-gamma)); % equation 36
c13 = 0.773000000000*w + (1+r)*at13-a13; % equation 37
at13 = 10^-10 + 0.009600000000*a12(-1)+ 0.980800000000*a13(-1)+ 0.009600000000*a14(-1); % equation 38
0.809608796309*(c13^-gamma) = beta*(1+r(+1))*(+ 0.008286174858*(c8(+1)^-gamma)+ 0.794064307420*(c13(+1)^-gamma)+ 0.009359594633*(c18(+1)^-gamma)); % equation 39
c14 = 0.773000000000*w + (1+r)*at14-a14; % equation 40
at14 = 10^-10 + 0.014400000000*a18(-1)+ 0.980800000000*a19(-1)+ 0.004800000000*a20(-1); % equation 41
1.088575212489*(c14^-gamma) = beta*(1+r(+1))*(+ 0.008286174858*(c8(+1)^-gamma)+ 0.794064307420*(c13(+1)^-gamma)+ 0.009359594633*(c18(+1)^-gamma)); % equation 42
c15 = 0.773000000000*w + (1+r)*at15-a15; % equation 43
at15 = 10^-10 ; % equation 44
a15 =  steady_state(a15); % equation 45
c16 = 1.597600000000*w + (1+r)*at16-a16; % equation 46
at16 = 10^-10 ; % equation 47
a16 =  steady_state(a16); % equation 48
c17 = 1.597600000000*w + (1+r)*at17-a17; % equation 49
at17 = 10^-10 ; % equation 50
a17 =  steady_state(a17); % equation 51
c18 = 1.597600000000*w + (1+r)*at18-a18; % equation 52
at18 = 10^-10 + 0.009600000000*a12(-1)+ 0.980800000000*a13(-1)+ 0.009600000000*a14(-1); % equation 53
0.974957774233*(c18^-gamma) = beta*(1+r(+1))*(+ 0.015675483060*(c14(+1)^-gamma)+ 1.158427774228*(c19(+1)^-gamma)+ 0.007212844108*(c24(+1)^-gamma)); % equation 54
c19 = 1.597600000000*w + (1+r)*at19-a19; % equation 55
at19 = 10^-10 + 0.014400000000*a18(-1)+ 0.980800000000*a19(-1)+ 0.004800000000*a20(-1); % equation 56
1.181104990037*(c19^-gamma) = beta*(1+r(+1))*(+ 0.015675483060*(c14(+1)^-gamma)+ 1.158427774228*(c19(+1)^-gamma)+ 0.007212844108*(c24(+1)^-gamma)); % equation 57
c20 = 1.597600000000*w + (1+r)*at20-a20; % equation 58
at20 = 10^-10 + 0.019300000000*a24(-1)+ 0.980700000000*a25(-1); % equation 59
1.599905549032*(c20^-gamma) = beta*(1+r(+1))*(+ 0.015675483060*(c14(+1)^-gamma)+ 1.158427774228*(c19(+1)^-gamma)+ 0.007212844108*(c24(+1)^-gamma)); % equation 60
c21 = 3.302300000000*w + (1+r)*at21-a21; % equation 61
at21 = 10^-10 ; % equation 62
a21 =  steady_state(a21); % equation 63
c22 = 3.302300000000*w + (1+r)*at22-a22; % equation 64
at22 = 10^-10 ; % equation 65
a22 =  steady_state(a22); % equation 66
c23 = 3.302300000000*w + (1+r)*at23-a23; % equation 67
at23 = 10^-10 ; % equation 68
a23 =  steady_state(a23); % equation 69
c24 = 3.302300000000*w + (1+r)*at24-a24; % equation 70
at24 = 10^-10 + 0.014400000000*a18(-1)+ 0.980800000000*a19(-1)+ 0.004800000000*a20(-1); % equation 71
1.502675855760*(c24^-gamma) = beta*(1+r(+1))*(+ 0.030878177096*(c20(+1)^-gamma)+ 1.820260607621*(c25(+1)^-gamma)); % equation 72
c25 = 3.302300000000*w + (1+r)*at25-a25; % equation 73
at25 = 10^-10 + 0.019300000000*a24(-1)+ 0.980700000000*a25(-1); % equation 74
1.856083009708*(c25^-gamma) = beta*(1+r(+1))*(+ 0.030878177096*(c20(+1)^-gamma)+ 1.820260607621*(c25(+1)^-gamma)); % equation 75
w = (1 - alpha)*(Z^(1/(1-alpha)))*((r + delta)/alpha)^(alpha/(alpha-1)); % equation 76
r  = (alpha*Z)*(K(-1)/1)^(alpha-1) - delta; % equation 77
Z = 1 + u; % equation 78
u = rho_z*u(-1) + eps; % equation 79
GDP = Z*K(-1)^alpha*1^(1-alpha); % equation 80

C =+0.061015683733*c1+0.001200777706*c2+0.000000000000*c3+0.000000000000*c4+0.000000000000*c5
+0.001200777706*c6+0.245358911212*c7+0.003602333117*c8+0.000000000000*c9+0.000000000000*c10
+0.000000000000*c11+0.003602333117*c12+0.368038366818*c13+0.003602333117*c14+0.000000000000*c15
+0.000000000000*c16+0.000000000000*c17+0.003602333117*c18+0.245358911212*c19+0.001200777706*c20
+0.000000000000*c21+0.000000000000*c22+0.000000000000*c23+0.001200777706*c24+0.061015683733*c25
; % equation 81
I =  K - (1 - delta)*K(-1); % equation 82
K = +0.061015683733*a1+0.001200777706*a2+0.000000000000*a3+0.000000000000*a4+0.000000000000*a5
+0.001200777706*a6+0.245358911212*a7+0.003602333117*a8+0.000000000000*a9+0.000000000000*a10
+0.000000000000*a11+0.003602333117*a12+0.368038366818*a13+0.003602333117*a14+0.000000000000*a15
+0.000000000000*a16+0.000000000000*a17+0.003602333117*a18+0.245358911212*a19+0.001200777706*a20
+0.000000000000*a21+0.000000000000*a22+0.000000000000*a23+0.001200777706*a24+0.061015683733*a25
; % equation 83
ut = 100*u; % equation 84
Kt = 100*(K/29.072362347015-1); % equation 85
Ct = 100*(C/steady_state(C)-1); % equation 86
wt = 100*(w/2.152944776511-1); % equation 87
rt = 100*(r - 0.016655763035); % equation 88
Zt = ut; % equation 89
GDPt = 100*(GDP/3.363976213299-1); % equation 90
It = 100*(I/0.726809058675-1); % equation 91

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initval;

a1  = 0.000000000000;
at1 = 0.167711210000;
c1  = 0.560187572721;
a2  = 8.689700000000;
at2 = 8.802382240000;
c2  = 0.648975637279;
a3  = 0.000000000000;
at3 = 0.000000000000;
c3  = 0.389683004549;
a4  = 0.000000000000;
at4 = 0.000000000000;
c4  = 0.389683004549;
a5  = 0.000000000000;
at5 = 0.000000000000;
c5  = 0.389683004549;
a6  = 0.000000000000;
at6 = 0.167711210000;
c6  = 0.975705914587;
a7  = 8.689900000000;
at7 = 8.802382240000;
c7  = 1.064293979146;
a8  = 19.397800000000;
at8 = 19.692585280000;
c8  = 1.427981660381;
a9  = 0.000000000000;
at9 = 0.000000000000;
c9  = 0.805201346415;
a10  = 0.000000000000;
at10 = 0.000000000000;
c10  = 0.805201346415;
a11  = 0.000000000000;
at11 = 0.000000000000;
c11  = 1.664226312243;
a12  = 8.761300000000;
at12 = 8.802382240000;
c12  = 1.851918944974;
a13  = 19.513600000000;
at13 = 19.692585280000;
c13  = 2.171206626209;
a14  = 48.91020