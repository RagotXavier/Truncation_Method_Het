var
r w Z K GDP u I C
rt wt Zt Kt GDPt ut It Ct  
a1 at1 c1 a2 at2 c2 a3 at3 c3 a4 at4 c4 a5 at5 c5 a6 at6 c6 a7 at7 c7 a8 at8 c8 a9 at9 c9 a10 at10 c10 
a11 at11 c11 a12 at12 c12 a13 at13 c13 a14 at14 c14 a15 at15 c15 a16 at16 c16 a17 at17 c17 a18 at18 c18 a19 at19 c19 a20 at20 c20 
a21 at21 c21 a22 at22 c22 a23 at23 c23 a24 at24 c24 a25 at25 c25 ; 

varexo eps; 

parameters
beta alpha abar delta gamma rho_u coef coef2;

alpha   = 0.36;
beta    = 0.98;
abar   = 0 ;
delta   = 0.025;
gamma   = 1;
rho_u   = 0.95;
coef  = 1.0;
coef2  = 1.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

c1 = 0.18100000*w + (1+r)*at1-a1 ; % equation 1
at1 = 10^-10 + 0.98070000*a1(-1)+ 0.01930000*a2(-1)+ 0.00000000*a3(-1)+ 0.00000000*a4(-1)+ 0.00000000*a5(-1); % equation 2
0.00000000*(c1)^-gamma = + beta*(1+r(+1))*(10^-10); % equation 3
c2 = 0.18100000*w + (1+r)*at2-a2 ; % equation 4
at2 = 10^-10 + 0.00480000*a6(-1)+ 0.98080000*a7(-1)+ 0.01440000*a8(-1)+ 0.00000000*a9(-1)+ 0.00000000*a10(-1); % equation 5
a2 =  steady_state(a2); % equation 6
c3 = 0.18100000*w + (1+r)*at3-a3 ; % equation 7
at3 = 10^-10 ; % equation 8
a3 =  steady_state(a3); % equation 9
c4 = 0.18100000*w + (1+r)*at4-a4 ; % equation 10
at4 = 10^-10 ; % equation 11
a4 =  steady_state(a4); % equation 12
c5 = 0.18100000*w + (1+r)*at5-a5 ; % equation 13
at5 = 10^-10 ; % equation 14
a5 =  steady_state(a5); % equation 15
c6 = 0.37400000*w + (1+r)*at6-a6 ; % equation 16
at6 = 10^-10 + 0.98070000*a1(-1)+ 0.01930000*a2(-1)+ 0.00000000*a3(-1)+ 0.00000000*a4(-1)+ 0.00000000*a5(-1); % equation 17
a6 =  steady_state(a6); % equation 18
c7 = 0.37400000*w + (1+r)*at7-a7 ; % equation 19
at7 = 10^-10 + 0.00480000*a6(-1)+ 0.98080000*a7(-1)+ 0.01440000*a8(-1)+ 0.00000000*a9(-1)+ 0.00000000*a10(-1); % equation 20
a7 =  steady_state(a7); % equation 21
c8 = 0.37400000*w + (1+r)*at8-a8 ; % equation 22
at8 = 10^-10 + 0.00000000*a11(-1)+ 0.00960000*a12(-1)+ 0.98080000*a13(-1)+ 0.00960000*a14(-1)+ 0.00000000*a15(-1); % equation 23
a8 =  steady_state(a8); % equation 24
c9 = 0.37400000*w + (1+r)*at9-a9 ; % equation 25
at9 = 10^-10 ; % equation 26
a9 =  steady_state(a9); % equation 27
c10 = 0.37400000*w + (1+r)*at10-a10 ; % equation 28
at10 = 10^-10 ; % equation 29
a10 =  steady_state(a10); % equation 30
c11 = 0.77300000*w + (1+r)*at11-a11 ; % equation 31
at11 = 10^-10 ; % equation 32
a11 =  steady_state(a11); % equation 33
c12 = 0.77300000*w + (1+r)*at12-a12 ; % equation 34
at12 = 10^-10 + 0.00480000*a6(-1)+ 0.98080000*a7(-1)+ 0.01440000*a8(-1)+ 0.00000000*a9(-1)+ 0.00000000*a10(-1); % equation 35
a12 =  steady_state(a12); % equation 36
c13 = 0.77300000*w + (1+r)*at13-a13 ; % equation 37
at13 = 10^-10 + 0.00000000*a11(-1)+ 0.00960000*a12(-1)+ 0.98080000*a13(-1)+ 0.00960000*a14(-1)+ 0.00000000*a15(-1); % equation 38
a13 =  steady_state(a13); % equation 39
c14 = 0.77300000*w + (1+r)*at14-a14 ; % equation 40
at14 = 10^-10 + 0.00000000*a16(-1)+ 0.00000000*a17(-1)+ 0.01440000*a18(-1)+ 0.98080000*a19(-1)+ 0.00480000*a20(-1); % equation 41
a14 =  steady_state(a14); % equation 42
c15 = 0.77300000*w + (1+r)*at15-a15 ; % equation 43
at15 = 10^-10 ; % equation 44
a15 =  steady_state(a15); % equation 45
c16 = 1.59760000*w + (1+r)*at16-a16 ; % equation 46
at16 = 10^-10 ; % equation 47
a16 =  steady_state(a16); % equation 48
c17 = 1.59760000*w + (1+r)*at17-a17 ; % equation 49
at17 = 10^-10 ; % equation 50
a17 =  steady_state(a17); % equation 51
c18 = 1.59760000*w + (1+r)*at18-a18 ; % equation 52
at18 = 10^-10 + 0.00000000*a11(-1)+ 0.00960000*a12(-1)+ 0.98080000*a13(-1)+ 0.00960000*a14(-1)+ 0.00000000*a15(-1); % equation 53
a18 =  steady_state(a18); % equation 54
c19 = 1.59760000*w + (1+r)*at19-a19 ; % equation 55
at19 = 10^-10 + 0.00000000*a16(-1)+ 0.00000000*a17(-1)+ 0.01440000*a18(-1)+ 0.98080000*a19(-1)+ 0.00480000*a20(-1); % equation 56
a19 =  steady_state(a19); % equation 57
c20 = 1.59760000*w + (1+r)*at20-a20 ; % equation 58
at20 = 10^-10 + 0.00000000*a21(-1)+ 0.00000000*a22(-1)+ 0.00000000*a23(-1)+ 0.01930000*a24(-1)+ 0.98070000*a25(-1); % equation 59
a20 =  steady_state(a20); % equation 60
c21 = 3.30230000*w + (1+r)*at21-a21 ; % equation 61
at21 = 10^-10 ; % equation 62
a21 =  steady_state(a21); % equation 63
c22 = 3.30230000*w + (1+r)*at22-a22 ; % equation 64
at22 = 10^-10 ; % equation 65
a22 =  steady_state(a22); % equation 66
c23 = 3.30230000*w + (1+r)*at23-a23 ; % equation 67
at23 = 10^-10 ; % equation 68
a23 =  steady_state(a23); % equation 69
c24 = 3.30230000*w + (1+r)*at24-a24 ; % equation 70
at24 = 10^-10 + 0.00000000*a16(-1)+ 0.00000000*a17(-1)+ 0.01440000*a18(-1)+ 0.98080000*a19(-1)+ 0.00480000*a20(-1); % equation 71
a24 =  steady_state(a24); % equation 72
c25 = 3.30230000*w + (1+r)*at25-a25 ; % equation 73
at25 = 10^-10 + 0.00000000*a21(-1)+ 0.00000000*a22(-1)+ 0.00000000*a23(-1)+ 0.01930000*a24(-1)+ 0.98070000*a25(-1); % equation 74
a25 =  steady_state(a25); % equation 75
w = (1 - alpha)*(Z^(1/(1-alpha)))*((r + delta)/alpha)^(alpha/(alpha-1)); % equation 76
(r + delta)/(alpha*Z) = (K/1)^(alpha-1); % equation 77
Z = 1 + u; % equation 78
u = rho_u*u(-1) + eps; % equation 79
GDP = Z*K(-1)^alpha*1^(1-alpha); % equation 80

C =+0.06101568*c1+0.00120078*c2+0.00000000*c3+0.00000000*c4+0.00000000*c5
+0.00120078*c6+0.24535891*c7+0.00360233*c8+0.00000000*c9+0.00000000*c10
+0.00000000*c11+0.00360233*c12+0.36803837*c13+0.00360233*c14+0.00000000*c15
+0.00000000*c16+0.00000000*c17+0.00360233*c18+0.24535891*c19+0.00120078*c20
+0.00000000*c21+0.00000000*c22+0.00000000*c23+0.00120078*c24+0.06101568*c25
; % equation 81
I =  K(+1) - ( 1 - delta)*K; % equation 82
K = +0.06101568*a1+0.00120078*a2+0.00000000*a3+0.00000000*a4+0.00000000*a5+0.00120078*a6+0.24535891*a7+0.00360233*a8+0.00000000*a9+0.00000000*a10+0.00000000*a11+0.00360233*a12+0.36803837*a13+0.00360233*a14+0.00000000*a15+0.00000000*a16+0.00000000*a17+0.00360233*a18+0.24535891*a19+0.00120078*a20+0.00000000*a21+0.00000000*a22+0.00000000*a23+0.00120078*a24+0.06101568*a25; % equation 83
ut = 100*u; % equation 84
Kt = 100*(K/29.36173135-1); % equation 85
Ct = 100*(C/steady_state(C)-1); % equation 86
wt = 100*(w/2.16063483-1); % equation 87
rt = 100*(r - 0.01639256); % equation 88
Zt = ut; % equation 89
GDPt = 100*(GDP/3.37599191-1); % equation 90
It = 100*(I/0.73404328-1); % equation 91

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;

a1 = 4.65100000;
at1 = 4.72894691;
c1 = 0.54654134;
a2 = 8.68970000;
at2 = 8.82470800;
c2 = 0.67074242;
a3 = 0.00000000;
at3 = 0.00000000;
c3 = 0.39107490;
a4 = 0.00000000;
at4 = 0.00000000;
c4 = 0.39107490;
a5 = 0.00000000;
at5 = 0.00000000;
c5 = 0.39107490;
a6 = 4.65120000;
at6 = 4.72894691;
c6 = 0.96334386;
a7 = 8.68990000;
at7 = 8.82470800;
c7 = 1.08754494;
a8 = 19.39780000;
at8 = 19.69258528;
c8 = 1.42567450;
a9 = 0.00000000;
at9 = 0.00000000;
c9 = 0.80807742;
a10 = 0.00000000;
at10 = 0.00000000;
c10 = 0.80807742;
a11 = 0.00000000;
at11 = 0.00000000;
c11 = 1.67017072;
a12 = 8.76130000;
at12 = 8.82470800;
c12 = 1.87823823;
a13 = 19.51360000;
at13 = 19.69258528;
c13 = 2.17196779;
a14 = 48.91020000;
at14 = 49.34345872;
c14 = 2.91229482;
a15 = 0.00000000;
at15 = 0.00000000;
c15 = 1.67017072;
a16 = 0.00000000;
at16 = 0.00000000;
c16 = 3.45183020;
a17 = 0.00000000;
at17 = 0.00000000;
c17 = 3.45183020;
a18 = 20.03540000;
at18 = 19.69258528;
c18 = 3.43182727;
a19 = 49.45600000;
at19 = 49.34345872;
c19 = 4.14815429;
a20 = 114.27170000;
at20 = 114.54430560;
c20 = 5.60210966;
a21 = 0.00000000;
at21 = 0.00000000;
c21 = 7.13506438;
a22 = 0.00000000;
at22 = 0.00000000;
c22 = 7.13506438;
a23 = 0.00000000;
at23 = 0.00000000;
c23 = 7.13506438;
a24 = 50.98710000;
at24 = 49.34345872;
c24 = 6.30028848;
a25 = 115.79510000;
at25 = 114.54430560;
c25 = 7.76194385;
r = 0.01639256;
w = 2.16063483;
K = 29.36173135;
Z = 1;
u = 0;
GDP = 3.37599191;
I = 0.73404328;
C = 2.64195059;
ut = 0;
Kt = 0;
Ct  = 0;
wt = 0;
rt = 0;
Zt = 0;
GDPt = 0;
It = 0;
end;

resid;

steady(nocheck);

check;

shocks;
 var eps; stderr 0.00312000;
 end;
options_.TeX=1;

stoch_simul (order=1,irf=200, periods=10000) ut Zt GDPt Ct Kt It rt wt u Z GDP C  K I r w;