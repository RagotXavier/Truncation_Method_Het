var
r w Z K GDP u I C  W 
rt wt Zt Kt GDPt ut It Ct Wt 
; 

varexo eps; 

parameters
beta alpha abar delta gamma rho_u coef coef2;

alpha   = 0.36;
beta    = 0.98;
abar   = 4.0945e-05;
delta   = 0.025;
gamma   = 1.0001;
rho_u   = 0.95;
coef  = 1.0;
coef2  = 1.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

C + K = (1+r)*K(-1) +w; % equation 1
C^-gamma  = beta*(1+r(+1))*C(+1)^-gamma; % equation 2
w = (1 - alpha)*Z*K^alpha; % equation 3
r + delta = (alpha*Z)*(K(-1))^(alpha-1); % equation 4
Z = 1 + u; % equation 5
u = rho_u*u(-1) + eps; % equation 6
GDP = Z*K(-1)^alpha*1^(1-alpha); % equation 7
I =  K(+1) - ( 1 - delta)*K; % equation 8

W = (C^(1- gamma) -1)/(1-gamma) + 10;ut = 100*u; % equation 9
Kt = 100*(K/25.40685351-1); % equation 10
Ct = 100*(C/steady_state(C)-1); % equation 11
wt = 100*(w/2.05098409-1); % equation 12
rt = 100*(r - 0.02040816); % equation 13
Zt = ut; % equation 14
GDPt = 100*(GDP/3.20466265-1); % equation 15
It = 100*(I/0.63517134-1); % equation 16
Wt = 100*(W-10.94366342); % equation 17

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;

r = 1/beta-1;
w = 2.05098409;
K = 25.40685351;
Z = 1;
u = 0;
GDP = 3.20466265;
I = 0.63517134;
C = 2.56949131;
W = 10.94366342;
ut = 0;
Kt = 0;
Ct  = 0;
Cbt  = 0;
Cmt  = 0;
Ctt  = 0;
wt = 0;
rt = 0;
Zt = 0;
GDPt = 0;
It = 0;
Wt = 0;
end;

resid;

steady(nocheck);

check;

shocks;
 var eps; stderr 0.00312000;
 end;
options_.TeX=1;

stoch_simul (order=1,irf=200, periods=10000) ut Zt GDPt Ct  Kt It Wt rt wt u Z GDP C K I W r w;