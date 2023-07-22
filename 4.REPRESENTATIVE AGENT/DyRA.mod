var
R w Z a K GDP u I c C W TT 
Rt wt Zt Kt GDPt ut It Ct Wt Tt;
varexo eps; 

parameters
beta alpha theta abar delta gamma rho_u tau coef coef2 coefz coefo coeft;


alpha   = 0.36;
beta    = 0.99001;
theta   = 0.23621;
tau   = 0.080085;
abar    = 6.2852e-11;
delta   = 0.025;
gamma   = 1.0001;
rho_u   = 0.95;
coefz     = -0.039898;
coefo     = 0.068017;
coeft     = 0.0020753;
coef   = 1.0;
coef2  = 1.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

c  = w + (R)*a(-1) - a - TT;
(c)^-gamma = beta*(R(+1))*((c(+1))^-gamma );
a = K;
w = (1 - alpha)*(Z^(1/(1-alpha)))*((R-1 + delta)/alpha)^(alpha/(alpha-1));
(R-1 + delta)/(alpha*Z) = (K/1)^(alpha-1);
u = rho_u*u(-1) + eps;
Z = 1 + u;
GDP = Z*K(-1)^alpha*1^(1-alpha);
I =  K(+1) - ( 1 - delta)*K;
C =  GDP - I - TT;
TT =  coefz + coefo*GDP + coeft*K;
W  =(((c^(1- gamma) -1)/(1-gamma) + 10)+ (TT)^(theta) );
ut = 100*u;
Kt = 100*(K/40.59010827-1);
Ct = 100*(c/2.47487278-1);
wt = 100*(w/2.42779073-1);
Rt = 100*(R - 1.00864446);
Zt = ut;
GDPt = 100*(GDP/3.79342302-1);
It = 100*(I/1.01475271-1);
Wt = 100*(exp(( W-11.46959010)));Tt = 100*(TT/0.30379753-1);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initval;

a = 40.59010827;
c = 2.47487278;
R = 1.00864446;
w = 2.42779073;
K = 40.59010827;
Z = 1;
u = 0;
GDP = 3.79342302;
I = 1.01475271;
C = 2.47540948;
TT = 0.30233731;
W = 11.46959010;
ut = 0;
Kt = 0;
Ct  = 0;
wt = 0;
Rt = 0;
Zt = 0;
GDPt = 0;
It = 0;
Wt = 0;
Tt = 0;
end;

resid;

options_.solve_tolf=10^-6;

steady(maxit = 100,solve_algo=3);

shocks;
 var eps; stderr 0.00312000;
 end;
options_.TeX=1;

stoch_simul (order=1,irf=200, periods=1000) Zt GDPt Ct Kt It Wt Rt wt Tt Z GDP C K I W R w TT;

