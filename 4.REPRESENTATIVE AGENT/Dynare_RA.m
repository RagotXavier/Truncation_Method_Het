
%% CALLING THE MAIN VARIABLES

clear
clc
warning('off','all');
load todynare_Truncation.mat
load tofigtruncation.mat coef TT_ss


%% CALCULATE THE STEADY STATE VALUES

R = eco.Agg.R;
w=eco.Agg.w;
L = 1;
K = L*(alpha/(R-1+delta))^(1/(1-alpha));
Y = K^alpha*L^(1-alpha);
C = K^alpha*L^(1-alpha)+(-delta)*K-G;
I = delta*K;
TT =eco.Agg.TT;

Uc = (C)^-gamma;
Ucc = -gamma*(C)^(-gamma-1);

%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u        = 0.95;
sigma_u     = 0.00312;

fid = fopen('DyRA.mod', 'w') ;
formatSpec = '%16.8f';

%% DEFINING THE ENDOGENOUS VARIABLES 

str = ['var\n'] ; fprintf(fid, str);
str = ['R w Z a K GDP u I c C W TT \n'] ; fprintf(fid, str);  % variables
str = ['Rt wt Zt Kt GDPt ut It Ct Wt Tt;\n'] ; fprintf(fid, str);

str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);
str = ['beta alpha theta abar delta gamma rho_u tau coef coef2 coefz coefo coeft;\n\n'] ; fprintf(fid, str);

str = ['\nalpha','   = ',num2str(alpha),';\n']; fprintf(fid, str);
str = ['beta','    = ',num2str(beta),';\n']; fprintf(fid, str);
str = ['theta','   = ',num2str(theta),';\n']; fprintf(fid, str);
str = ['tau','   = ',num2str(tau),';\n']; fprintf(fid, str);
str = ['abar','    = ',num2str(min(eco.aep)),';\n']; fprintf(fid, str);
str = ['delta','   = ',num2str(delta),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(gamma),';\n']; fprintf(fid, str);
str = ['rho_u','   = ',num2str(rho_u),';\n']; fprintf(fid, str); 

str         = ['coefz  ','   = ',num2str(coef(1)),';\n']; fprintf(fid, str);
str         = ['coefo  ','   = ',num2str(coef(2)),';\n']; fprintf(fid, str);
str         = ['coeft  ','   = ',num2str(coef(3)),';\n']; fprintf(fid, str);

str = ['coef','   = 1.0;\n']; fprintf(fid, str);
str = ['coef2','  = 1.;\n']; fprintf(fid, str);

%% EQUATIONS

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);
str = ['model;\n\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%% BUDGET CONSTRAINT %%%%%%%%%%%%%%%
str = ['c  = w + (R)*a(-1) - a - TT;\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%% EULER EQUATION %%%%%%%%%%%%%%%
str = ['(c)^-gamma = beta*(R(+1))*((c(+1))^-gamma );\n';]; fprintf(fid, str);

%% AGGREGATIONS
str = ['a = K;\n'] ; fprintf(fid, str); 
str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((R-1 + delta)/alpha)^(alpha/(alpha-1));\n']; fprintf(fid, str); 
str = ['(R-1 + delta)/(alpha*Z) = (K/1)^(alpha-1);\n']; fprintf(fid, str);
str = ['u = rho_u*u(-1) + eps;\n'] ; fprintf(fid, str);   
str = ['Z = 1 + u;\n'] ; fprintf(fid, str);  
str = ['GDP = Z*K(-1)^alpha*1^(1-alpha);\n'] ; fprintf(fid, str);   
str = ['I =  K(+1) - ( 1 - delta)*K;\n']; fprintf(fid, str); 
str = ['C =  GDP - I - TT;\n']; fprintf(fid, str);
str = ['TT =  coefz + coefo*GDP + coeft*K;\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%% WELFARE %%%%%%%%%%%%%%%
str = ['W  =(((c^(1- gamma) -1)/(1-gamma) + 10)+ (TT)^(theta) );\n']; fprintf(fid, str);

%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;\n'] ; fprintf(fid, str);  
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);\n'] ; fprintf(fid, str);
str = ['Ct = 100*(c/',num2str(C,formatSpec),'-1);\n'] ; fprintf(fid, str);
str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);\n'] ; fprintf(fid, str);  
str = ['Rt = 100*(R - ',num2str(eco.R,formatSpec),');\n'] ; fprintf(fid, str); 
str = ['Zt = ut;\n'] ; fprintf(fid, str);
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);\n'] ; fprintf(fid, str);  
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);\n'] ; fprintf(fid, str);  
%str = ['Wt = 100*(exp((1-beta)*( W-',num2str(Ws,formatSpec),')) -1);'] ; fprintf(fid, str);  
str = ['Wt = 100*(exp(( W-',num2str(Ws,formatSpec),')));'] ; fprintf(fid, str);  
str = ['Tt = 100*(TT/',num2str(eco.Agg.TT,formatSpec),'-1);\n'] ; fprintf(fid, str);  

str = ['end;\n\n'] ; fprintf(fid, str);


%% STEADY STATE

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

%str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
str = ['initval;\n\n'] ; fprintf(fid, str);

str = ['a = ',num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['c = ',num2str(C,formatSpec),';\n'] ; fprintf(fid, str);
str = ['R = ', num2str(eco.R,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);
str = ['I = ', num2str(I,formatSpec),';\n'] ; fprintf(fid, str);
str = ['C = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['TT = ', num2str(TT_ss,formatSpec),';\n'] ; fprintf(fid, str);
str = ['W = ', num2str(Ws,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['Rt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['It = 0;\n'] ; fprintf(fid, str);
str = ['Wt = 0;\n'] ; fprintf(fid, str);
str = ['Tt = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
%str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
%str = ['check;\n\n'] ; fprintf(fid, str);

str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
%str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=1000) Zt GDPt Ct Kt It Wt Rt wt Tt Z GDP C K I W R w TT;\n\n'] ; fprintf(fid, str);

%% CALLING DYNARE

dynare DyRA noclearall

%% GENERATING THE RESULTS

save tofigRA Zt_eps GDPt_eps Ct_eps Kt_eps It_eps Wt_eps Rt_eps wt_eps Tt_eps

