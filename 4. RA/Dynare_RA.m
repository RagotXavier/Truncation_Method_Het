
%% CALLING THE MAIN VARIABLES

clear
clc
warning('off','all');
load '../1. Truncation/todynare_Truncation.mat'



Nbin  = eco.Nbin;    % total number of bins
PI    = eco.Matab;   % transition matrix 
Sp    = eco.Sp; % size of each bin
K     = eco.A-eco.B; % capital
L  = 1; % labor
Y     = K^alpha*L^(1-alpha); % GDP
I = delta*K; % investment

%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u        = 0.95;
sigma_u     = 0.00312;

P            = ones(Nbin,1);
P(eco.indcc) = 0;

fid = fopen('DyRa.mod', 'w') ;
formatSpec = '%16.8f';

%% DEFINING THE ENDOGENOUS VARIABLES 

str = ['var\n'] ; fprintf(fid, str);
str = ['r w Z K GDP u I C  W \n'] ; fprintf(fid, str);  % variables
str = ['rt wt Zt Kt GDPt ut It Ct Wt \n'] ; fprintf(fid, str);

str = ['; \n\n'] ; fprintf(fid, str);

str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);
str = ['beta alpha abar delta gamma rho_u coef coef2;\n\n'] ; fprintf(fid, str);

str = ['alpha','   = ',num2str(alpha),';\n']; fprintf(fid, str);
str = ['beta','    = ',num2str(beta),';\n']; fprintf(fid, str);
str = ['abar','   = ',num2str(min(eco.aep)),';\n']; fprintf(fid, str);
str = ['delta','   = ',num2str(delta),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(gamma),';\n']; fprintf(fid, str);
str = ['rho_u','   = ',num2str(rho_u),';\n']; fprintf(fid, str); 

str = ['coef','  = 1.0;\n']; fprintf(fid, str);
str = ['coef2','  = 1.;\n']; fprintf(fid, str);


%%% computing the steady state in the comple-market economy
Kss = (alpha*beta/(1-beta+beta*delta))^(1/(1-alpha));
Css = Kss^alpha - delta*Kss;
rss = 1/beta-1;
wss = (1-alpha)*Kss^alpha;
Iss = delta*Kss;
GDPss = Kss^alpha;
Wss = (Css^(1- gamma) -1)/(1-gamma) + 10;

%% EQUATIONS

equa = 1;
tol = 10^-8;  % threshold for transition, to avoid to consider very small transitions

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

str = ['model;\n\n']; fprintf(fid, str);
%%%%%%%%%%%%%%%%% economy %%%%%%%%%%%%%%%
        str = ['C + K = (1+r)*K(-1) +w;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['C^-gamma  = beta*(1+r(+1))*C(+1)^-gamma;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['w = (1 - alpha)*Z*K^alpha;']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
        str = ['r + delta = (alpha*Z)*(K(-1))^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['Z = 1 + u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['u = rho_u*u(-1) + eps;'] ; fprintf(fid, str);    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['GDP = Z*K(-1)^alpha*1^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        str = ['I =  K(+1) - ( 1 - delta)*K;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
%%%%%%%%%% TOTAL WELFARE %%%%%%%%%%%%%%%%%%%%%%%%
    str = ['\n'] ;         fprintf(fid, str);
    str = ['W = (C^(1- gamma) -1)/(1-gamma) + 10;'] ; fprintf(fid, str);
%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
    str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Kt = 100*(K/',num2str(Kss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Ct = 100*(C/steady_state(C)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['wt = 100*(w/',num2str(wss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['rt = 100*(r - ',num2str(1/beta-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Zt = ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['GDPt = 100*(GDP/',num2str(GDPss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['It = 100*(I/',num2str(Iss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Wt = 100*(W-',num2str(Wss,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    str = ['\n'] ;         fprintf(fid, str);
    str = ['end;\n\n'] ; fprintf(fid, str);

%% STEADY STATE

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
str = ['r = 1/beta-1;\n'] ; fprintf(fid, str);
str = ['w = ', num2str(wss,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(Kss,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(GDPss,formatSpec),';\n'] ; fprintf(fid, str);

str = ['I = ', num2str(Iss,formatSpec),';\n'] ; fprintf(fid, str);
str = ['C = ', num2str(Css,formatSpec),';\n'] ; fprintf(fid, str);
str = ['W = ', num2str(Wss,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['Cbt  = 0;\n'] ; fprintf(fid, str);
str = ['Cmt  = 0;\n'] ; fprintf(fid, str);
str = ['Ctt  = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['rt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['It = 0;\n'] ; fprintf(fid, str);
str = ['Wt = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
%str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
%str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
%str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
str = ['check;\n\n'] ; fprintf(fid, str);
str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

%str = ['extended_path (periods=100, order=0);'] ; fprintf(fid, str);   
%str = ['stoch_simul (order=1, hp_filter =1600)  u Z GDP C Cb Cm Ctp K I W r w;'] ; fprintf(fid, str);    
%str = ['stoch_simul (order=1,irf=100, hp_filter =1600)  u Z GDP C I w;'] ; fprintf(fid, str); 
%str = ['stoch_simul (order=1,irf=200, periods = 10000)  u Z GDP C K I r w;'] ; fprintf(fid, str);
str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=10000) ut Zt GDPt Ct  Kt It Wt rt wt u Z GDP C K I W r w;'] ; fprintf(fid, str);

%% CALLING DYNARE

dynare DyRa noclearall

%% GENERATING THE RESULTS

save tofigRa ut_eps Zt_eps GDPt_eps Ct_eps Kt_eps It_eps Wt_eps rt_eps wt_eps 
%save totalRa oo_

save tofigRasim u Z GDP C K I W r w 
save tofigRasimt ut Zt GDPt Ct Kt It Wt rt wt 


