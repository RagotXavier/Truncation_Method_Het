
%% CALLING THE MAIN VARIABLES

clear
clc
warning('off','all');
load ../todynare_Truncation.mat



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

fid = fopen('DyTruncation.mod', 'w') ;
formatSpec = '%16.8f';

%% DEFINING THE ENDOGENOUS VARIABLES 

str = ['var\n'] ; fprintf(fid, str);
str = ['r w Z A K GDP u I C Cb Cm Ctp W TT \n'] ; fprintf(fid, str);  % variables
str = ['rt wt Zt Kt GDPt ut It Ct Cbt Cmt Ctt Wt \n'] ; fprintf(fid, str);

for p = 1:Nbin
    str = ['a',num2str(p),' '] ; fprintf(fid,str);     %  end-of-period
    str = ['at',num2str(p),' '] ; fprintf(fid, str);   % beginning-of-period
    str = ['c',num2str(p),' '] ; fprintf(fid, str);
    if mod(p,10)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;

end;

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

%% EQUATIONS

equa = 1;
tol = 10^-8;  % threshold for transition, to avoid to consider very small transitions

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

str = ['model;\n\n']; fprintf(fid, str);
for h = 1:Nbin
    nsi = eco.ytype(h);

%%%%%%%%%%%%%%%%% BUDGET CONSTRAINT %%%%%%%%%%%%%%%
        str = ['c',num2str(h),' = ',num2str((eco.states(nsi)),formatSpec), '*w + (1+r)*','at',num2str(h),'-a', num2str(h),'+ TT;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        
%%%%%%%%%%%%%%%%% DEFIITION OF A CONSTRAINT %%%%%%%%%%%%%%%
        str = ['at',num2str(h),' = 10^-10 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if (PI(h,hi))>tol   % hi to h
              str = ['+ ',num2str(Sp(hi)*PI(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        
%%%%%%%%%% EULER EQUATION %%%%%%%%%%%%%%%%%%%%%%%%
        if P(h)==1 % not constrained
            str = [num2str(eco.xsip(h),formatSpec),'*(c',num2str(h),' - ','0)^-gamma =', num2str(eco.Res(h)/Sp(h),formatSpec),' + beta*(1+r(+1))*(10^-10']; fprintf(fid, str);
            for hp = 1:Nbin
                    if (PI(hp,h)*eco.xsip(hp))>tol
                        str = ['+ ',num2str((PI(hp,h))*eco.xsip(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                    end
            end
            str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['a',num2str(h),' =  steady_state(a',num2str(h),');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end
end
  
%% AGGREGATIONS
str = ['GDP -delta*K(-1) = TT +r*A(-1) + w;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((r + delta)/alpha)^(alpha/(alpha-1));']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['(r + delta)/(alpha*Z) = (K/1)^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Z = 1 + u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_u*u(-1) + eps;'] ; fprintf(fid, str);    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDP = Z*K(-1)^alpha*1^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['A = K;'] ; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%% TOTAL CONSUMPTION %%%%%%%%%%%%%%%%%%%%%%%%
str = ['\n'] ;         fprintf(fid, str);
str = ['C ='] ; fprintf(fid, str);
for h = 1:Nbin
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%% TOTAL WELFARE %%%%%%%%%%%%%%%%%%%%%%%%
str = ['\n'] ;         fprintf(fid, str);
str = ['W ='] ; fprintf(fid, str);
for h = 1:Nbin
        %str = ['+',num2str(Sp(h),formatSpec),'*', num2str(eco.xsyn(h),formatSpec) ,'*(((c',num2str(h),'/', num2str(Sp(h),formatSpec) ')^(1- gamma) -1)/(1-gamma) + 10)'] ; fprintf(fid, str);
        str = ['+',num2str(Sp(h),formatSpec),'*', num2str(eco.xsyn(h),formatSpec) ,'*(((c',num2str(h),')^(1- gamma) -1)/(1-gamma) + 10)'] ; fprintf(fid, str);
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['I =  K(+1) - ( 1 - delta)*K;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%%%%%%%%% CONSUMPTION BOTTOM AGENTS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['\n'] ;         fprintf(fid, str);
str = ['#Cba =('] ; fprintf(fid, str);
for h = [b']
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['#ba = ']; fprintf(fid, str);
for h = [b']
     str    = ['+',num2str(Sp(h),formatSpec),'*1']; fprintf(fid, str);
      if mod(h,5)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  


str = ['Cb = Cba/ba;']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%%%%%%%%% CONSUMPTION MEDIAN AGENTS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['\n'] ;         fprintf(fid, str);
str = ['#Cma =('] ; fprintf(fid, str);
for h = [m']
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['#ma = ']; fprintf(fid, str);
for h = [m']
     str    = ['+',num2str(Sp(h),formatSpec),'*1']; fprintf(fid, str);
      if mod(h,5)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  


str = ['Cm = Cma/ma;']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%%%%%%%%% CONSUMPTION TOP AGENTS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['\n'] ;         fprintf(fid, str);
str = ['#Ctpa =('] ; fprintf(fid, str);
for h = [t']
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['#ta = ']; fprintf(fid, str);
for h = [t']
     str    = ['+',num2str(Sp(h),formatSpec),'*1']; fprintf(fid, str);
      if mod(h,5)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  


str = ['Ctp = Ctpa/ta;']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   


%%%%%%%%%% TOTAL SAVINGS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['A = '] ; fprintf(fid, str);
for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*a',num2str(h)] ; fprintf(fid, str);
    if mod(p,10)==10; tr = ['\n'] ;fprintf(fid, str);end;

end;
 str = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(C/steady_state(C)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Cbt = 100*(Cb/',num2str(Cbs,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Cmt = 100*(Cm/',num2str(Cms,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ctt = 100*(Ctp/',num2str(Cts,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = 100*(r - ',num2str(eco.R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Wt = 100*(W-',num2str(Ws,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['\n'] ;         fprintf(fid, str);
str = ['end;\n\n'] ; fprintf(fid, str);

%% STEADY STATE

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
%str = ['initval;\n\n'] ; fprintf(fid, str);
for h = 1:Nbin
    str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['at',num2str(h),' = ',num2str(eco.abp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(h),' = ',num2str(eco.cp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
 end;
str = ['r = ', num2str(eco.R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['A = ', num2str(eco.A,formatSpec),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);

str = ['I = ', num2str(I,formatSpec),';\n'] ; fprintf(fid, str);
str = ['C = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Cb = ', num2str(Cbs,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Cm = ', num2str(Cms,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Ctp = ', num2str(Cts,formatSpec),';\n'] ; fprintf(fid, str);
str = ['TT = ', num2str(eco.Agg.TT,formatSpec),';\n'] ; fprintf(fid, str);
str = ['W = ', num2str(Ws,formatSpec),';\n'] ; fprintf(fid, str);

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

str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
str = ['check;\n\n'] ; fprintf(fid, str);
str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);


str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=10000) ut Zt GDPt Ct Cbt Cmt Ctt Kt It Wt rt wt u Z GDP C Cb Cm Ctp K I W r w;'] ; fprintf(fid, str);


%% CALLING DYNARE

dynare DyTruncation noclearall



%% GENERATING THE RESULTS

save tofigtruncation ut_eps Zt_eps GDPt_eps Ct_eps Cbt_eps Cmt_eps Ctt_eps Kt_eps It_eps Wt_eps rt_eps wt_eps 
%save totaltruncation oo_

save tofigtruncationsim u Z GDP C Cb Cm Ctp K I W r w 
save tofigtruncationsimt ut Zt GDPt Ct Cbt Cmt Ctt Kt It Wt rt wt 


