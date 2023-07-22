%% CALLING THE MAIN VARIABLES

tic
clear
clc
warning('off','all');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

isOctave && warning('off', 'Octave:array-as-logical');


bool_IRF = true;
%bool_Trans = !bool_IRF;

load ../todynare_Truncation_N_2.mat

%~ alpha = alpha(1);
Nbin  = eco.Nbin(1);    % total number of bins
PI    = eco.Matab;   % transition matrix 
Sp    = eco.Sp; % size of each bin
K     = eco.A-0*eco.B; % capital
L     = 1; % laborC_eps
FKK   = alpha*(alpha-1)*K^(alpha-2)*L^(1-alpha);
FLK   = alpha*(1-alpha)*K^(alpha-1)*L^(-alpha);
Y     = K^alpha*L^(1-alpha); % GDP
I     = delta*K; % investment


%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u        = 0.95;
sigma_u      = 0.00312;

P            = ones(Nbin,1);
P(eco.indcc) = 0;

fid = fopen('DyTruncation.mod', 'w') ;
formatSpec = '%16.12f';

%% DEFINING THE ENDOGENOUS VARIABLES 
if bool_IRF
    str = ['@#define IRF = true\n']; fprintf(fid, str);
    str = ['@#define TRANSITION = false\n']; fprintf(fid, str);
else
    str = ['@#define IRF = false\n']; fprintf(fid, str);
    str = ['@#define TRANSITION = true\n']; fprintf(fid, str);
end

str = ['@#define EXTENDED_PATH_SIMULATION = false\n\n']; fprintf(fid, str);

str = ['var\n'] ; fprintf(fid, str);
str = ['R w FK FL FKK FLK u Z A K GDP I C TT \n'] ; fprintf(fid, str);  % variables
str = ['Rt wt ut Zt Kt GDPt It Ct Tt \n'] ; fprintf(fid, str);

isOctave && fflush(fid);

for h = 1:Nbin
    str = ['a',num2str(h),' '] ; fprintf(fid,str);     %  end-of-period
    str = ['at',num2str(h),' '] ; fprintf(fid, str);   % beginning-of-period
    str = ['c',num2str(h),' '] ; fprintf(fid, str);
    str = ['psih',num2str(h),' '] ; fprintf(fid, str);    
    str = ['lambda',num2str(h),' '] ; fprintf(fid, str);    
    str = ['lambdat',num2str(h),' '] ; fprintf(fid, str);    % last period 
    if mod(h,10)==0; % multiple of 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
isOctave && fflush(fid);

str = ['; \n\n'] ; fprintf(fid, str);

str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);
str = ['beta alpha theta abar delta gamma rho_u coef coef2;\n\n'] ; fprintf(fid, str);

str = ['alpha','   = ',num2str(alpha,formatSpec),';\n']; fprintf(fid, str);
str = ['beta','    = ',num2str(beta,formatSpec),';\n']; fprintf(fid, str);
str = ['theta','   = ',num2str(theta,formatSpec),';\n']; fprintf(fid, str);
str = ['abar','    = ',num2str(min(eco.aep),formatSpec),';\n']; fprintf(fid, str);
str = ['delta','   = ',num2str(delta,formatSpec),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(gamma,formatSpec),';\n']; fprintf(fid, str);
str = ['rho_u','   = ',num2str(rho_u),';\n']; fprintf(fid, str); 

str = ['coef','   = 1.0;\n']; fprintf(fid, str);
str = ['coef2','  = 1.;\n']; fprintf(fid, str);


%% EQUATIONS

equa = 1;
tol  = 10^-18;  % threshold for transition, to avoid to consider very small transitions

str  = ['%%%%%%%%%%%%%%%%%%%%%%%%']; fprintf(fid, str);
str  = ['\n\n'] ; fprintf(fid, str);

str  = ['model;\n\n']; fprintf(fid, str);
for h = 1:Nbin
    nsi = eco.ytype(h);

%%%%%%%%%%%%%%%%% BUDGET CONSTRAINT %%%%%%%%%%%%%%%
        str = ['c',num2str(h),' = ',num2str((eco.ytype(h)),formatSpec), '*w + R*','at',num2str(h),'-a', num2str(h),'- TT;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        
%%%%%%%%%%%%%%%%% DEFIITION OF A CONSTRAINT %%%%%%%%%%%%%%%
        str = ['at',num2str(h),' = 10^-18 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if (PI(h,hi))>tol   % hi to h
              str = ['+ ',num2str(Sp(hi)*PI(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        
%%%%%%%%%% EULER EQUATION %%%%%%%%%%%%%%%%%%%%%%%%
        if P(h)==1 % not constrained
            str = [num2str(eco.xsyn1(h),formatSpec),'*(c',num2str(h),')^(-gamma) =', num2str(eco.Res(h),formatSpec),' + beta*(R(+1))*(']; fprintf(fid, str);
            for hp = 1:Nbin
                    if (PI(hp,h)*eco.xsyn1(hp))>tol
                        str = ['+ ',num2str((PI(hp,h))*eco.xsyn1(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                    end
            end
            str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['a',num2str(h),' = ', num2str(eco.aep(h)/Sp(h),formatSpec),';\n']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end
        
%%%%%%%%%% DEFINITION OF PSI AND F.O.C PLANNER a %%%%%%%%%%%%%%%%%%%%%%%%
        
         
         if P(h)==1
            str = ['psih',num2str(h),'=',num2str(eco.xsyn1(h),formatSpec) ,'*(c',num2str(h),')^(-gamma) - (lambda',num2str(h),...
             '- R*lambdat',num2str(h),')*(',num2str(-gamma*eco.xsyn2(h),formatSpec),')*((c',num2str(h),')^(-gamma-1));'];fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
            str = ['psih',num2str(h),'  = beta*(R(+1))*(']; fprintf(fid, str);
            for hp = 1:Nbin
                if (PI(hp,h))>tol   % be careful transpose from h to hp
                    str = ['+ ',num2str((PI(hp,h)),formatSpec),'*psih',num2str(hp),'(+1)']; fprintf(fid, str);
                    % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
                end
            end 
            str = [') + ']; fprintf(fid, str);
            str = ['beta*(10^-18']; fprintf(fid, str);
            for hp = 1:Nbin
                %if (PI(hp,h))>tol   % be careful transpose from h to hp %%ERREUR: il faut sommer toutes les histoires
                    str = ['+ (psih',num2str(hp),'(+1))*',num2str(Sp(hp),formatSpec),'*((at',num2str(hp),'(+1))*FKK(+1) + FLK(+1)*' ,num2str((eco.ytype(hp)),formatSpec),')',]; fprintf(fid, str);
                    % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
                %end
            end 
            str = [') + ']; fprintf(fid, str);
            str = ['beta*FKK(+1)*(10^-18']; fprintf(fid, str);
            for hp = 1:Nbin
                %if (PI(hp,h))>tol   % be careful transpose from h to hp %%ERREUR: il faut sommer toutes les histoires
                    str = ['+',num2str(Sp(hp),formatSpec),'*(lambdat',num2str(hp),'(+1))*' , num2str(eco.xsyn1(hp),formatSpec) ,'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                    % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
                %end
             end 
            str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
                                      
        else
            str = ['psih',num2str(h),'=',num2str(eco.xsyn1(h),formatSpec) ,'*(c',num2str(h),')^(-gamma) + ',...
             'R*lambdat',num2str(h),'*(',num2str(-gamma*eco.xsyn2(h),formatSpec),')*((c',num2str(h),')^(-gamma-1));'];fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
            str = ['lambda',num2str(h),' = ',num2str(Pl.lambda(h),formatSpec),';\n']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
         end
        
%%%%%%%%%% DEFINITION OF LAMBDAT %%%%%%%%%%%%%%%%%%%%%%%%
        str = ['lambdat',num2str(h),' =']; fprintf(fid, str);
        for hi = 1:Nbin
            if PI(h,hi)>tol  % be careful transpose
              str = ['+ ',num2str((PI(h,hi))*Sp(hi)/Sp(h),formatSpec),'*lambda',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
        end
        str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        isOctave && fflush(fid);
end
          
%%%%%%%%%%%%%%%%% FOC PLANNER T %%%%%%%%%%%%%%%
 str = ['theta*(TT)^(theta-1) =     ']; fprintf(fid, str); 
    for hi = 1:Nbin
          if (Sp(hi))>tol
              str = [' +',num2str(Sp(hi),formatSpec) ,'*psih',num2str(hi) ]; fprintf(fid, str);
          end
         if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
   end             
 str = [';']; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
 
  
%% AGGREGATIONS
str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((R- 1+ delta)/alpha)^(alpha/(alpha-1));']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['(R-1 + delta)/(alpha*Z) = (K(-1)/1)^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['C +I + TT = GDP;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['FL = (1 - alpha)*Z*(K(-1)/1)^alpha;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['FK = alpha*Z*(K(-1)/1)^(alpha-1);'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['FLK = (1-alpha)*alpha*Z*(K(-1))^(alpha-1);'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['FKK = alpha*(alpha-1)*Z*(K(-1))^(alpha-2);'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['Z = 1 + u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_u*u(-1) + eps;'] ; fprintf(fid, str);    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%str = ['GDP = Z*K(-1)^alpha*1^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
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

%%%%%%%%%% TOTAL INVESTMENT %%%%%%%%%%%%%%%%%%%%%%%%
str = ['I =  K - ( 1 - delta)*K(-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%%%%%%%%% TOTAL SAVINGS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['A = '] ; fprintf(fid, str);
for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*a',num2str(h)] ; fprintf(fid, str);
    if mod(h,10)==10; tr = ['\n'] ;fprintf(fid, str);end;

end;
 str = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(C/',num2str(Ctot,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Rt = 100*(R - ',num2str(eco.R,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Tt = 100*(TT/',num2str(eco.TT(1),formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['\n'] ;         fprintf(fid, str);
str = ['end;\n\n'] ; fprintf(fid, str);
isOctave && fflush(fid);

%% STEADY STATE

%~ str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);

%str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
str = ['initval;\n\n'] ; fprintf(fid, str);

%fflush(fid);

for h = 1:Nbin
    %disp(h)
    str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['at',num2str(h),' = ',num2str(eco.abp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(h),' = ',num2str(eco.cp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['lambda',num2str(h),' = ',num2str(1*Pl.lambda(h),formatSpec),';\n'] ; fprintf(fid, str);    
    str = ['lambdat',num2str(h),' = ',num2str(1*Pl.lambdat(h),formatSpec),';\n'] ; fprintf(fid, str);    
    str = ['psih',num2str(h),' = ',num2str(Pl.psih(h),formatSpec),';\n'] ; fprintf(fid, str);     
end;

str = ['R = ', num2str(eco.R,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['FK = ', num2str((eco.R-1 + delta),formatSpec),';\n'] ; fprintf(fid, str);
str = ['FL = ', num2str(1*eco.w,formatSpec),';\n'] ; fprintf(fid, str);

str = ['FKK = ', num2str(FKK,formatSpec),';\n'] ; fprintf(fid, str);
str = ['FLK = ', num2str(FLK,formatSpec),';\n'] ; fprintf(fid, str);

str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['A = ', num2str(eco.A,formatSpec),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);

str = ['I = ', num2str(I,formatSpec),';\n'] ; fprintf(fid, str);
str = ['C = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['TT = ', num2str(eco.TT,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['Rt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['It = 0;\n'] ; fprintf(fid, str);
str = ['Tt = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SIMULATE IMPULSE RESPONSE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['@#if IRF \n\n']; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
%~ str = ['disp("here");\n\n'] ; fprintf(fid, str);

str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
%~ str = ['disp("there");\n\n'] ; fprintf(fid, str);

%str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
%~ str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
%~ str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

%~ str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=10000) Zt GDPt Ct Kt It Rt wt Tt Z GDP C K I R w TT;\n\n'] ; fprintf(fid, str);

str = ['@#endif\n\n']; fprintf(fid, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TRANSITION FROM A SHOCK IN PERIOD 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['@#if TRANSITION\n\n']; fprintf(fid, str);

str = ['endval;\n\n'] ; fprintf(fid, str);
for h = 1:Nbin   
    str = ['lambda',num2str(h),' = ',num2str(Pl.lambda(h),formatSpec),';\n'] ; fprintf(fid, str);         
end;
str = ['end;\n\n'] ; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
%str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
%str = ['check;\n\n'] ; fprintf(fid, str);

%str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
%~ str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['perfect_foresight_setup (periods=400);\n\n']; fprintf(fid, str);
str = ['perfect_foresight_solver(robust_lin_solve);\n\n']; fprintf(fid, str);

str = ['rplot TT ;\n'] ; fprintf(fid, str);
str = ['plot (C./GDP) ;\n'] ; fprintf(fid, str);
str = ['plot (TT./GDP) ;\n'] ; fprintf(fid, str);

str = ['@#endif\n\n']; fprintf(fid, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EXTEND PATH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['@#if EXTENDED_PATH_SIMULATION\n\n']; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
%str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
%str = ['check;\n\n'] ; fprintf(fid, str);

%str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

str = ['extended_path(periods=150, order=0);\n\n']; fprintf(fid, str);

str = ['rplot TT ;\n'] ; fprintf(fid, str);



str = ['@#endif\n\n']; fprintf(fid, str);

%fflush(fid);


dynare DyTruncation noclearall
if bool_IRF
    run post_dynare_stoch_simul_XR.m
else
    run post_dynare_perfect_foresight.m
end
toc
%~ save pathwelfare TT GDP C


%~ %%% RUNNING REGRESSION OF TT TO OBTAIN THE COEFFICIENTS

%~ %tbl = table (TT(1:200), ones(size(C(1:200))), C(1:200), GDP(1:200), Z(1:200),K(1:200) , ...
%~ %    'VariableNames',{'TT','Interc','C', 'GDP', 'Z', 'K'});

%~ %mdl = fitlm(tbl, 'TT ~ GDP + K ');
%~ X = [ones(size(C(1:200))) GDP(1:200) K(1:200)];
%~ coef= 1%regress(TT(1:200),X);

%~ TT_ss = oo_.steady_state(15);

%~ save tofigtruncation Zt_eps GDPt_eps Ct_eps Kt_eps It_eps Wt_eps Rt_eps wt_eps Tt_eps coef TT_ss






