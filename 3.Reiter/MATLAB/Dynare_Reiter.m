
%% CALLING THE MAIN VARIABLES

clear
clc
warning('off','all');
warning('off', 'Octave:array-as-logical');

disp('Writing the mod file...')

load ../todynare_Reiter_FLG.mat
load ../rule_TT.mat
Ntot = na.*ns;


Ny = ns; % idiosyncratic shocks
thres       = 0; % threshold to consider 0 (not to have too many equations)
Y = K.^alpha*Ltot.^(1-alpha); % GDP
I = delta*K; % investment

%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u       = 0.95;
sigma_u     = 0.00312;

fid         = fopen('DyReiter.mod', 'w') ;
formatSpec  = '%16.12f';


%% DEFINING THE ENDOGENOUS VARIABLES

str = ['var\n'] ; fprintf(fid, str);
str = ['R w Z K GDP u Pop I C  W TT\n'] ; fprintf(fid, str);  % variables
str = ['Rt wt Zt Kt GDPt ut It Ct Wt Tt\n'] ; fprintf(fid, str);

for i = 1:Ntot
    str     = ['ap',num2str(i),' '] ; fprintf(fid,str);
    str     = ['c',num2str(i),' '] ; fprintf(fid, str);
    str     = ['S',num2str(i),' '] ; fprintf(fid, str);
    str     = ['ww',num2str(i),' '] ; fprintf(fid, str);
    if mod(i,10)==0; % multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;

for i = 1:Ny
    str     = ['l',num2str(i),' '] ; fprintf(fid, str);
    if mod(i,10)==0; % multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;
str         = ['; \n\n'] ; fprintf(fid, str);fflush(fid);

str         = ['varexo eps; \n\n'] ; fprintf(fid, str);
str         = ['parameters\n'] ; fprintf(fid, str);
str         = ['beta alpha theta abar delta gamma rho_u tau a0 aC aGDP TT_ss\n'] ; fprintf(fid, str);
str = ['; \n\n'] ; fprintf(fid, str);

l0          = Ltot;
abar        = abar;

str         = ['alpha','   = ',num2str(alpha,formatSpec),';\n']; fprintf(fid, str);
str         = ['beta ','   = ',num2str(beta,formatSpec),';\n']; fprintf(fid, str);
str         = ['theta ','   = ',num2str(theta,formatSpec),';\n']; fprintf(fid, str);
str         = ['abar ','   = ',num2str(abar,formatSpec),';\n']; fprintf(fid, str);
str         = ['gamma','   = ',num2str(gamma,formatSpec),';\n']; fprintf(fid, str);
str         = ['rho_u','   = ',num2str(rho_u,formatSpec),';\n']; fprintf(fid, str);
str         = ['delta','   = ',num2str(delta,formatSpec),';\n']; fprintf(fid, str);
str         = ['l    ','   = ',num2str(l0,formatSpec),';\n']; fprintf(fid, str);
str         = ['tau  ','   = ',num2str(tau,formatSpec),';\n']; fprintf(fid, str);
str         = ['TT_ss  ','   = ',num2str(TT_ss,formatSpec),';\n']; fprintf(fid, str);
str         = ['a0  ','   = ',num2str(a0,formatSpec),';\n']; fprintf(fid, str);
str         = ['aC  ','   = ',num2str(aC,formatSpec),';\n']; fprintf(fid, str);
str         = ['aGDP  ','   = ',num2str(aGDP,formatSpec),';\n']; fprintf(fid, str);

str = ['\n\n'] ; fprintf(fid, str);
yv = ytype; % vector of productivity type from 1 to ns*na (for budget constraint)
yvv = [1:1:ns]; % vector of productivity type from 1 to ns
Th = 10^-8;
indc    = find(apol <= Th); % credit constraind agents
indnc = find(apol >Th);

%% EQUATIONS

%~ str     = ['%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str     = ['\n\n'] ; fprintf(fid, str);

str     = ['model;\n\n']; fprintf(fid, str);
equa    = 1;

%%%%%%%%%%%%%%%%% LABOR SUPPLY %%%%%%%%%%%%%%%%%%%%
for i = 1:Ny
    str = ['l',num2str(i),' =  (',num2str(l0,formatSpec),')^',num2str(1),';']; fprintf(fid, str);
    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end;

%%%%%%%%%%%%%%%%% BUDGET CONSTRAINT %%%%%%%%%%%%%%%
for i = 1:Ntot
    str = ['c',num2str(i),'+ap',num2str(i), ' = w*',num2str(states(yv(i)),formatSpec),'*l',num2str(yv(i)),'+ (R)*',num2str(aGrid(i),formatSpec) ,'- TT;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end;
str = ['\n\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%  WEALTH TRANSITION  %%%%%%%%%%%%%%
for i = 1:Ntot
     str = ['ap',num2str(i),' =  ww',num2str(i),'*',num2str(aGrid(Vind(i)),formatSpec),'+ (1 -ww',num2str(i),')*',num2str(aGrid(Vind(i)+1),formatSpec),';']; fprintf(fid, str);
     str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

end;
str = ['\n\n']; fprintf(fid, str);

%%%%%%%%%% EULER EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str     = ['%%%%%%%% equation ',num2str(equa),': Euler\n\n']; fprintf(fid, str);

for i = 1:length(indc)
     str = ['ap',num2str(indc(i)),' = ',num2str(abar) ,';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end;

for k = 1:length(indnc)
    i = indnc(k);
    yi = ytype(i); % current productivity status
    threse = 0;
          str = ['(c',num2str(i),')^-gamma  = beta*(R(+1))*('];fprintf(fid, str);
          for j = 1:ns % next period productivity status
                  if (Trans(yi,j)>threse) % does agent i reach agent j ? effective transition toward the bin j
                     str = ['+',num2str(Trans(yi,j),formatSpec),'*( ww',num2str(i),'*c',num2str(Vind(i)+(j-1)*na),'(+1)+(1-ww',num2str(i),')*c',num2str(Vind(i)+1+(j-1)*na),'(+1) )^-gamma']; fprintf(fid, str);
                  end;
          end;
            str = [')+',num2str(resE(i),formatSpec),';']; fprintf(fid, str);
            str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%  SHARES  %%%%%%%%%%%%%%%%
 str    = [' %%%%%%%% equation ',num2str(equa),': Shares \n\n']; fprintf(fid, str);
for i = 1:na
    for j=1:ns
        pass=0;
        ind =   i + (j-1)*na;
        str    = ['S',num2str(ind),' = ']; fprintf(fid, str);

        for ip = 1:na
            for jp=1:ns
             indp =   ip + (jp-1)*na ;
               if (Vind(indp)==i)&(Trans(jp,j)>0) % goes from indp to ip (i+1 appears later)
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(jp,j),formatSpec),'*(',num2str(aGrid(Vind(indp)+1),formatSpec) ,'-ap',num2str(indp),'(-1))/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
               end
               if (Vind(indp)==i-1)&(Trans(jp,j)>0) % goes from indp to ip (i+1 appears later)
                    if aGrid(Vind(indp))>0
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(jp,j),formatSpec),'*(ap',num2str(indp),'(-1)-',num2str(aGrid(Vind(indp)),formatSpec) ,')/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
                    else
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(jp,j),formatSpec),'*(ap',num2str(indp),'(-1)' ,')/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
                    end
              end
            end
       end
         if pass==0
          str    = [num2str(D_ss(ind),formatSpec)]; fprintf(fid, str);
        end
       str    = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    end
end

%% AGGREGATIONS
str = [' %%%%%%%% equation ',num2str(equa),': Aggretations \n\n']; fprintf(fid, str);
str = ['Pop = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i)]; fprintf(fid, str);
      if mod(i,10)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;

str = ['K = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i),'(-1)*ap',num2str(i),'(-1)']; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%%%%%%%%%%%%%%%%%  TOTAL WELFARE  %%%%%%%%%%%%%%%%
str = ['W = (TT)^(theta) + ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i),'*(((c',num2str(i), ')^(1-gamma) -1)/(1-gamma) ) ']; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;

str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((R-1 + delta)/alpha)^(alpha/(alpha-1));']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['(R-1 + delta)/(alpha*Z) = (K/1)^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['C +I + TT = GDP;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Z = 1+u;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_u*u(-1)+eps;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDP = Z*(K^alpha)*(1^(1-alpha));']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;
str = ['I =  K(+1) - ( 1 - delta)*K;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['TT = TT_ss + aGDP*(GDP - ',num2str(Y,formatSpec),') + aC*(C - ',num2str(Ctot,formatSpec),');']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(C/steady_state(C)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Rt = 100*(R - ',num2str(R,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Tt = 100*(TT/',num2str(TT_ss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Zt = ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%str = ['Wt = 100*(exp((1-beta)*( W-',num2str(Ws,formatSpec),')) -1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Wt = 100*(exp(( W-',num2str(Ws,formatSpec),'))-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['end;\n\n'] ; fprintf(fid, str);

%% STEADY STATE

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);
str = ['\n\n'] ; fprintf(fid, str);
%str = ['initval;\n\n'] ; fprintf(fid, str);

str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
for i = 1:Ntot
    str = ['ap',num2str(i),' = ',num2str(apol(i),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(i), ' = ',num2str(Cpol(i),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['S',num2str(i), ' = ',num2str(D_ss(i),formatSpec),';\n'] ; fprintf(fid, str);
     str = ['ww',num2str(i), ' = ',num2str(Wind(i),formatSpec),';\n'] ; fprintf(fid, str);
end;
for i = 1:Ny
    str = ['l',num2str(i),' =  (',num2str(l0,formatSpec),')^',num2str(1),';']; fprintf(fid, str);
end;

str = ['R = ', num2str(R,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['TT = ', num2str(TT_ss,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z  = 1;\n\n']; fprintf(fid, str);
str = ['u  = 0;\n\n']; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Pop = 1;\n']; fprintf(fid, str);

str = ['I =  delta*K;\n']; fprintf(fid, str);
str = ['C = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['W = ', num2str(Ws,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['Rt = 0;\n'] ; fprintf(fid, str);
str = ['Tt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['It = 0;\n'] ; fprintf(fid, str);
str = ['Wt = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);

str = ['resid;\n\n'] ; fprintf(fid, str);
str = ['steady(nocheck);\n\n'] ; fprintf(fid, str);
str = ['check;\n\n'] ; fprintf(fid, str);

%str = ['options_.solve_tolf=10^-6;\n\n'] ; fprintf(fid, str);
%str = ['steady(maxit = 100,solve_algo=3);\n\n'] ; fprintf(fid, str);
%str = ['steady;\n\n'] ; fprintf(fid, str);

str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=10000) Zt GDPt Ct Kt It Wt Rt wt Tt Z GDP C K I W R w TT;\n\n'] ; fprintf(fid, str);
fflush(fid);

%% CALLING DYNARE
disp('Dynare file written. Lauching Dynare...')
dynare DyReiter.mod nostrict


%% SAVING RESULTS 
disp('Outcomes')
run post_dynare_stochsimul.m

%~ %% GENERATING THE RESULTS

%~ save tofigreiter Zt_eps GDPt_eps Ct_eps Kt_eps It_eps Wt_eps Rt_eps wt_eps Tt_eps
