
%% CALLING THE MAIN VARIABLES

clear
clc
warning('off','all');
load ../todynare_Reiter.mat 


Ntot = na*ns; 
Ny = ns; % idiosyncratic shocks
thres       = 0; % threshold to consider 0 (not to have too many equations)
Y = K^alpha*Ltot^(1-alpha); % GDP
I = delta*K; % investment

%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u       = 0.95;
sigma_u     = 0.00312;

fid         = fopen('DyReiter.mod', 'w') ;
formatSpec  = '%16.12f';

%% DEFINING THE ENDOGENOUS VARIABLES 

str = ['var\n'] ; fprintf(fid, str);
str = ['r w Z K GDP u Pop I C Cb W\n'] ; fprintf(fid, str);  % variables
str = ['rt wt Zt Kt GDPt ut It Ct Cbt Wt\n'] ; fprintf(fid, str);

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
str         = ['; \n\n'] ; fprintf(fid, str);

str         = ['varexo eps; \n\n'] ; fprintf(fid, str);
str         = ['parameters\n'] ; fprintf(fid, str);
str         = ['beta alpha abar delta gamma rho_u \n'] ; fprintf(fid, str);
str = ['; \n\n'] ; fprintf(fid, str);

l0          = Ltot;
abar        = abar;
 
str         = ['alpha','   = ',num2str(alpha),';\n']; fprintf(fid, str); 
str         = ['beta ','   = ',num2str(beta),';\n']; fprintf(fid, str); 
str         = ['abar ','   = ',num2str(abar),';\n']; fprintf(fid, str); 
str         = ['gamma','   = ',num2str(gamma),';\n']; fprintf(fid, str); 
str         = ['rho_u','   = ',num2str(rho_u),';\n']; fprintf(fid, str); 
str         = ['delta','   = ',num2str(delta),';\n']; fprintf(fid, str);
str         = ['l    ','   = ',num2str(l0),';\n']; fprintf(fid, str);
%str     = ['TT  ','   = ',num2str(TT),';\n']; fprintf(fid, str);

str = ['\n\n'] ; fprintf(fid, str);
yv = ytype; % vector of productivity type from 1 to ns*na (for budget constraint)
yvv = [1:1:ns]; % vector of productivity type from 1 to ns 
Th = 10^-8;
indc    = find(apol <= Th); % credit constraind agents
indnc = find(apol >Th); 

%% EQUATIONS

str     = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str); 
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
    str = ['c',num2str(i),'+ap',num2str(i), ' = w*',num2str(states(yv(i)),formatSpec),'*l',num2str(yv(i)),'+ (1+r)*',num2str(aGrid(i)) ,';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;     
end;
str = ['\n\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%  WEALTH TRANSITION  %%%%%%%%%%%%%%
for i = 1:Ntot 
     str = ['ap',num2str(i),' =  ww',num2str(i),'*',num2str(aGrid(Vind(i))),'+ (1 -ww',num2str(i),')*',num2str(aGrid(Vind(i)+1)),';']; fprintf(fid, str);
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
          str = ['(c',num2str(i),' - ','0)^-gamma  = beta*(1+r(+1))*('];fprintf(fid, str);
          for j = 1:ns % next period productivity status
                  if (Trans(j,yi)>threse) % does agent i reach agent j ? effective transition toward the bin j
                     str = ['+',num2str(Trans(j,yi)/1,formatSpec),'*( ww',num2str(i),'*c',num2str(Vind(i)+(j-1)*na ),'(+1)+(1-ww',num2str(i),')*c',num2str(Vind(i)+1+(j-1)*na),'(+1) )^-gamma']; fprintf(fid, str);                               
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
               if (Vind(indp)==i)&(Trans(j,jp)>0) % goes from indp to ip (i+1 appears later)
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(j,jp),formatSpec),'*(',num2str(aGrid(Vind(indp)+1),formatSpec) ,'-ap',num2str(indp),'(-1))/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
               end
               if (Vind(indp)==i-1)&(Trans(j,jp)>0) % goes from indp to ip (i+1 appears later)
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(j,jp),formatSpec),'*(ap',num2str(indp),'(-1)-',num2str(aGrid(Vind(indp)),formatSpec) ,')/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
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
str = ['W = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i),'*(((c',num2str(i), ')^(1- gamma) -1)/(1-gamma) + 10)']; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  
 
%%%%%%%%%%%%%%%%%%%%%%%%%  CONSUMPTION BOTTOM 10 PERCENT  %%%%%%%%%%%%%%%%
str = ['#Cba = (']; fprintf(fid, str);
for i = [b']
      str = ['+S',num2str(i),'(1)*c',num2str(i),'(1)']; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;
end;
str = [');'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  

str = ['#ba = ']; fprintf(fid, str);
for i = [b']
     str    = ['+S',num2str(i)]; fprintf(fid, str);
      if mod(i,10)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  

str = ['Cb = Cba/ba;']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
 

str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((r + delta)/alpha)^(alpha/(alpha-1));']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['(r + delta)/(alpha*Z) = (K/1)^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['Z = 1+u;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['u = rho_u*u(-1)+eps;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['GDP = Z*(K^alpha)*(1^(1-alpha));']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;   
str = ['I =  K(+1) - ( 1 - delta)*K;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['C  =  GDP - I;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(C/steady_state(C)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Cbt = 100*(Cb/',num2str(Cbs,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = 100*(r - ',num2str(R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Wt = 100*(W-',num2str(Ws,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

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
str = ['r = ', num2str(R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z  = 1;\n\n']; fprintf(fid, str);
str = ['u  = 0;\n\n']; fprintf(fid, str);
str = ['GDP = Z*K^alpha*1^(1-alpha);\n']; fprintf(fid, str);
str = ['Pop = 1;\n']; fprintf(fid, str);

str = ['I =  delta*K;\n']; fprintf(fid, str); 
str = ['C =  GDP - I;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Cb = ', num2str(Cbs,formatSpec),';\n'] ; fprintf(fid, str);
str = ['W = ', num2str(Ws,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['Cbt  = 0;\n'] ; fprintf(fid, str);
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
str = ['shocks;\n var eps; stderr ',num2str(sigma_u,formatSpec),';\n end;\n'] ; fprintf(fid, str);

%str = ['extended_path (periods=100, order=0);'] ; fprintf(fid, str);   
%str = ['stoch_simul (order=1,irf=100, hp_filter =1600)  u Z GDP C Cb Cm Ctp K I W r w;'] ; fprintf(fid, str);    
%str = ['stoch_simul (order=1,irf=100, hp_filter =1600)  u Z GDP C I w;'] ; fprintf(fid, str); 
%str = ['stoch_simul (order=1,irf=200, periods = 10000)  u Z GDP C K I r w;'] ; fprintf(fid, str);
str = ['options_.TeX=1;\n\n'] ; fprintf(fid, str);
str = ['stoch_simul (order=1,irf=200, periods=10000) ut Zt GDPt Ct Cbt Kt It Wt rt wt u Z GDP C Cb  K I W r w;'] ; fprintf(fid, str);

%% CALLING DYNARE

dynare DyReiter.mod nostrict

%% GENERATING THE RESULTS

save tofigreiter ut_eps Zt_eps GDPt_eps Ct_eps Cbt_eps Kt_eps It_eps Wt_eps rt_eps wt_eps 
%save totalreiter oo_

save tofigreitersim u Z GDP C Cb K I W r w 
save tofigreitersimt ut Zt GDPt Ct Cbt Kt It Wt rt wt 



