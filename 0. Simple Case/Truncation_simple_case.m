% Presentation

% This code simulates the dynamics of a simple heterogeneous-agent model 
% with production. There are wto risks: an aggregate technology shock and 
% individual productivity risk. Agents inelastically supply one unit of 
% labor. They face a standard consumption-saving trade-off. The production 
% sector is standard and the production function is Cobb-Douglas.
% A detailed acccount of the model can be found in the paper "The 
% truncation method to solve heterogeneous-agent models" by François Le 
% Grand and Xavier Ragot. 


% Prequisites: 

% The following code requires Matlab/Octave + the Dynare package.
% Dynare can be downloaded here: https://www.dynare.org/

% # How to use it?

% This file is self-contained. Once run in Matlab/Octave, it produces the
% the output for the truncation method in roughly one second.
 
% # How to use it ? 
% Starting from an initial distribution of wealth, we first derive the 
% steady-state allocation and the truncated model for N=2.  

% We then write a Dynare file, which countaineq the equations of the
% dynamic model. This file in then run to get the outcome.

% Obviously, the algorithm could be accelerated, at the cost of a less readable
% code. For higher values of N, this in done in Julia in another
% file (Please see the README file for details).

clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters of the model       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Production function 
alpha   = 0.36;                     % Cobb-Dougls parametr
delta   = 0.025;                    % capital depreciation

% Preferences
beta    = 0.98;                     % discount factor
gamma   = 1.0;                      % inverse of IES

% Asset market
abar    = 0;                        % credit constraint

% Aggregate risk (AR(1))
rho_z   = 0.95;                     % Persistence
sigma_z = 0.00312;                  % innovation standard error

% Idiosyncratic risk

% NB: The following parameters are taken from a Rouwenhorst
% procedure. Roughly consistent model of US labor market.

%   # PI: Transition Matrix across states
PI = [0.9807    0.0048      0.0       0.0         0.0
      0.0193    0.9808      0.0096    0.0         0.0
      0.0       0.0144      0.9808    0.0144      0.0
      0.0       0.0         0.0096    0.9808      0.0193
      0.0       0.0         0.0       0.0048      0.9807];
%NB: PI is a Markov transition matrix. The columns of PI sum to 1. 
%CHECK: max(abs(sum(PI) - ones(1,ns))) < 1e-12
  
%   # Productivity by states
states = [0.1810
    0.3740
    0.7730
    1.5976
    3.3023];

%   # Vector of indices per history
ns = length(states);
y0 = kron([1:ns],ones(1,ns)); 
% NB: y0 is a vector of current productivity index by history. 
% y0(h) gives the productivity level index of history h
% The order is such that the % 1st five histories have the lowest 
% productivity level, the 2nd set of five has the second lowest 
% productivity level and so on and so forth.


% Truncation parameters
N    = 2; % Truncation length
Ntot = ns^N; %Total number of histories

 
%Initial distribution of wealth
ae = [0 8.6897 0 0 0 0 8.6899 19.3978 0 0 0 8.7613 19.5136 48.9102 0 0 0 20.0354 49.4560  114.2717 0  0 0 50.9871 115.7951]';
% NB: at least one history (with positive size) must be credit-constrained.
% CHECK: (length(S(ae==0.0)) > 1) && (max(S(ae==0.0)) > 0.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing the truncated model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIh: transition matrix across histories.
PIh = zeros(Ntot,Ntot);  
for h = 1:Ntot
    for y = 1:ns
        hp = y0(h) + (y-1)*ns;      % this is the index of a possible continuation history
                                    % Indeed any possible continuation history is of the form (y0(h), y) 
                                    %   with any possible y -> index = y0(h) + (y-1)*ns
        PIh(hp,h) = PI(y,y0(h));    % Proba of h to hp is proba of y0(h) to y
    end
end
  
% Steady-state distribution across histories
S = (PIh^10000);
S=S(:,1); 

% NB: S verifies S = PIh*S 
% CHECK: max(abs(S-PIh*S) < 1e-12)

% NB: aggregate labor supply sums to 1 efficient unit/
% CHECK: abs(S'*states(y0)-1) < 1e-12

%real interest rates and wage
K = S'*ae;                          %financial market equilibrium
rss = alpha*K^(alpha-1)-delta;      %real interest rate (recall that total labor = 1)
wss = (1-alpha)*K^alpha;            % real wage rate

% Computing the allocation
a_tilde                 = (1./S).*(PIh*(S.*ae));    
                                    % per-capita beginning of period wealth
                                    % See paper for formula
a_tilde(abs(S) < 1e-12) = 0;        % remove the history vith size 0
c                       = (1+rss)*a_tilde - ae + wss*states(y0); 
                                    % vector of consumption
Ctot                    = S'*c;     % total consumption
Y                       = K^alpha;  % GDP (as labor is 1)
I                       = delta*K;  % investment

% NB: S is the steady-state distribution across histories
% CHECK: max(abs(S.*a_tilde - PIh*(S.*ae))) < 1e-12 

% Computing the xsis

P   = diag(ae > 0.0);               % diagnonal matrix with 1 on index h is h is unconstrained
Id  = eye(Ntot);                    % Identity matrix
xsi = (P*(Id - beta*(1+rss)*PIh')*diag(c.^-gamma)+Id-P)\((Id-P)*ones(Ntot,1));
                                    % Computation of the xsis 
% NB: We use left division rather than matrix inversion since it is more efficient.
% CHECK: max(abs(P*(xsi.*(c.^-gamma) - beta*(1+rss)*PIh'*(xsi.*(c.^-gamma))))) < 1e-12
% CHECK: max(abs((P*(Id - beta*(1+rss)*PIh')*diag(c.^-gamma)+Id-P)*xsi - (Id-P)*ones(Ntot,1))) < 1e-12


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing the Dynare file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this part, we write a Dynare file 'DyTruncation.mod' to use Dynare
% to solve for the model with perturbation. Other aspects of this solver
% could be used (estimation, high-order perturbation, etc.).

% To write the file, we define strings that are printed in the file. 

fid = fopen('DyTruncation.mod', 'w');  
                                    % opening the file with wirting permissions
formatSpec = '%20.12f';             % format for writing doubles

% Defining endogeneous variables
str = ['var\n']; fprintf(fid, str);
str = ['r w Z K GDP u I C\n']; fprintf(fid, str);  
str = ['rt wt Zt Kt GDPt ut It Ct  \n']; fprintf(fid, str);

for p = 1:Ntot % loop to cover all history-dependent variables
    str = ['a',num2str(p),' ']; fprintf(fid,str);      % end-of-period savings
    str = ['at',num2str(p),' ']; fprintf(fid, str);    % beginning-of-period savings
    str = ['c',num2str(p),' ']; fprintf(fid, str);     % consumption
    if mod(p,10)==0
    % every 10 history, input a newline (10 is an arbitrary choice). This is only 
    % for readibility purpose of the mod file.
        str = ['\n'];         fprintf(fid, str);
    end
end
str = ['; \n\n']; fprintf(fid, str);% Close the definition of endogeneous vars

% Defining exogeneous variables
str = ['varexo eps; \n\n']; fprintf(fid, str);
                                    % eps is the aggregate shock innovation

% Defining parameters
    % They have exactly the same name / same values as those in the Matlab file
str = ['parameters\n']; fprintf(fid, str);
    % Declaration
str = ['beta alpha abar delta gamma rho_z;\n\n']; fprintf(fid, str); 
    % Initialisation
str = ['alpha','   = ',num2str(alpha),';\n']; fprintf(fid, str);
str = ['beta','    = ',num2str(beta),';\n']; fprintf(fid, str);
str = ['abar','    = ',num2str(abar),';\n']; fprintf(fid, str);
str = ['delta','   = ',num2str(delta),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(gamma),';\n']; fprintf(fid, str);
str = ['rho_z','   = ',num2str(rho_z),';\n']; fprintf(fid, str); 

% Writing the model equations
equa    = 1; %index for the numebr of equation
tol     = 10^-10;  % hreshold for transitions. Avoids  considering very small transitions

%~ str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']; fprintf(fid, str);
str = ['\n\n']; fprintf(fid, str);
str = ['model;\n\n']; fprintf(fid, str);

for h = 1:Ntot   % For each history...

%%%%%%%%%%%%%%%%% Budget Constraint 
        str = ['c',num2str(h),' = ',num2str((states(y0(h))),formatSpec), '*w + (1+r)*','at',num2str(h),'-a', num2str(h),';']; fprintf(fid, str);
        str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; % PRinting the equation number to ease references in the mod file

%%%%%%%%%%%%%%%%% Definition of a_tilde (begiining-of-period wealth)
        str = ['at',num2str(h),' = 10^-10 ']; fprintf(fid, str);
        for hi = 1:Ntot
            if (S(hi)*PIh(h,hi)) > tol   % hi to h
              str = ['+ ',num2str(S(hi)*PIh(h,hi)/S(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;        

%%%%%%%%%%%%%%%%% Euler equation
        if P(h,h)==1 % not constrained
            str = [num2str(xsi(h),formatSpec),'*(c',num2str(h),'^-gamma) = beta*(1+r(+1))*(']; fprintf(fid, str);
            for hp = 1:Ntot
                if (PIh(hp,h)*xsi(hp))>tol % h to hp
                    str = ['+ ',num2str((PIh(hp,h))*xsi(hp),formatSpec),'*(c',num2str(hp),'(+1)^-gamma)']; fprintf(fid, str);
                end
            end
            str = [');']; fprintf(fid, str);
            str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['a',num2str(h),' =  steady_state(a',num2str(h),');']; fprintf(fid, str);
            str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end
end
  
% Aggregate quantities

%%% Wage rate
str = ['w = (1 - alpha)*(Z^(1/(1-alpha)))*((r + delta)/alpha)^(alpha/(alpha-1));']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%% Interest rate
str = ['r  = (alpha*Z)*(K(-1)/1)^(alpha-1) - delta;']; fprintf(fid, str); 
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%% Aggregate shock
str = ['Z = 1 + u;']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_z*u(-1) + eps;']; fprintf(fid, str);    
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%% GDP
str = ['GDP = Z*K(-1)^alpha*1^(1-alpha);']; fprintf(fid, str);      
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%% Aggregate consumption
str = ['\n']; fprintf(fid, str);
str = ['C =']; fprintf(fid, str);
for h = 1:Ntot % For all histories
    str = ['+',num2str(S(h),formatSpec),'*c',num2str(h)]; fprintf(fid, str);
    if mod(h,5) == 0 % break line once in a while (5 is arbitrary)
        fprintf(fid, '\n');
    end
end
str = [';']; fprintf(fid, str);
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%% Investment
str = ['I =  K - (1 - delta)*K(-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   

%%% Capital
str = ['K = ']; fprintf(fid, str);
for h = 1:Ntot
    str = ['+', num2str(S(h),formatSpec),'*a',num2str(h)]; fprintf(fid, str);
    if mod(h,5) == 0 % break line once in a while (5 is arbitrary)
        fprintf(fid, '\n');
    end
end
str = [';']; fprintf(fid, str); 
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%% Defining variables in relative deviation from steday-state value
% Same notation as before but with a terminal 't'
str = ['ut = 100*u;']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Ct = 100*(C/steady_state(C)-1);']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['wt = 100*(w/',num2str(wss,formatSpec),'-1);']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['rt = 100*(r - ',num2str(rss,formatSpec),');']; fprintf(fid, str); 
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Zt = ut;']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['It = 100*(I/',num2str(I,formatSpec),'-1);']; fprintf(fid, str);  
str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['\n']; fprintf(fid, str);
str = ['end;\n\n']; fprintf(fid, str);

% Steady State

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']; fprintf(fid, str);
str = ['\n\n']; fprintf(fid, str);

% Defining initial values
str = ['initval;\n\n']; fprintf(fid, str);
for h = 1:Ntot
    str = ['a',num2str(h),'  = ', num2str(ae(h),formatSpec),';\n']; fprintf(fid, str);
    str = ['at',num2str(h),' = ', num2str(a_tilde(h),formatSpec),';\n']; fprintf(fid, str);
    str = ['c',num2str(h),'  = ', num2str(c(h),formatSpec),';\n']; fprintf(fid, str);
end
str = ['r = ', num2str(rss,formatSpec),';\n']; fprintf(fid, str);
str = ['w = ', num2str(wss,formatSpec),';\n']; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n']; fprintf(fid, str);
str = ['Z = ', num2str(1),';\n']; fprintf(fid, str);
str = ['u = ', num2str(0),';\n']; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n']; fprintf(fid, str);
str = ['I = ', num2str(I,formatSpec),';\n']; fprintf(fid, str);
str = ['C = ', num2str(Ctot,formatSpec),';\n']; fprintf(fid, str);
str = ['ut = 0;\n']; fprintf(fid, str);
str = ['Kt = 0;\n']; fprintf(fid, str);
str = ['Ct  = 0;\n']; fprintf(fid, str);
str = ['wt = 0;\n']; fprintf(fid, str);
str = ['rt = 0;\n']; fprintf(fid, str);
str = ['Zt = 0;\n']; fprintf(fid, str);
str = ['GDPt = 0;\n']; fprintf(fid, str);
str = ['It = 0;\n']; fprintf(fid, str);
str = ['end;\n\n']; fprintf(fid, str);

% Check that the residuals of all equations are close to 0
str = ['resid;\n\n']; fprintf(fid, str);

% Defining the shocks and lauching stoch_simul for 10000 periods and 200 for IRFs
str = ['shocks;\n var eps; stderr ',num2str(sigma_z,formatSpec),';\n end;\n']; fprintf(fid, str);
str = ['stoch_simul (order=1, irf=200, periods=10000) ut Zt GDPt Ct Kt It rt wt u Z GDP C  K I r w;']; fprintf(fid, str);

% Calling Dynare
dynare DyTruncation noclearall

% NB: The code is a simplified version of the one used in the paper (to get the initial vector ae). So, the results 
% might slightly differ from the ones in the paper.
