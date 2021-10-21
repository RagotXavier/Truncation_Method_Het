function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = DyRa.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(18, 24);
g1(1,3)=(-y(1));
g1(1,4)=(-1);
g1(1,1)=(-(1+y(3)));
g1(1,6)=1;
g1(1,10)=1;
g1(2,21)=(-(params(1)*T(1)));
g1(2,10)=getPowerDeriv(y(10),(-params(5)),1);
g1(2,23)=(-(params(1)*(1+y(21))*getPowerDeriv(y(23),(-params(5)),1)));
g1(3,4)=1;
g1(3,5)=(-((1-params(2))*T(2)));
g1(3,6)=(-((1-params(2))*y(5)*getPowerDeriv(y(6),params(2),1)));
g1(4,3)=1;
g1(4,5)=(-(params(2)*T(3)));
g1(4,1)=(-(params(2)*y(5)*getPowerDeriv(y(1),params(2)-1,1)));
g1(5,5)=1;
g1(5,8)=(-1);
g1(6,2)=(-params(6));
g1(6,8)=1;
g1(6,24)=(-1);
g1(7,5)=(-T(4));
g1(7,1)=(-(y(5)*getPowerDeriv(y(1),params(2),1)));
g1(7,7)=1;
g1(8,6)=1-params(4);
g1(8,22)=(-1);
g1(8,9)=1;
g1(9,10)=(-(getPowerDeriv(y(10),1-params(5),1)/(1-params(5))));
g1(9,11)=1;
g1(10,8)=(-100);
g1(10,17)=1;
g1(11,6)=(-3.93594586439602);
g1(11,15)=1;
g1(12,10)=(-(100*1/(steady_state(8))));
g1(12,19)=1;
g1(13,4)=(-48.75708226483609);
g1(13,13)=1;
g1(14,3)=(-100);
g1(14,12)=1;
g1(15,14)=1;
g1(15,17)=(-1);
g1(16,7)=(-31.20453255820859);
g1(16,16)=1;
g1(17,9)=(-157.4378340181407);
g1(17,18)=1;
g1(18,11)=(-100);
g1(18,20)=1;

end
