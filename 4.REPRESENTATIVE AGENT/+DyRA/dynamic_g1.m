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
    T = DyRA.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(22, 29);
g1(1,4)=(-y(1));
g1(1,5)=(-1);
g1(1,1)=(-y(4));
g1(1,7)=1;
g1(1,12)=1;
g1(1,15)=1;
g1(2,26)=(-(params(1)*T(1)));
g1(2,12)=getPowerDeriv(y(12),(-params(6)),1);
g1(2,28)=(-(params(1)*y(26)*getPowerDeriv(y(28),(-params(6)),1)));
g1(3,7)=1;
g1(3,8)=(-1);
g1(4,4)=(-(T(2)*1/params(2)*getPowerDeriv((y(4)-1+params(5))/params(2),params(2)/(params(2)-1),1)));
g1(4,5)=1;
g1(4,6)=(-(T(3)*(1-params(2))*getPowerDeriv(y(6),1/(1-params(2)),1)));
g1(5,4)=1/(params(2)*y(6));
g1(5,6)=(-(params(2)*(y(4)-1+params(5))))/(params(2)*y(6)*params(2)*y(6));
g1(5,8)=(-(getPowerDeriv(y(8),params(2)-1,1)));
g1(6,3)=(-params(7));
g1(6,10)=1;
g1(6,29)=(-1);
g1(7,6)=1;
g1(7,10)=(-1);
g1(8,6)=(-T(4));
g1(8,2)=(-(y(6)*getPowerDeriv(y(2),params(2),1)));
g1(8,9)=1;
g1(9,8)=1-params(5);
g1(9,27)=(-1);
g1(9,11)=1;
g1(10,9)=(-1);
g1(10,11)=1;
g1(10,13)=1;
g1(10,15)=1;
g1(11,8)=(-params(13));
g1(11,9)=(-params(12));
g1(11,15)=1;
g1(12,12)=(-(getPowerDeriv(y(12),1-params(6),1)/(1-params(6))));
g1(12,14)=1;
g1(12,15)=(-(getPowerDeriv(y(15),params(3),1)));
g1(13,10)=(-100);
g1(13,21)=1;
g1(14,8)=(-2.46365442868034);
g1(14,19)=1;
g1(15,12)=(-40.40611736010123);
g1(15,23)=1;
g1(16,5)=(-41.18971160253174);
g1(16,17)=1;
g1(17,4)=(-100);
g1(17,16)=1;
g1(18,18)=1;
g1(18,21)=(-1);
g1(19,9)=(-26.36141539521738);
g1(19,20)=1;
g1(20,11)=(-98.54617683159476);
g1(20,22)=1;
g1(21,14)=(-(100*exp(y(14)-11.46959010)));
g1(21,24)=1;
g1(22,15)=(-329.1665998732775);
g1(22,25)=1;

end
