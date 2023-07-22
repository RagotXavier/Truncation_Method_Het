function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
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
%   residual
%

if T_flag
    T = DyRA.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(22, 1);
lhs = y(12);
rhs = y(5)+y(4)*y(1)-y(7)-y(15);
residual(1) = lhs - rhs;
lhs = y(12)^(-params(6));
rhs = params(1)*y(26)*T(1);
residual(2) = lhs - rhs;
lhs = y(7);
rhs = y(8);
residual(3) = lhs - rhs;
lhs = y(5);
rhs = T(2)*T(3);
residual(4) = lhs - rhs;
lhs = (y(4)-1+params(5))/(params(2)*y(6));
rhs = y(8)^(params(2)-1);
residual(5) = lhs - rhs;
lhs = y(10);
rhs = params(7)*y(3)+x(it_, 1);
residual(6) = lhs - rhs;
lhs = y(6);
rhs = 1+y(10);
residual(7) = lhs - rhs;
lhs = y(9);
rhs = y(6)*T(4);
residual(8) = lhs - rhs;
lhs = y(11);
rhs = y(27)-y(8)*(1-params(5));
residual(9) = lhs - rhs;
lhs = y(13);
rhs = y(9)-y(11)-y(15);
residual(10) = lhs - rhs;
lhs = y(15);
rhs = params(11)+y(9)*params(12)+y(8)*params(13);
residual(11) = lhs - rhs;
lhs = y(14);
rhs = (y(12)^(1-params(6))-1)/(1-params(6))+10+y(15)^params(3);
residual(12) = lhs - rhs;
lhs = y(21);
rhs = y(10)*100;
residual(13) = lhs - rhs;
lhs = y(19);
rhs = 100*(y(8)/40.59010827-1);
residual(14) = lhs - rhs;
lhs = y(23);
rhs = 100*(y(12)/2.47487278-1);
residual(15) = lhs - rhs;
lhs = y(17);
rhs = 100*(y(5)/2.42779073-1);
residual(16) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(4)-1.00864446);
residual(17) = lhs - rhs;
lhs = y(18);
rhs = y(21);
residual(18) = lhs - rhs;
lhs = y(20);
rhs = 100*(y(9)/3.79342302-1);
residual(19) = lhs - rhs;
lhs = y(22);
rhs = 100*(y(11)/1.01475271-1);
residual(20) = lhs - rhs;
lhs = y(24);
rhs = 100*exp(y(14)-11.46959010);
residual(21) = lhs - rhs;
lhs = y(25);
rhs = 100*(y(15)/0.30379753-1);
residual(22) = lhs - rhs;

end
