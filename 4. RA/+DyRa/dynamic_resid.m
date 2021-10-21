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
    T = DyRa.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(18, 1);
lhs = y(10)+y(6);
rhs = (1+y(3))*y(1)+y(4);
residual(1) = lhs - rhs;
lhs = y(10)^(-params(5));
rhs = params(1)*(1+y(21))*T(1);
residual(2) = lhs - rhs;
lhs = y(4);
rhs = (1-params(2))*y(5)*T(2);
residual(3) = lhs - rhs;
lhs = y(3)+params(4);
rhs = params(2)*y(5)*T(3);
residual(4) = lhs - rhs;
lhs = y(5);
rhs = 1+y(8);
residual(5) = lhs - rhs;
lhs = y(8);
rhs = params(6)*y(2)+x(it_, 1);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = y(5)*T(4);
residual(7) = lhs - rhs;
lhs = y(9);
rhs = y(22)-y(6)*(1-params(4));
residual(8) = lhs - rhs;
lhs = y(11);
rhs = (y(10)^(1-params(5))-1)/(1-params(5))+10;
residual(9) = lhs - rhs;
lhs = y(17);
rhs = y(8)*100;
residual(10) = lhs - rhs;
lhs = y(15);
rhs = 100*(y(6)/25.40685351-1);
residual(11) = lhs - rhs;
lhs = y(19);
rhs = 100*(y(10)/(steady_state(8))-1);
residual(12) = lhs - rhs;
lhs = y(13);
rhs = 100*(y(4)/2.05098409-1);
residual(13) = lhs - rhs;
lhs = y(12);
rhs = 100*(y(3)-0.02040816);
residual(14) = lhs - rhs;
lhs = y(14);
rhs = y(17);
residual(15) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(7)/3.20466265-1);
residual(16) = lhs - rhs;
lhs = y(18);
rhs = 100*(y(9)/0.63517134-1);
residual(17) = lhs - rhs;
lhs = y(20);
rhs = 100*(y(11)-10.94366342);
residual(18) = lhs - rhs;

end
