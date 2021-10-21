function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = DyRa.static_resid_tt(T, y, x, params);
end
residual = zeros(18, 1);
lhs = y(8)+y(4);
rhs = y(4)*(1+y(1))+y(2);
residual(1) = lhs - rhs;
lhs = T(1);
rhs = T(1)*(1+y(1))*params(1);
residual(2) = lhs - rhs;
lhs = y(2);
rhs = (1-params(2))*y(3)*T(2);
residual(3) = lhs - rhs;
lhs = y(1)+params(4);
rhs = params(2)*y(3)*T(3);
residual(4) = lhs - rhs;
lhs = y(3);
rhs = 1+y(6);
residual(5) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(6)+x(1);
residual(6) = lhs - rhs;
lhs = y(5);
rhs = y(3)*T(2);
residual(7) = lhs - rhs;
lhs = y(7);
rhs = y(4)-y(4)*(1-params(4));
residual(8) = lhs - rhs;
lhs = y(9);
rhs = (y(8)^(1-params(5))-1)/(1-params(5))+10;
residual(9) = lhs - rhs;
lhs = y(15);
rhs = y(6)*100;
residual(10) = lhs - rhs;
lhs = y(13);
rhs = 100*(y(4)/25.40685351-1);
residual(11) = lhs - rhs;
lhs = y(17);
rhs = 100*(y(8)/(y(8))-1);
residual(12) = lhs - rhs;
lhs = y(11);
rhs = 100*(y(2)/2.05098409-1);
residual(13) = lhs - rhs;
lhs = y(10);
rhs = 100*(y(1)-0.02040816);
residual(14) = lhs - rhs;
lhs = y(12);
rhs = y(15);
residual(15) = lhs - rhs;
lhs = y(14);
rhs = 100*(y(5)/3.20466265-1);
residual(16) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(7)/0.63517134-1);
residual(17) = lhs - rhs;
lhs = y(18);
rhs = 100*(y(9)-10.94366342);
residual(18) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
