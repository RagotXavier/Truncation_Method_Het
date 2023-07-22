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
    T = DyRA.static_resid_tt(T, y, x, params);
end
residual = zeros(22, 1);
lhs = y(9);
rhs = y(2)+y(1)*y(4)-y(4)-y(12);
residual(1) = lhs - rhs;
lhs = T(1);
rhs = T(1)*y(1)*params(1);
residual(2) = lhs - rhs;
lhs = y(4);
rhs = y(5);
residual(3) = lhs - rhs;
lhs = y(2);
rhs = T(2)*T(3);
residual(4) = lhs - rhs;
lhs = (y(1)-1+params(5))/(params(2)*y(3));
rhs = y(5)^(params(2)-1);
residual(5) = lhs - rhs;
lhs = y(7);
rhs = y(7)*params(7)+x(1);
residual(6) = lhs - rhs;
lhs = y(3);
rhs = 1+y(7);
residual(7) = lhs - rhs;
lhs = y(6);
rhs = y(3)*T(4);
residual(8) = lhs - rhs;
lhs = y(8);
rhs = y(5)-y(5)*(1-params(5));
residual(9) = lhs - rhs;
lhs = y(10);
rhs = y(6)-y(8)-y(12);
residual(10) = lhs - rhs;
lhs = y(12);
rhs = params(11)+y(6)*params(12)+y(5)*params(13);
residual(11) = lhs - rhs;
lhs = y(11);
rhs = (y(9)^(1-params(6))-1)/(1-params(6))+10+y(12)^params(3);
residual(12) = lhs - rhs;
lhs = y(18);
rhs = y(7)*100;
residual(13) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(5)/40.59010827-1);
residual(14) = lhs - rhs;
lhs = y(20);
rhs = 100*(y(9)/2.47487278-1);
residual(15) = lhs - rhs;
lhs = y(14);
rhs = 100*(y(2)/2.42779073-1);
residual(16) = lhs - rhs;
lhs = y(13);
rhs = 100*(y(1)-1.00864446);
residual(17) = lhs - rhs;
lhs = y(15);
rhs = y(18);
residual(18) = lhs - rhs;
lhs = y(17);
rhs = 100*(y(6)/3.79342302-1);
residual(19) = lhs - rhs;
lhs = y(19);
rhs = 100*(y(8)/1.01475271-1);
residual(20) = lhs - rhs;
lhs = y(21);
rhs = 100*exp(y(11)-11.46959010);
residual(21) = lhs - rhs;
lhs = y(22);
rhs = 100*(y(12)/0.30379753-1);
residual(22) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
