function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 19);

T(1) = y(22)^(-params(5));
T(2) = 0.980700000000*y(19)^(-params(5))+0.019300000000*y(34)^(-params(5));
T(3) = y(37)^(-params(5));
T(4) = y(52)^(-params(5));
T(5) = (1+y(1))*params(1)*(T(1)*0.005494799726+T(3)*0.630960676273+0.009943940773*T(4));
T(6) = y(40)^(-params(5));
T(7) = y(55)^(-params(5));
T(8) = y(70)^(-params(5));
T(9) = (1+y(1))*params(1)*(T(6)*0.008286174858+0.794064307420*T(7)+0.009359594633*T(8));
T(10) = y(58)^(-params(5));
T(11) = y(73)^(-params(5));
T(12) = y(88)^(-params(5));
T(13) = (1+y(1))*params(1)*(T(10)*0.015675483060+1.158427774228*T(11)+0.007212844108*T(12));
T(14) = y(76)^(-params(5));
T(15) = y(91)^(-params(5));
T(16) = (1-params(2))*y(3)^(1/(1-params(2)));
T(17) = ((y(1)+params(4))/params(2))^(params(2)/(params(2)-1));
T(18) = y(4)^(params(2)-1);
T(19) = y(4)^params(2);

end
