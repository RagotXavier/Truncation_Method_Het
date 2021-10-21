function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 33);

T = DyTruncation.static_resid_tt(T, y, x, params);

T(20) = getPowerDeriv(y(22),(-params(5)),1);
T(21) = getPowerDeriv(y(37),(-params(5)),1);
T(22) = getPowerDeriv(y(40),(-params(5)),1);
T(23) = getPowerDeriv(y(52),(-params(5)),1);
T(24) = (-((1+y(1))*params(1)*0.009943940773*T(23)));
T(25) = getPowerDeriv(y(55),(-params(5)),1);
T(26) = getPowerDeriv(y(58),(-params(5)),1);
T(27) = getPowerDeriv(y(70),(-params(5)),1);
T(28) = (-((1+y(1))*params(1)*0.009359594633*T(27)));
T(29) = getPowerDeriv(y(73),(-params(5)),1);
T(30) = getPowerDeriv(y(76),(-params(5)),1);
T(31) = getPowerDeriv(y(88),(-params(5)),1);
T(32) = (-((1+y(1))*params(1)*0.007212844108*T(31)));
T(33) = getPowerDeriv(y(91),(-params(5)),1);

end
