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

assert(length(T) >= 44);

T = DyTruncation.static_resid_tt(T, y, x, params);

T(25) = (-((1+y(1))*params(1)*0.00478739*getPowerDeriv(y(28),(-params(5)),1)));
T(26) = getPowerDeriv(y(34),(-params(5)),1);
T(27) = getPowerDeriv(y(37),(-params(5)),1);
T(28) = getPowerDeriv(y(40),(-params(5)),1);
T(29) = getPowerDeriv(y(46),(-params(5)),1);
T(30) = (-((1+y(1))*params(1)*0.00802544*T(29)));
T(31) = getPowerDeriv(y(49),(-params(5)),1);
T(32) = (-((1+y(1))*params(1)*0.63307748*T(31)));
T(33) = getPowerDeriv(y(52),(-params(5)),1);
T(34) = getPowerDeriv(y(55),(-params(5)),1);
T(35) = getPowerDeriv(y(58),(-params(5)),1);
T(36) = (-((1+y(1))*params(1)*0.00004976*T(35)));
T(37) = getPowerDeriv(y(61),(-params(5)),1);
T(38) = (-((1+y(1))*params(1)*0.00737431*T(37)));
T(39) = getPowerDeriv(y(64),(-params(5)),1);
T(40) = (-((1+y(1))*params(1)*0.91152873*T(39)));
T(41) = getPowerDeriv(y(67),(-params(5)),1);
T(42) = getPowerDeriv(y(70),(-params(5)),1);
T(43) = (-((1+y(1))*params(1)*0.00562809*T(42)));
T(44) = getPowerDeriv(y(73),(-params(5)),1);

end
