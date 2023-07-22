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

assert(length(T) >= 24);

T(1) = y(46)*0.00359063+y(58)*0.00001753+y(28)*0.00119685+y(43)*0.00000870+y(25)*0.06084884+y(31)*0.00118815;
T(2) = 0.0668507;
T(3) = y(34)^(-params(5));
T(4) = y(46)^(-params(5));
T(5) = y(58)^(-params(5));
T(6) = 1e-10+0.00478739*y(28)^(-params(5))+T(3)*0.52603191+0.00802544*T(4)+0.00004976*T(5);
T(7) = (1+y(1))*params(1)*T(6);
T(8) = y(37)^(-params(5));
T(9) = y(40)^(-params(5));
T(10) = y(49)^(-params(5));
T(11) = y(61)^(-params(5));
T(12) = (1+y(1))*params(1)*(1e-10+T(8)*0.00671861+0.63307748*T(10)+0.00737431*T(11));
T(13) = y(52)^(-params(5));
T(14) = y(55)^(-params(5));
T(15) = y(64)^(-params(5));
T(16) = y(70)^(-params(5));
T(17) = (1+y(1))*params(1)*(1e-10+T(9)*0.00007537+T(13)*0.01240083+0.91152873*T(15)+0.00562809*T(16));
T(18) = 0.00000000+T(17);
T(19) = y(67)^(-params(5));
T(20) = y(73)^(-params(5));
T(21) = (1-params(2))*y(3)^(1/(1-params(2)));
T(22) = ((y(1)+params(4))/params(2))^(params(2)/(params(2)-1));
T(23) = y(5)^params(2);
T(24) = 1-params(5);

end
