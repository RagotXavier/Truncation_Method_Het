function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = DyRa.static_g1_tt(T, y, x, params);
end
g1 = zeros(18, 18);
g1(1,1)=(-y(4));
g1(1,2)=(-1);
g1(1,4)=1-(1+y(1));
g1(1,8)=1;
g1(2,1)=(-(T(1)*params(1)));
g1(2,8)=T(5)-(1+y(1))*params(1)*T(5);
g1(3,2)=1;
g1(3,3)=(-((1-params(2))*T(2)));
g1(3,4)=(-((1-params(2))*y(3)*T(4)));
g1(4,1)=1;
g1(4,3)=(-(params(2)*T(3)));
g1(4,4)=(-(params(2)*y(3)*getPowerDeriv(y(4),params(2)-1,1)));
g1(5,3)=1;
g1(5,6)=(-1);
g1(6,6)=1-params(6);
g1(7,3)=(-T(2));
g1(7,4)=(-(y(3)*T(4)));
g1(7,5)=1;
g1(8,4)=(-(1-(1-params(4))));
g1(8,7)=1;
g1(9,8)=(-(getPowerDeriv(y(8),1-params(5),1)/(1-params(5))));
g1(9,9)=1;
g1(10,6)=(-100);
g1(10,15)=1;
g1(11,4)=(-3.93594586439602);
g1(11,13)=1;
g1(12,8)=(-(100*((y(8))-y(8))/((y(8))*(y(8)))));
g1(12,17)=1;
g1(13,2)=(-48.75708226483609);
g1(13,11)=1;
g1(14,1)=(-100);
g1(14,10)=1;
g1(15,12)=1;
g1(15,15)=(-1);
g1(16,5)=(-31.20453255820859);
g1(16,14)=1;
g1(17,7)=(-157.4378340181407);
g1(17,16)=1;
g1(18,9)=(-100);
g1(18,18)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
