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
    T = DyRA.static_g1_tt(T, y, x, params);
end
g1 = zeros(22, 22);
g1(1,1)=(-y(4));
g1(1,2)=(-1);
g1(1,4)=(-(y(1)-1));
g1(1,9)=1;
g1(1,12)=1;
g1(2,1)=(-(T(1)*params(1)));
g1(2,9)=T(5)-y(1)*params(1)*T(5);
g1(3,4)=1;
g1(3,5)=(-1);
g1(4,1)=(-(T(2)*1/params(2)*getPowerDeriv((y(1)-1+params(5))/params(2),params(2)/(params(2)-1),1)));
g1(4,2)=1;
g1(4,3)=(-(T(3)*(1-params(2))*getPowerDeriv(y(3),1/(1-params(2)),1)));
g1(5,1)=1/(params(2)*y(3));
g1(5,3)=(-(params(2)*(y(1)-1+params(5))))/(params(2)*y(3)*params(2)*y(3));
g1(5,5)=(-(getPowerDeriv(y(5),params(2)-1,1)));
g1(6,7)=1-params(7);
g1(7,3)=1;
g1(7,7)=(-1);
g1(8,3)=(-T(4));
g1(8,5)=(-(y(3)*getPowerDeriv(y(5),params(2),1)));
g1(8,6)=1;
g1(9,5)=(-(1-(1-params(5))));
g1(9,8)=1;
g1(10,6)=(-1);
g1(10,8)=1;
g1(10,10)=1;
g1(10,12)=1;
g1(11,5)=(-params(13));
g1(11,6)=(-params(12));
g1(11,12)=1;
g1(12,9)=(-(getPowerDeriv(y(9),1-params(6),1)/(1-params(6))));
g1(12,11)=1;
g1(12,12)=(-(getPowerDeriv(y(12),params(3),1)));
g1(13,7)=(-100);
g1(13,18)=1;
g1(14,5)=(-2.46365442868034);
g1(14,16)=1;
g1(15,9)=(-40.40611736010123);
g1(15,20)=1;
g1(16,2)=(-41.18971160253174);
g1(16,14)=1;
g1(17,1)=(-100);
g1(17,13)=1;
g1(18,15)=1;
g1(18,18)=(-1);
g1(19,6)=(-26.36141539521738);
g1(19,17)=1;
g1(20,8)=(-98.54617683159476);
g1(20,19)=1;
g1(21,11)=(-(100*exp(y(11)-11.46959010)));
g1(21,21)=1;
g1(22,12)=(-329.1665998732775);
g1(22,22)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
