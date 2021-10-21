function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 13);

T(1) = 0.980700000000*y(108)^(-params(5))+0.019300000000*y(110)^(-params(5));
T(2) = 0.005494799726*y(109)^(-params(5))+0.630960676273*y(111)^(-params(5))+0.009943940773*y(113)^(-params(5));
T(3) = params(1)*(1+y(107))*T(2);
T(4) = 0.008286174858*y(112)^(-params(5))+0.794064307420*y(114)^(-params(5))+0.009359594633*y(116)^(-params(5));
T(5) = params(1)*(1+y(107))*T(4);
T(6) = 0.015675483060*y(115)^(-params(5))+1.158427774228*y(117)^(-params(5))+0.007212844108*y(119)^(-params(5));
T(7) = params(1)*(1+y(107))*T(6);
T(8) = 0.030878177096*y(118)^(-params(5))+1.820260607621*y(120)^(-params(5));
T(9) = params(1)*(1+y(107))*T(8);
T(10) = (1-params(2))*y(18)^(1/(1-params(2)));
T(11) = ((y(16)+params(4))/params(2))^(params(2)/(params(2)-1));
T(12) = y(1)^(params(2)-1);
T(13) = y(1)^params(2);

end
