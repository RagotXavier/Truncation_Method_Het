function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 24);

T = DyTruncation.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(14) = (-(params(1)*(1+y(107))*0.005494799726*getPowerDeriv(y(109),(-params(5)),1)));
T(15) = (-(params(1)*(1+y(107))*0.630960676273*getPowerDeriv(y(111),(-params(5)),1)));
T(16) = (-(params(1)*(1+y(107))*0.008286174858*getPowerDeriv(y(112),(-params(5)),1)));
T(17) = (-(params(1)*(1+y(107))*0.009943940773*getPowerDeriv(y(113),(-params(5)),1)));
T(18) = (-(params(1)*(1+y(107))*0.794064307420*getPowerDeriv(y(114),(-params(5)),1)));
T(19) = (-(params(1)*(1+y(107))*0.015675483060*getPowerDeriv(y(115),(-params(5)),1)));
T(20) = (-(params(1)*(1+y(107))*0.009359594633*getPowerDeriv(y(116),(-params(5)),1)));
T(21) = (-(params(1)*(1+y(107))*1.158427774228*getPowerDeriv(y(117),(-params(5)),1)));
T(22) = (-(params(1)*(1+y(107))*0.030878177096*getPowerDeriv(y(118),(-params(5)),1)));
T(23) = (-(params(1)*(1+y(107))*0.007212844108*getPowerDeriv(y(119),(-params(5)),1)));
T(24) = (-(params(1)*(1+y(107))*1.820260607621*getPowerDeriv(y(120),(-params(5)),1)));

end
