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

assert(length(T) >= 30);

T = DyTruncation.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(17) = (-(params(1)*(1+y(94))*0.00478739*getPowerDeriv(y(96),T(3),1)));
T(18) = (-(params(1)*(1+y(94))*0.52603191*getPowerDeriv(y(97),T(3),1)));
T(19) = (-(params(1)*(1+y(94))*0.00671861*getPowerDeriv(y(98),T(3),1)));
T(20) = (-(params(1)*(1+y(94))*0.00007537*getPowerDeriv(y(99),T(3),1)));
T(21) = (-(params(1)*(1+y(94))*0.00802544*getPowerDeriv(y(100),T(3),1)));
T(22) = (-(params(1)*(1+y(94))*0.63307748*getPowerDeriv(y(101),T(3),1)));
T(23) = (-(params(1)*(1+y(94))*0.01240083*getPowerDeriv(y(102),T(3),1)));
T(24) = (-(params(1)*(1+y(94))*0.00018170*getPowerDeriv(y(103),T(3),1)));
T(25) = (-(params(1)*(1+y(94))*0.00004976*getPowerDeriv(y(104),T(3),1)));
T(26) = (-(params(1)*(1+y(94))*0.00737431*getPowerDeriv(y(105),T(3),1)));
T(27) = (-(params(1)*(1+y(94))*0.91152873*getPowerDeriv(y(106),T(3),1)));
T(28) = (-(params(1)*(1+y(94))*0.02407548*getPowerDeriv(y(107),T(3),1)));
T(29) = (-(params(1)*(1+y(94))*0.00562809*getPowerDeriv(y(108),T(3),1)));
T(30) = (-(params(1)*(1+y(94))*1.42199333*getPowerDeriv(y(109),T(3),1)));

end
