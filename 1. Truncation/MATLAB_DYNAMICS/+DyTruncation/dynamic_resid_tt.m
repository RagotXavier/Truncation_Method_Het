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

assert(length(T) >= 16);

T(1) = y(66)*0.00359063+y(78)*0.00001753+y(48)*0.00119685+y(63)*0.00000870+y(45)*0.06084884+y(51)*0.00118815;
T(2) = 0.0668507;
T(3) = (-params(5));
T(4) = 1e-10+0.00478739*y(96)^T(3)+0.52603191*y(97)^T(3)+0.00802544*y(100)^T(3)+0.00004976*y(104)^T(3);
T(5) = params(1)*(1+y(94))*T(4);
T(6) = 1e-10+0.00671861*y(98)^T(3)+0.63307748*y(101)^T(3)+0.00737431*y(105)^T(3);
T(7) = params(1)*(1+y(94))*T(6);
T(8) = 1e-10+0.00007537*y(99)^T(3)+0.01240083*y(102)^T(3)+0.91152873*y(106)^T(3)+0.00562809*y(108)^T(3);
T(9) = params(1)*(1+y(94))*T(8);
T(10) = 0.00000000+T(9);
T(11) = 1e-10+0.00018170*y(103)^T(3)+0.02407548*y(107)^T(3)+1.42199333*y(109)^T(3);
T(12) = 0.00000000+params(1)*(1+y(94))*T(11);
T(13) = (1-params(2))*y(23)^(1/(1-params(2)));
T(14) = ((y(21)+params(4))/params(2))^(params(2)/(params(2)-1));
T(15) = y(2)^params(2);
T(16) = 1-params(5);

end
