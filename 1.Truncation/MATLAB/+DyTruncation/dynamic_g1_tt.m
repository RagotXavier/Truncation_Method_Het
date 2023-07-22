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

assert(length(T) >= 89);

T = DyTruncation.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(42) = (-(params(1)*(0.062028126707*y(171)*y(169)+0.000471873293*y(167)*y(165)+0.000471873293*y(163)*y(161)+0.248112501707*y(159)*y(157)+0.001415625000*y(155)*y(153)+0.001415625000*y(151)*y(149)+y(147)*0.372168750000*y(145)+y(143)*0.001415625000*y(141)+y(139)*0.001415625000*y(137)+y(135)*0.248112501707*y(133)+y(130)*0.000471873293*y(128)+y(127)*0.062028126707*5.116037217023345+0.000471873293*y(132)*5.116037217023345)+params(1)*T(39)));
T(43) = (-(params(1)*(2.414874053818*0.062028126707*y(171)+0.000471873293*y(167)*2.414874053818+1.459493862568*0.000471873293*y(163)+1.459493862568*0.248112501707*y(159)+0.001415625000*y(155)*1.459493862568+0.882084236031*0.001415625000*y(151)+0.882084236031*y(147)*0.372168750000+y(143)*0.001415625000*0.882084236031+0.533111251380*y(139)*0.001415625000+0.533111251380*y(135)*0.248112501707+0.322200074254*y(130)*0.000471873293+0.322200074254*y(127)*0.062028126707+0.533111251380*0.000471873293*y(132))));
T(44) = getPowerDeriv(y(2),params(2)-1,1);
T(45) = (-(params(1)*0.062028126707*(0.322200074254*y(126)+y(125)*5.116037217023345)));
T(46) = getPowerDeriv(y(129),(-params(6)),1);
T(47) = (-(params(1)*y(125)*1.176900140040*0.000471873293*y(131)*T(46)));
T(48) = (-(params(1)*0.000471873293*(0.322200074254*y(126)+y(125)*y(128))));
T(49) = (-(params(1)*0.000471873293*(0.533111251380*y(126)+y(125)*5.116037217023345)));
T(50) = 1.063205188402*getPowerDeriv(y(66),(-params(6)),1);
T(51) = getPowerDeriv(y(134),(-params(6)),1);
T(52) = (-(params(1)*y(125)*1.063205188402*0.248112501707*y(136)*T(51)));
T(53) = (-(params(1)*0.248112501707*(0.533111251380*y(126)+y(125)*y(133))));
T(54) = 1.105333749474*getPowerDeriv(y(72),(-params(6)),1);
T(55) = getPowerDeriv(y(138),(-params(6)),1);
T(56) = (-(params(1)*y(125)*0.001415625000*y(140)*1.105333749474*T(55)));
T(57) = (-(params(1)*0.001415625000*(0.533111251380*y(126)+y(125)*y(137))));
T(58) = 1.021196458903*getPowerDeriv(y(78),(-params(6)),1);
T(59) = getPowerDeriv(y(142),(-params(6)),1);
T(60) = (-(params(1)*y(125)*0.001415625000*y(144)*1.021196458903*T(59)));
T(61) = (-(params(1)*0.001415625000*(y(125)*y(141)+y(126)*0.882084236031)));
T(62) = 1.042236122064*getPowerDeriv(y(84),(-params(6)),1);
T(63) = getPowerDeriv(y(146),(-params(6)),1);
T(64) = (-(params(1)*y(125)*0.372168750000*y(148)*1.042236122064*T(63)));
T(65) = (-(params(1)*0.372168750000*(y(126)*0.882084236031+y(125)*y(145))));
T(66) = 1.058263406894*getPowerDeriv(y(90),(-params(6)),1);
T(67) = getPowerDeriv(y(150),(-params(6)),1);
T(68) = (-(params(1)*y(125)*0.001415625000*y(152)*1.058263406894*T(67)));
T(69) = (-(params(1)*0.001415625000*(y(126)*0.882084236031+y(125)*y(149))));
T(70) = 1.018670598592*getPowerDeriv(y(96),(-params(6)),1);
T(71) = getPowerDeriv(y(154),(-params(6)),1);
T(72) = (-(params(1)*y(125)*0.001415625000*y(156)*1.018670598592*T(71)));
T(73) = (-(params(1)*0.001415625000*(y(125)*y(153)+y(126)*1.459493862568)));
T(74) = 1.029944405894*getPowerDeriv(y(102),(-params(6)),1);
T(75) = getPowerDeriv(y(158),(-params(6)),1);
T(76) = (-(params(1)*y(125)*0.248112501707*y(160)*1.029944405894*T(75)));
T(77) = (-(params(1)*0.248112501707*(y(126)*1.459493862568+y(125)*y(157))));
T(78) = 1.039669027202*getPowerDeriv(y(108),(-params(6)),1);
T(79) = getPowerDeriv(y(162),(-params(6)),1);
T(80) = (-(params(1)*y(125)*0.000471873293*y(164)*1.039669027202*T(79)));
T(81) = (-(params(1)*0.000471873293*(y(126)*1.459493862568+y(125)*y(161))));
T(82) = 1.015225620707*getPowerDeriv(y(114),(-params(6)),1);
T(83) = getPowerDeriv(y(166),(-params(6)),1);
T(84) = (-(params(1)*y(125)*0.000471873293*y(168)*1.015225620707*T(83)));
T(85) = (-(params(1)*0.000471873293*(y(125)*y(165)+y(126)*2.414874053818)));
T(86) = 1.022110952269*getPowerDeriv(y(120),(-params(6)),1);
T(87) = getPowerDeriv(y(170),(-params(6)),1);
T(88) = (-(params(1)*y(125)*0.062028126707*y(172)*1.022110952269*T(87)));
T(89) = (-(params(1)*0.062028126707*(y(126)*2.414874053818+y(125)*y(169))));

end
