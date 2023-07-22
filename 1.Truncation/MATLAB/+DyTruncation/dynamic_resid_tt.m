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

assert(length(T) >= 41);

T(1) = y(54)^((-params(6))-1);
T(2) = 1.063205188402*y(66)^(-params(6));
T(3) = y(129)^(-params(6));
T(4) = y(134)^(-params(6));
T(5) = y(142)^(-params(6));
T(6) = y(66)^((-params(6))-1);
T(7) = y(138)^(-params(6));
T(8) = y(146)^(-params(6));
T(9) = y(150)^(-params(6));
T(10) = y(154)^(-params(6));
T(11) = y(158)^(-params(6));
T(12) = y(162)^(-params(6));
T(13) = y(166)^(-params(6));
T(14) = y(170)^(-params(6));
T(15) = 1.105333749474*y(72)^(-params(6));
T(16) = y(72)^((-params(6))-1);
T(17) = 1.021196458903*y(78)^(-params(6));
T(18) = y(78)^((-params(6))-1);
T(19) = 1.042236122064*y(84)^(-params(6));
T(20) = y(84)^((-params(6))-1);
T(21) = 1.058263406894*y(90)^(-params(6));
T(22) = y(90)^((-params(6))-1);
T(23) = 1.018670598592*y(96)^(-params(6));
T(24) = y(96)^((-params(6))-1);
T(25) = 1.029944405894*y(102)^(-params(6));
T(26) = y(102)^((-params(6))-1);
T(27) = 1.039669027202*y(108)^(-params(6));
T(28) = y(108)^((-params(6))-1);
T(29) = 1.015225620707*y(114)^(-params(6));
T(30) = y(114)^((-params(6))-1);
T(31) = 1.022110952269*y(120)^(-params(6));
T(32) = y(120)^((-params(6))-1);
T(33) = (1-params(2))*y(30)^(1/(1-params(2)));
T(34) = ((y(23)-1+params(5))/params(2))^(params(2)/(params(2)-1));
T(35) = y(2)^(params(2)-1);
T(36) = y(2)^params(2);
T(37) = y(2)^(params(2)-2);
T(38) = params(1)*(0.062028126707*y(171)*(y(126)*2.414874053818+y(125)*y(169))+0.000471873293*y(167)*(y(125)*y(165)+y(126)*2.414874053818)+0.000471873293*y(163)*(y(126)*1.459493862568+y(125)*y(161))+0.248112501707*y(159)*(y(126)*1.459493862568+y(125)*y(157))+0.001415625000*y(155)*(y(125)*y(153)+y(126)*1.459493862568)+0.001415625000*y(151)*(y(126)*0.882084236031+y(125)*y(149))+y(147)*0.372168750000*(y(126)*0.882084236031+y(125)*y(145))+y(143)*0.001415625000*(y(125)*y(141)+y(126)*0.882084236031)+y(139)*0.001415625000*(0.533111251380*y(126)+y(125)*y(137))+y(135)*0.248112501707*(0.533111251380*y(126)+y(125)*y(133))+y(130)*0.000471873293*(0.322200074254*y(126)+y(125)*y(128))+1e-18+y(127)*0.062028126707*(0.322200074254*y(126)+y(125)*5.116037217023345)+0.000471873293*y(132)*(0.533111251380*y(126)+y(125)*5.116037217023345));
T(39) = 0.062028126707*y(172)*1.022110952269*T(14)+0.000471873293*y(168)*1.015225620707*T(13)+0.000471873293*y(164)*1.039669027202*T(12)+0.248112501707*y(160)*1.029944405894*T(11)+0.001415625000*y(156)*1.018670598592*T(10)+0.001415625000*y(152)*1.058263406894*T(9)+0.372168750000*y(148)*1.042236122064*T(8)+T(5)*0.001415625000*y(144)*1.021196458903+0.001415625000*y(140)*1.105333749474*T(7)+T(4)*1.063205188402*0.248112501707*y(136)+1e-18+T(3)*1.176900140040*0.000471873293*y(131);
T(40) = params(1)*y(125)*T(39);
T(41) = params(1)*y(124)*(0.001887493173*y(130)+0.992450006827*y(135)+0.005662500000*y(143))+T(38)+T(40);

end
