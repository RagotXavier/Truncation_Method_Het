function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = DyTruncation.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(91, 121);
g1(1,16)=(-y(33));
g1(1,17)=(-0.181000000000);
g1(1,32)=1;
g1(1,33)=(-(1+y(16)));
g1(1,34)=1;
g1(2,3)=(-0.980700000000);
g1(2,33)=1;
g1(2,4)=(-0.019300000000);
g1(3,32)=1;
g1(4,16)=(-y(36));
g1(4,17)=(-0.181000000000);
g1(4,35)=1;
g1(4,36)=(-(1+y(16)));
g1(4,37)=1;
g1(5,36)=1;
g1(5,5)=(-0.004800000000);
g1(5,6)=(-0.980800000000);
g1(5,7)=(-0.014400000000);
g1(6,107)=(-(params(1)*T(1)));
g1(6,108)=(-(params(1)*(1+y(107))*0.980700000000*getPowerDeriv(y(108),(-params(5)),1)));
g1(6,37)=1.144749942929*getPowerDeriv(y(37),(-params(5)),1);
g1(6,110)=(-(params(1)*(1+y(107))*0.019300000000*getPowerDeriv(y(110),(-params(5)),1)));
g1(7,16)=(-1e-10);
g1(7,17)=(-0.181000000000);
g1(7,38)=1;
g1(7,40)=1;
g1(8,39)=1;
g1(9,38)=1;
g1(10,16)=(-1e-10);
g1(10,17)=(-0.181000000000);
g1(10,41)=1;
g1(10,43)=1;
g1(11,42)=1;
g1(12,41)=1;
g1(13,16)=(-1e-10);
g1(13,17)=(-0.181000000000);
g1(13,44)=1;
g1(13,46)=1;
g1(14,45)=1;
g1(15,44)=1;
g1(16,16)=(-y(48));
g1(16,17)=(-0.374000000000);
g1(16,47)=1;
g1(16,48)=(-(1+y(16)));
g1(16,49)=1;
g1(17,3)=(-0.980700000000);
g1(17,4)=(-0.019300000000);
g1(17,48)=1;
g1(18,47)=1;
g1(19,16)=(-y(51));
g1(19,17)=(-0.374000000000);
g1(19,50)=1;
g1(19,51)=(-(1+y(16)));
g1(19,52)=1;
g1(20,5)=(-0.004800000000);
g1(20,6)=(-0.980800000000);
g1(20,51)=1;
g1(20,7)=(-0.014400000000);
g1(21,107)=(-(params(1)*T(2)));
g1(21,109)=T(14);
g1(21,52)=0.643312271893*getPowerDeriv(y(52),(-params(5)),1);
g1(21,111)=T(15);
g1(21,113)=T(17);
g1(22,16)=(-y(54));
g1(22,17)=(-0.374000000000);
g1(22,53)=1;
g1(22,54)=(-(1+y(16)));
g1(22,55)=1;
g1(23,54)=1;
g1(23,8)=(-0.009600000000);
g1(23,9)=(-0.980800000000);
g1(23,10)=(-0.009600000000);
g1(24,107)=(-(params(1)*T(2)));
g1(24,109)=T(14);
g1(24,111)=T(15);
g1(24,55)=0.863143214338*getPowerDeriv(y(55),(-params(5)),1);
g1(24,113)=T(17);
g1(25,16)=(-1e-10);
g1(25,17)=(-0.374000000000);
g1(25,56)=1;
g1(25,58)=1;
g1(26,57)=1;
g1(27,56)=1;
g1(28,16)=(-1e-10);
g1(28,17)=(-0.374000000000);
g1(28,59)=1;
g1(28,61)=1;
g1(29,60)=1;
g1(30,59)=1;
g1(31,16)=(-1e-10);
g1(31,17)=(-0.773000000000);
g1(31,62)=1;
g1(31,64)=1;
g1(32,63)=1;
g1(33,62)=1;
g1(34,16)=(-y(66));
g1(34,17)=(-0.773000000000);
g1(34,65)=1;
g1(34,66)=(-(1+y(16)));
g1(34,67)=1;
g1(35,5)=(-0.004800000000);
g1(35,6)=(-0.980800000000);
g1(35,7)=(-0.014400000000);
g1(35,66)=1;
g1(36,107)=(-(params(1)*T(4)));
g1(36,112)=T(16);
g1(36,67)=0.690551442596*getPowerDeriv(y(67),(-params(5)),1);
g1(36,114)=T(18);
g1(36,116)=T(20);
g1(37,16)=(-y(69));
g1(37,17)=(-0.773000000000);
g1(37,68)=1;
g1(37,69)=(-(1+y(16)));
g1(37,70)=1;
g1(38,8)=(-0.009600000000);
g1(38,9)=(-0.980800000000);
g1(38,69)=1;
g1(38,10)=(-0.009600000000);
g1(39,107)=(-(params(1)*T(4)));
g1(39,112)=T(16);
g1(39,70)=0.809608796309*getPowerDeriv(y(70),(-params(5)),1);
g1(39,114)=T(18);
g1(39,116)=T(20);
g1(40,16)=(-y(72));
g1(40,17)=(-0.773000000000);
g1(40,71)=1;
g1(40,72)=(-(1+y(16)));
g1(40,73)=1;
g1(41,72)=1;
g1(41,11)=(-0.014400000000);
g1(41,12)=(-0.980800000000);
g1(41,13)=(-0.004800000000);
g1(42,107)=(-(params(1)*T(4)));
g1(42,112)=T(16);
g1(42,114)=T(18);
g1(42,73)=1.088575212489*getPowerDeriv(y(73),(-params(5)),1);
g1(42,116)=T(20);
g1(43,16)=(-1e-10);
g1(43,17)=(-0.773000000000);
g1(43,74)=1;
g1(43,76)=1;
g1(44,75)=1;
g1(45,74)=1;
g1(46,16)=(-1e-10);
g1(46,17)=(-1.597600000000);
g1(46,77)=1;
g1(46,79)=1;
g1(47,78)=1;
g1(48,77)=1;
g1(49,16)=(-1e-10);
g1(49,17)=(-1.597600000000);
g1(49,80)=1;
g1(49,82)=1;
g1(50,81)=1;
g1(51,80)=1;
g1(52,16)=(-y(84));
g1(52,17)=(-1.597600000000);
g1(52,83)=1;
g1(52,84)=(-(1+y(16)));
g1(52,85)=1;
g1(53,8)=(-0.009600000000);
g1(53,9)=(-0.980800000000);
g1(53,10)=(-0.009600000000);
g1(53,84)=1;
g1(54,107)=(-(params(1)*T(6)));
g1(54,115)=T(19);
g1(54,85)=0.974957774233*getPowerDeriv(y(85),(-params(5)),1);
g1(54,117)=T(21);
g1(54,119)=T(23);
g1(55,16)=(-y(87));
g1(55,17)=(-1.597600000000);
g1(55,86)=1;
g1(55,87)=(-(1+y(16)));
g1(55,88)=1;
g1(56,11)=(-0.014400000000);
g1(56,12)=(-0.980800000000);
g1(56,87)=1;
g1(56,13)=(-0.004800000000);
g1(57,107)=(-(params(1)*T(6)));
g1(57,115)=T(19);
g1(57,88)=1.181104990036*getPowerDeriv(y(88),(-params(5)),1);
g1(57,117)=T(21);
g1(57,119)=T(23);
g1(58,16)=(-y(90));
g1(58,17)=(-1.597600000000);
g1(58,89)=1;
g1(58,90)=(-(1+y(16)));
g1(58,91)=1;
g1(59,90)=1;
g1(59,14)=(-0.019300000000);
g1(59,15)=(-0.980700000000);
g1(60,107)=(-(params(1)*T(6)));
g1(60,115)=T(19);
g1(60,117)=T(21);
g1(60,91)=1.599905549031*getPowerDeriv(y(91),(-params(5)),1);
g1(60,119)=T(23);
g1(61,16)=(-1e-10);
g1(61,17)=(-3.302300000000);
g1(61,92)=1;
g1(61,94)=1;
g1(62,93)=1;
g1(63,92)=1;
g1(64,16)=(-1e-10);
g1(64,17)=(-3.302300000000);
g1(64,95)=1;
g1(64,97)=1;
g1(65,96)=1;
g1(66,95)=1;
g1(67,16)=(-1e-10);
g1(67,17)=(-3.302300000000);
g1(67,98)=1;
g1(67,100)=1;
g1(68,99)=1;
g1(69,98)=1;
g1(70,16)=(-y(102));
g1(70,17)=(-3.302300000000);
g1(70,101)=1;
g1(70,102)=(-(1+y(16)));
g1(70,103)=1;
g1(71,11)=(-0.014400000000);
g1(71,12)=(-0.980800000000);
g1(71,13)=(-0.004800000000);
g1(71,102)=1;
g1(72,107)=(-(params(1)*T(8)));
g1(72,118)=T(22);
g1(72,103)=1.502675855760*getPowerDeriv(y(103),(-params(5)),1);
g1(72,120)=T(24);
g1(73,16)=(-y(105));
g1(73,17)=(-3.302300000000);
g1(73,104)=1;
g1(73,105)=(-(1+y(16)));
g1(73,106)=1;
g1(74,14)=(-0.019300000000);
g1(74,15)=(-0.980700000000);
g1(74,105)=1;
g1(75,107)=(-(params(1)*T(8)));
g1(75,118)=T(22);
g1(75,106)=1.856083009708*getPowerDeriv(y(106),(-params(5)),1);
g1(75,120)=T(24);
g1(76,16)=(-(T(10)*1/params(2)*getPowerDeriv((y(16)+params(4))/params(2),params(2)/(params(2)-1),1)));
g1(76,17)=1;
g1(76,18)=(-(T(11)*(1-params(2))*getPowerDeriv(y(18),1/(1-params(2)),1)));
g1(77,16)=1;
g1(77,18)=(-(params(2)*T(12)));
g1(77,1)=(-(params(2)*y(18)*getPowerDeriv(y(1),params(2)-1,1)));
g1(78,18)=1;
g1(78,21)=(-1);
g1(79,2)=(-params(6));
g1(79,21)=1;
g1(79,121)=(-1);
g1(80,18)=(-T(13));
g1(80,1)=(-(y(18)*getPowerDeriv(y(1),params(2),1)));
g1(80,20)=1;
g1(81,23)=1;
g1(81,34)=(-0.061015683733);
g1(81,37)=(-0.001200777706);
g1(81,40)=(-0.000000000000);
g1(81,43)=(-0.000000000000);
g1(81,46)=(-0.000000000000);
g1(81,49)=(-0.001200777706);
g1(81,52)=(-0.245358911212);
g1(81,55)=(-0.003602333117);
g1(81,58)=(-0.000000000000);
g1(81,61)=(-0.000000000000);
g1(81,64)=(-0.000000000000);
g1(81,67)=(-0.003602333117);
g1(81,70)=(-0.368038366818);
g1(81,73)=(-0.003602333117);
g1(81,76)=(-0.000000000000);
g1(81,79)=(-0.000000000000);
g1(81,82)=(-0.000000000000);
g1(81,85)=(-0.003602333117);
g1(81,88)=(-0.245358911212);
g1(81,91)=(-0.001200777706);
g1(81,94)=(-0.000000000000);
g1(81,97)=(-0.000000000000);
g1(81,100)=(-0.000000000000);
g1(81,103)=(-0.001200777706);
g1(81,106)=(-0.061015683733);
g1(82,1)=1-params(4);
g1(82,19)=(-1);
g1(82,22)=1;
g1(83,19)=1;
g1(83,32)=(-0.061015683733);
g1(83,35)=(-0.001200777706);
g1(83,38)=(-0.000000000000);
g1(83,41)=(-0.000000000000);
g1(83,44)=(-0.000000000000);
g1(83,47)=(-0.001200777706);
g1(83,50)=(-0.245358911212);
g1(83,53)=(-0.003602333117);
g1(83,56)=(-0.000000000000);
g1(83,59)=(-0.000000000000);
g1(83,62)=(-0.000000000000);
g1(83,65)=(-0.003602333117);
g1(83,68)=(-0.368038366818);
g1(83,71)=(-0.003602333117);
g1(83,74)=(-0.000000000000);
g1(83,77)=(-0.000000000000);
g1(83,80)=(-0.000000000000);
g1(83,83)=(-0.003602333117);
g1(83,86)=(-0.245358911212);
g1(83,89)=(-0.001200777706);
g1(83,92)=(-0.000000000000);
g1(83,95)=(-0.000000000000);
g1(83,98)=(-0.000000000000);
g1(83,101)=(-0.001200777706);
g1(83,104)=(-0.061015683733);
g1(84,21)=(-100);
g1(84,29)=1;
g1(85,19)=(-3.439692956711771);
g1(85,27)=1;
g1(86,23)=(-(100*1/(steady_state(8))));
g1(86,31)=1;
g1(87,17)=(-46.44800976365827);
g1(87,25)=1;
g1(88,16)=(-100);
g1(88,24)=1;
g1(89,26)=1;
g1(89,29)=(-1);
g1(90,20)=(-29.72672624873632);
g1(90,28)=1;
g1(91,22)=(-137.5877182685418);
g1(91,30)=1;

end
