function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = DyTruncation.static_resid_tt(T, y, x, params);
end
residual = zeros(91, 1);
lhs = y(19);
rhs = 0.181000000000*y(2)+(1+y(1))*y(18)-y(17);
residual(1) = lhs - rhs;
lhs = y(18);
rhs = 1e-10+y(17)*0.980700000000+0.019300000000*y(20);
residual(2) = lhs - rhs;
lhs = y(17);
rhs = (y(17));
residual(3) = lhs - rhs;
lhs = y(22);
rhs = 0.181000000000*y(2)+(1+y(1))*y(21)-y(20);
residual(4) = lhs - rhs;
lhs = y(21);
rhs = 1e-10+0.004800000000*y(32)+0.980800000000*y(35)+0.014400000000*y(38);
residual(5) = lhs - rhs;
lhs = 1.144749942929*T(1);
rhs = (1+y(1))*params(1)*T(2);
residual(6) = lhs - rhs;
lhs = y(25);
rhs = 0.181000000000*y(2)+(1+y(1))*1e-10-y(23);
residual(7) = lhs - rhs;
lhs = y(24);
rhs = 1e-10;
residual(8) = lhs - rhs;
lhs = y(23);
rhs = (y(23));
residual(9) = lhs - rhs;
lhs = y(28);
rhs = 0.181000000000*y(2)+(1+y(1))*1e-10-y(26);
residual(10) = lhs - rhs;
lhs = y(27);
rhs = 1e-10;
residual(11) = lhs - rhs;
lhs = y(26);
rhs = (y(26));
residual(12) = lhs - rhs;
lhs = y(31);
rhs = 0.181000000000*y(2)+(1+y(1))*1e-10-y(29);
residual(13) = lhs - rhs;
lhs = y(30);
rhs = 1e-10;
residual(14) = lhs - rhs;
lhs = y(29);
rhs = (y(29));
residual(15) = lhs - rhs;
lhs = y(34);
rhs = y(2)*0.374000000000+(1+y(1))*y(33)-y(32);
residual(16) = lhs - rhs;
lhs = y(33);
rhs = 1e-10+y(17)*0.980700000000+0.019300000000*y(20);
residual(17) = lhs - rhs;
lhs = y(32);
rhs = (y(32));
residual(18) = lhs - rhs;
lhs = y(37);
rhs = y(2)*0.374000000000+(1+y(1))*y(36)-y(35);
residual(19) = lhs - rhs;
lhs = y(36);
rhs = 1e-10+0.004800000000*y(32)+0.980800000000*y(35)+0.014400000000*y(38);
residual(20) = lhs - rhs;
lhs = 0.643312271893*T(3);
rhs = T(5);
residual(21) = lhs - rhs;
lhs = y(40);
rhs = y(2)*0.374000000000+(1+y(1))*y(39)-y(38);
residual(22) = lhs - rhs;
lhs = y(39);
rhs = 1e-10+0.009600000000*y(50)+0.980800000000*y(53)+0.009600000000*y(56);
residual(23) = lhs - rhs;
lhs = 0.863143214338*T(6);
rhs = T(5);
residual(24) = lhs - rhs;
lhs = y(43);
rhs = (1+y(1))*1e-10+y(2)*0.374000000000-y(41);
residual(25) = lhs - rhs;
lhs = y(42);
rhs = 1e-10;
residual(26) = lhs - rhs;
lhs = y(41);
rhs = (y(41));
residual(27) = lhs - rhs;
lhs = y(46);
rhs = (1+y(1))*1e-10+y(2)*0.374000000000-y(44);
residual(28) = lhs - rhs;
lhs = y(45);
rhs = 1e-10;
residual(29) = lhs - rhs;
lhs = y(44);
rhs = (y(44));
residual(30) = lhs - rhs;
lhs = y(49);
rhs = (1+y(1))*1e-10+y(2)*0.773000000000-y(47);
residual(31) = lhs - rhs;
lhs = y(48);
rhs = 1e-10;
residual(32) = lhs - rhs;
lhs = y(47);
rhs = (y(47));
residual(33) = lhs - rhs;
lhs = y(52);
rhs = y(2)*0.773000000000+(1+y(1))*y(51)-y(50);
residual(34) = lhs - rhs;
lhs = y(51);
rhs = 1e-10+0.004800000000*y(32)+0.980800000000*y(35)+0.014400000000*y(38);
residual(35) = lhs - rhs;
lhs = T(4)*0.690551442596;
rhs = T(9);
residual(36) = lhs - rhs;
lhs = y(55);
rhs = y(2)*0.773000000000+(1+y(1))*y(54)-y(53);
residual(37) = lhs - rhs;
lhs = y(54);
rhs = 1e-10+0.009600000000*y(50)+0.980800000000*y(53)+0.009600000000*y(56);
residual(38) = lhs - rhs;
lhs = T(7)*0.809608796309;
rhs = T(9);
residual(39) = lhs - rhs;
lhs = y(58);
rhs = y(2)*0.773000000000+(1+y(1))*y(57)-y(56);
residual(40) = lhs - rhs;
lhs = y(57);
rhs = 1e-10+0.014400000000*y(68)+0.980800000000*y(71)+0.004800000000*y(74);
residual(41) = lhs - rhs;
lhs = 1.088575212489*T(10);
rhs = T(9);
residual(42) = lhs - rhs;
lhs = y(61);
rhs = (1+y(1))*1e-10+y(2)*0.773000000000-y(59);
residual(43) = lhs - rhs;
lhs = y(60);
rhs = 1e-10;
residual(44) = lhs - rhs;
lhs = y(59);
rhs = (y(59));
residual(45) = lhs - rhs;
lhs = y(64);
rhs = (1+y(1))*1e-10+y(2)*1.597600000000-y(62);
residual(46) = lhs - rhs;
lhs = y(63);
rhs = 1e-10;
residual(47) = lhs - rhs;
lhs = y(62);
rhs = (y(62));
residual(48) = lhs - rhs;
lhs = y(67);
rhs = (1+y(1))*1e-10+y(2)*1.597600000000-y(65);
residual(49) = lhs - rhs;
lhs = y(66);
rhs = 1e-10;
residual(50) = lhs - rhs;
lhs = y(65);
rhs = (y(65));
residual(51) = lhs - rhs;
lhs = y(70);
rhs = y(2)*1.597600000000+(1+y(1))*y(69)-y(68);
residual(52) = lhs - rhs;
lhs = y(69);
rhs = 1e-10+0.009600000000*y(50)+0.980800000000*y(53)+0.009600000000*y(56);
residual(53) = lhs - rhs;
lhs = T(8)*0.974957774233;
rhs = T(13);
residual(54) = lhs - rhs;
lhs = y(73);
rhs = y(2)*1.597600000000+(1+y(1))*y(72)-y(71);
residual(55) = lhs - rhs;
lhs = y(72);
rhs = 1e-10+0.014400000000*y(68)+0.980800000000*y(71)+0.004800000000*y(74);
residual(56) = lhs - rhs;
lhs = T(11)*1.181104990036;
rhs = T(13);
residual(57) = lhs - rhs;
lhs = y(76);
rhs = y(2)*1.597600000000+(1+y(1))*y(75)-y(74);
residual(58) = lhs - rhs;
lhs = y(75);
rhs = 1e-10+0.019300000000*y(86)+0.980700000000*y(89);
residual(59) = lhs - rhs;
lhs = 1.599905549031*T(14);
rhs = T(13);
residual(60) = lhs - rhs;
lhs = y(79);
rhs = (1+y(1))*1e-10+y(2)*3.302300000000-y(77);
residual(61) = lhs - rhs;
lhs = y(78);
rhs = 1e-10;
residual(62) = lhs - rhs;
lhs = y(77);
rhs = (y(77));
residual(63) = lhs - rhs;
lhs = y(82);
rhs = (1+y(1))*1e-10+y(2)*3.302300000000-y(80);
residual(64) = lhs - rhs;
lhs = y(81);
rhs = 1e-10;
residual(65) = lhs - rhs;
lhs = y(80);
rhs = (y(80));
residual(66) = lhs - rhs;
lhs = y(85);
rhs = (1+y(1))*1e-10+y(2)*3.302300000000-y(83);
residual(67) = lhs - rhs;
lhs = y(84);
rhs = 1e-10;
residual(68) = lhs - rhs;
lhs = y(83);
rhs = (y(83));
residual(69) = lhs - rhs;
lhs = y(88);
rhs = y(2)*3.302300000000+(1+y(1))*y(87)-y(86);
residual(70) = lhs - rhs;
lhs = y(87);
rhs = 1e-10+0.014400000000*y(68)+0.980800000000*y(71)+0.004800000000*y(74);
residual(71) = lhs - rhs;
lhs = T(12)*1.502675855760;
rhs = (1+y(1))*params(1)*(T(14)*0.030878177096+1.820260607621*T(15));
residual(72) = lhs - rhs;
lhs = y(91);
rhs = y(2)*3.302300000000+(1+y(1))*y(90)-y(89);
residual(73) = lhs - rhs;
lhs = y(90);
rhs = 1e-10+0.019300000000*y(86)+0.980700000000*y(89);
residual(74) = lhs - rhs;
lhs = T(15)*1.856083009708;
rhs = (1+y(1))*params(1)*(T(14)*0.030878177096+1.820260607621*T(15));
residual(75) = lhs - rhs;
lhs = y(2);
rhs = T(16)*T(17);
residual(76) = lhs - rhs;
lhs = y(1);
rhs = params(2)*y(3)*T(18)-params(4);
residual(77) = lhs - rhs;
lhs = y(3);
rhs = 1+y(6);
residual(78) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(6)+x(1);
residual(79) = lhs - rhs;
lhs = y(5);
rhs = y(3)*T(19);
residual(80) = lhs - rhs;
lhs = y(8);
rhs = y(19)*0.061015683733+y(22)*0.001200777706+y(25)*0.000000000000+y(28)*0.000000000000+y(31)*0.000000000000+y(34)*0.001200777706+y(37)*0.245358911212+y(40)*0.003602333117+y(43)*0.000000000000+y(46)*0.000000000000+y(49)*0.000000000000+y(52)*0.003602333117+y(55)*0.368038366818+y(58)*0.003602333117+y(61)*0.000000000000+y(64)*0.000000000000+y(67)*0.000000000000+y(70)*0.003602333117+y(73)*0.245358911212+y(76)*0.001200777706+y(79)*0.000000000000+y(82)*0.000000000000+y(85)*0.000000000000+y(88)*0.001200777706+y(91)*0.061015683733;
residual(81) = lhs - rhs;
lhs = y(7);
rhs = y(4)-y(4)*(1-params(4));
residual(82) = lhs - rhs;
lhs = y(4);
rhs = y(17)*0.061015683733+y(20)*0.001200777706+y(23)*0.000000000000+y(26)*0.000000000000+y(29)*0.000000000000+y(32)*0.001200777706+y(35)*0.245358911212+y(38)*0.003602333117+y(41)*0.000000000000+y(44)*0.000000000000+y(47)*0.000000000000+y(50)*0.003602333117+y(53)*0.368038366818+y(56)*0.003602333117+y(59)*0.000000000000+y(62)*0.000000000000+y(65)*0.000000000000+y(68)*0.003602333117+y(71)*0.245358911212+y(74)*0.001200777706+y(77)*0.000000000000+y(80)*0.000000000000+y(83)*0.000000000000+y(86)*0.001200777706+y(89)*0.061015683733;
residual(83) = lhs - rhs;
lhs = y(14);
rhs = y(6)*100;
residual(84) = lhs - rhs;
lhs = y(12);
rhs = 100*(y(4)/29.072362347015-1);
residual(85) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(8)/(y(8))-1);
residual(86) = lhs - rhs;
lhs = y(10);
rhs = 100*(y(2)/2.152944776511-1);
residual(87) = lhs - rhs;
lhs = y(9);
rhs = 100*(y(1)-0.016655763035);
residual(88) = lhs - rhs;
lhs = y(11);
rhs = y(14);
residual(89) = lhs - rhs;
lhs = y(13);
rhs = 100*(y(5)/3.363976213299-1);
residual(90) = lhs - rhs;
lhs = y(15);
rhs = 100*(y(7)/0.726809058675-1);
residual(91) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
