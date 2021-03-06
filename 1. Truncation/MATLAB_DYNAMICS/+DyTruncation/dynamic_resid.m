function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
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
%   residual
%

if T_flag
    T = DyTruncation.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(91, 1);
lhs = y(46);
rhs = 0.18100000*y(29)+(1+y(28))*y(45)-y(44);
residual(1) = lhs - rhs;
lhs = y(45);
rhs = 1e-10+0.98070000*y(3)+0.01930000*y(4)+0.00000000*y(5)+0.00000000*y(6)+0.00000000*y(7);
residual(2) = lhs - rhs;
lhs = 0.00000000*y(46)^(-params(5));
rhs = 1e-10*params(1)*(1+y(119));
residual(3) = lhs - rhs;
lhs = y(49);
rhs = 0.18100000*y(29)+(1+y(28))*y(48)-y(47);
residual(4) = lhs - rhs;
lhs = y(48);
rhs = 1e-10+0.00480000*y(8)+0.98080000*y(9)+0.01440000*y(10)+0.00000000*y(11)+0.00000000*y(12);
residual(5) = lhs - rhs;
lhs = y(47);
rhs = (steady_state(20));
residual(6) = lhs - rhs;
lhs = y(52);
rhs = 0.18100000*y(29)+(1+y(28))*1e-10-y(50);
residual(7) = lhs - rhs;
lhs = y(51);
rhs = 1e-10;
residual(8) = lhs - rhs;
lhs = y(50);
rhs = (steady_state(23));
residual(9) = lhs - rhs;
lhs = y(55);
rhs = 0.18100000*y(29)+(1+y(28))*1e-10-y(53);
residual(10) = lhs - rhs;
lhs = y(54);
rhs = 1e-10;
residual(11) = lhs - rhs;
lhs = y(53);
rhs = (steady_state(26));
residual(12) = lhs - rhs;
lhs = y(58);
rhs = 0.18100000*y(29)+(1+y(28))*1e-10-y(56);
residual(13) = lhs - rhs;
lhs = y(57);
rhs = 1e-10;
residual(14) = lhs - rhs;
lhs = y(56);
rhs = (steady_state(29));
residual(15) = lhs - rhs;
lhs = y(61);
rhs = y(29)*0.37400000+(1+y(28))*y(60)-y(59);
residual(16) = lhs - rhs;
lhs = y(60);
rhs = 1e-10+0.98070000*y(3)+0.01930000*y(4)+0.00000000*y(5)+0.00000000*y(6)+0.00000000*y(7);
residual(17) = lhs - rhs;
lhs = y(59);
rhs = (steady_state(32));
residual(18) = lhs - rhs;
lhs = y(64);
rhs = y(29)*0.37400000+(1+y(28))*y(63)-y(62);
residual(19) = lhs - rhs;
lhs = y(63);
rhs = 1e-10+0.00480000*y(8)+0.98080000*y(9)+0.01440000*y(10)+0.00000000*y(11)+0.00000000*y(12);
residual(20) = lhs - rhs;
lhs = y(62);
rhs = (steady_state(35));
residual(21) = lhs - rhs;
lhs = y(67);
rhs = y(29)*0.37400000+(1+y(28))*y(66)-y(65);
residual(22) = lhs - rhs;
lhs = y(66);
rhs = 1e-10+0.00000000*y(13)+0.00960000*y(14)+0.98080000*y(15)+0.00960000*y(16)+0.00000000*y(17);
residual(23) = lhs - rhs;
lhs = y(65);
rhs = (steady_state(38));
residual(24) = lhs - rhs;
lhs = y(70);
rhs = y(29)*0.37400000+(1+y(28))*1e-10-y(68);
residual(25) = lhs - rhs;
lhs = y(69);
rhs = 1e-10;
residual(26) = lhs - rhs;
lhs = y(68);
rhs = (steady_state(41));
residual(27) = lhs - rhs;
lhs = y(73);
rhs = y(29)*0.37400000+(1+y(28))*1e-10-y(71);
residual(28) = lhs - rhs;
lhs = y(72);
rhs = 1e-10;
residual(29) = lhs - rhs;
lhs = y(71);
rhs = (steady_state(44));
residual(30) = lhs - rhs;
lhs = y(76);
rhs = y(29)*0.77300000+(1+y(28))*1e-10-y(74);
residual(31) = lhs - rhs;
lhs = y(75);
rhs = 1e-10;
residual(32) = lhs - rhs;
lhs = y(74);
rhs = (steady_state(47));
residual(33) = lhs - rhs;
lhs = y(79);
rhs = y(29)*0.77300000+(1+y(28))*y(78)-y(77);
residual(34) = lhs - rhs;
lhs = y(78);
rhs = 1e-10+0.00480000*y(8)+0.98080000*y(9)+0.01440000*y(10)+0.00000000*y(11)+0.00000000*y(12);
residual(35) = lhs - rhs;
lhs = y(77);
rhs = (steady_state(50));
residual(36) = lhs - rhs;
lhs = y(82);
rhs = y(29)*0.77300000+(1+y(28))*y(81)-y(80);
residual(37) = lhs - rhs;
lhs = y(81);
rhs = 1e-10+0.00000000*y(13)+0.00960000*y(14)+0.98080000*y(15)+0.00960000*y(16)+0.00000000*y(17);
residual(38) = lhs - rhs;
lhs = y(80);
rhs = (steady_state(53));
residual(39) = lhs - rhs;
lhs = y(85);
rhs = y(29)*0.77300000+(1+y(28))*y(84)-y(83);
residual(40) = lhs - rhs;
lhs = y(84);
rhs = 1e-10+0.00000000*y(18)+0.00000000*y(19)+0.01440000*y(20)+0.98080000*y(21)+0.00480000*y(22);
residual(41) = lhs - rhs;
lhs = y(83);
rhs = (steady_state(56));
residual(42) = lhs - rhs;
lhs = y(88);
rhs = y(29)*0.77300000+(1+y(28))*1e-10-y(86);
residual(43) = lhs - rhs;
lhs = y(87);
rhs = 1e-10;
residual(44) = lhs - rhs;
lhs = y(86);
rhs = (steady_state(59));
residual(45) = lhs - rhs;
lhs = y(91);
rhs = y(29)*1.59760000+(1+y(28))*1e-10-y(89);
residual(46) = lhs - rhs;
lhs = y(90);
rhs = 1e-10;
residual(47) = lhs - rhs;
lhs = y(89);
rhs = (steady_state(62));
residual(48) = lhs - rhs;
lhs = y(94);
rhs = y(29)*1.59760000+(1+y(28))*1e-10-y(92);
residual(49) = lhs - rhs;
lhs = y(93);
rhs = 1e-10;
residual(50) = lhs - rhs;
lhs = y(92);
rhs = (steady_state(65));
residual(51) = lhs - rhs;
lhs = y(97);
rhs = y(29)*1.59760000+(1+y(28))*y(96)-y(95);
residual(52) = lhs - rhs;
lhs = y(96);
rhs = 1e-10+0.00000000*y(13)+0.00960000*y(14)+0.98080000*y(15)+0.00960000*y(16)+0.00000000*y(17);
residual(53) = lhs - rhs;
lhs = y(95);
rhs = (steady_state(68));
residual(54) = lhs - rhs;
lhs = y(100);
rhs = y(29)*1.59760000+(1+y(28))*y(99)-y(98);
residual(55) = lhs - rhs;
lhs = y(99);
rhs = 1e-10+0.00000000*y(18)+0.00000000*y(19)+0.01440000*y(20)+0.98080000*y(21)+0.00480000*y(22);
residual(56) = lhs - rhs;
lhs = y(98);
rhs = (steady_state(71));
residual(57) = lhs - rhs;
lhs = y(103);
rhs = y(29)*1.59760000+(1+y(28))*y(102)-y(101);
residual(58) = lhs - rhs;
lhs = y(102);
rhs = 1e-10+0.00000000*y(23)+0.00000000*y(24)+0.00000000*y(25)+0.01930000*y(26)+0.98070000*y(27);
residual(59) = lhs - rhs;
lhs = y(101);
rhs = (steady_state(74));
residual(60) = lhs - rhs;
lhs = y(106);
rhs = y(29)*3.30230000+(1+y(28))*1e-10-y(104);
residual(61) = lhs - rhs;
lhs = y(105);
rhs = 1e-10;
residual(62) = lhs - rhs;
lhs = y(104);
rhs = (steady_state(77));
residual(63) = lhs - rhs;
lhs = y(109);
rhs = y(29)*3.30230000+(1+y(28))*1e-10-y(107);
residual(64) = lhs - rhs;
lhs = y(108);
rhs = 1e-10;
residual(65) = lhs - rhs;
lhs = y(107);
rhs = (steady_state(80));
residual(66) = lhs - rhs;
lhs = y(112);
rhs = y(29)*3.30230000+(1+y(28))*1e-10-y(110);
residual(67) = lhs - rhs;
lhs = y(111);
rhs = 1e-10;
residual(68) = lhs - rhs;
lhs = y(110);
rhs = (steady_state(83));
residual(69) = lhs - rhs;
lhs = y(115);
rhs = y(29)*3.30230000+(1+y(28))*y(114)-y(113);
residual(70) = lhs - rhs;
lhs = y(114);
rhs = 1e-10+0.00000000*y(18)+0.00000000*y(19)+0.01440000*y(20)+0.98080000*y(21)+0.00480000*y(22);
residual(71) = lhs - rhs;
lhs = y(113);
rhs = (steady_state(86));
residual(72) = lhs - rhs;
lhs = y(118);
rhs = y(29)*3.30230000+(1+y(28))*y(117)-y(116);
residual(73) = lhs - rhs;
lhs = y(117);
rhs = 1e-10+0.00000000*y(23)+0.00000000*y(24)+0.00000000*y(25)+0.01930000*y(26)+0.98070000*y(27);
residual(74) = lhs - rhs;
lhs = y(116);
rhs = (steady_state(89));
residual(75) = lhs - rhs;
lhs = y(29);
rhs = T(1)*T(2);
residual(76) = lhs - rhs;
lhs = (y(28)+params(4))/(params(2)*y(30));
rhs = y(31)^(params(2)-1);
residual(77) = lhs - rhs;
lhs = y(30);
rhs = 1+y(33);
residual(78) = lhs - rhs;
lhs = y(33);
rhs = params(6)*y(2)+x(it_, 1);
residual(79) = lhs - rhs;
lhs = y(32);
rhs = y(30)*T(3);
residual(80) = lhs - rhs;
lhs = y(35);
rhs = y(46)*0.06101568+y(49)*0.00120078+0.00000000*y(52)+0.00000000*y(55)+0.00000000*y(58)+y(61)*0.00120078+y(64)*0.24535891+y(67)*0.00360233+0.00000000*y(70)+0.00000000*y(73)+0.00000000*y(76)+y(79)*0.00360233+y(82)*0.36803837+y(85)*0.00360233+0.00000000*y(88)+0.00000000*y(91)+0.00000000*y(94)+y(97)*0.00360233+y(100)*0.24535891+y(103)*0.00120078+0.00000000*y(106)+0.00000000*y(109)+0.00000000*y(112)+y(115)*0.00120078+y(118)*0.06101568;
residual(81) = lhs - rhs;
lhs = y(34);
rhs = y(120)-y(31)*(1-params(4));
residual(82) = lhs - rhs;
lhs = y(31);
rhs = y(44)*0.06101568+y(47)*0.00120078+0.00000000*y(50)+0.00000000*y(53)+0.00000000*y(56)+y(59)*0.00120078+y(62)*0.24535891+y(65)*0.00360233+0.00000000*y(68)+0.00000000*y(71)+0.00000000*y(74)+y(77)*0.00360233+y(80)*0.36803837+y(83)*0.00360233+0.00000000*y(86)+0.00000000*y(89)+0.00000000*y(92)+y(95)*0.00360233+y(98)*0.24535891+y(101)*0.00120078+0.00000000*y(104)+0.00000000*y(107)+0.00000000*y(110)+y(113)*0.00120078+y(116)*0.06101568;
residual(83) = lhs - rhs;
lhs = y(41);
rhs = y(33)*100;
residual(84) = lhs - rhs;
lhs = y(39);
rhs = 100*(y(31)/29.36173135-1);
residual(85) = lhs - rhs;
lhs = y(43);
rhs = 100*(y(35)/(steady_state(8))-1);
residual(86) = lhs - rhs;
lhs = y(37);
rhs = 100*(y(29)/2.16063483-1);
residual(87) = lhs - rhs;
lhs = y(36);
rhs = 100*(y(28)-0.01639256);
residual(88) = lhs - rhs;
lhs = y(38);
rhs = y(41);
residual(89) = lhs - rhs;
lhs = y(40);
rhs = 100*(y(32)/3.37599191-1);
residual(90) = lhs - rhs;
lhs = y(42);
rhs = 100*(y(34)/0.73404328-1);
residual(91) = lhs - rhs;

end
