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
residual = zeros(101, 1);
lhs = y(26);
rhs = 0.322200074254*y(2)+y(1)*5.116037217023345-5.072882123885-y(14);
residual(1) = lhs - rhs;
lhs = y(25);
rhs = 5.116037217023345;
residual(2) = lhs - rhs;
lhs = y(24);
rhs = 5.072882123885;
residual(3) = lhs - rhs;
lhs = y(27);
rhs = 1.095562480837*y(26)^(-params(6));
residual(4) = lhs - rhs;
residual(5) = y(28);
residual(6) = y(29);
lhs = y(32);
rhs = 0.322200074254*y(2)+y(1)*y(31)-10.788809175582-y(14);
residual(7) = lhs - rhs;
lhs = y(31);
rhs = 0.005662500000*y(48)+0.992450006827*y(42)+0.009581417989538045;
residual(8) = lhs - rhs;
lhs = y(30);
rhs = 10.788809175582;
residual(9) = lhs - rhs;
lhs = y(33);
rhs = 1.176900140040*T(1)+y(1)*y(35)*(-1.491436810105)*T(2);
residual(10) = lhs - rhs;
residual(11) = y(34);
lhs = y(35);
rhs = 0.992450006827*y(46)+0.005662500000*y(52);
residual(12) = lhs - rhs;
lhs = y(38);
rhs = y(1)*5.116037217023345+y(2)*0.533111251380-5.076266302097-y(14);
residual(13) = lhs - rhs;
lhs = y(37);
rhs = 5.116037217023345;
residual(14) = lhs - rhs;
lhs = y(36);
rhs = 5.076266302097;
residual(15) = lhs - rhs;
lhs = y(39);
rhs = 1.032942242255*y(38)^(-params(6));
residual(16) = lhs - rhs;
residual(17) = y(40);
residual(18) = y(41);
lhs = y(44);
rhs = y(2)*0.533111251380+y(1)*y(43)-y(42)-y(14);
residual(19) = lhs - rhs;
lhs = y(43);
rhs = 0.005662500000*y(48)+0.992450006827*y(42)+0.009581417989538045;
residual(20) = lhs - rhs;
lhs = 1.063205188402*T(3);
rhs = 0.001887341574+y(1)*params(1)*(T(1)*0.002221390979+T(3)*1.055177996489+0.005782524949*T(4));
residual(21) = lhs - rhs;
lhs = y(45);
rhs = 1.063205188402*T(3)-(y(46)-y(1)*y(47))*(-1.170362773073)*T(5);
residual(22) = lhs - rhs;
lhs = y(45);
rhs = T(17);
residual(23) = lhs - rhs;
lhs = y(47);
rhs = 0.992450006827*y(46)+0.005662500000*y(52);
residual(24) = lhs - rhs;
lhs = y(50);
rhs = y(2)*0.533111251380+y(1)*y(49)-y(48)-y(14);
residual(25) = lhs - rhs;
lhs = y(49);
rhs = 1e-18+0.003775000000*y(54)+0.992450000000*y(60)+0.003775000000*y(66);
residual(26) = lhs - rhs;
lhs = 1.105333749474*T(14);
rhs = y(1)*params(1)*(T(1)*0.002221390979+T(3)*1.055177996489+0.005782524949*T(4))-0.130033653512;
residual(27) = lhs - rhs;
lhs = y(51);
rhs = 1.105333749474*T(14)-(y(52)-y(1)*y(53))*(-1.303759795716)*T(18);
residual(28) = lhs - rhs;
lhs = y(51);
rhs = T(17);
residual(29) = lhs - rhs;
lhs = y(53);
rhs = 0.003775000000*y(58)+0.992450000000*y(64)+0.003775000000*y(70);
residual(30) = lhs - rhs;
lhs = y(56);
rhs = y(2)*0.882084236031+y(1)*y(55)-y(54)-y(14);
residual(31) = lhs - rhs;
lhs = y(55);
rhs = 0.005662500000*y(48)+0.992450006827*y(42)+0.009581417989538045;
residual(32) = lhs - rhs;
lhs = T(4)*1.021196458903;
rhs = 0.040063355768+y(1)*params(1)*(T(14)*0.004172634904+T(13)*1.034367239342+T(11)*0.003845481510);
residual(33) = lhs - rhs;
lhs = y(57);
rhs = T(4)*1.021196458903-(y(58)-y(1)*y(59))*(-1.057406720804)*T(19);
residual(34) = lhs - rhs;
lhs = y(57);
rhs = T(16)+T(6)+y(1)*params(1)*(y(51)*0.003775000000+y(63)*0.992450000000+y(75)*0.003775000000);
residual(35) = lhs - rhs;
lhs = y(59);
rhs = 0.992450006827*y(46)+0.005662500000*y(52);
residual(36) = lhs - rhs;
lhs = y(62);
rhs = y(2)*0.882084236031+y(1)*y(61)-y(60)-y(14);
residual(37) = lhs - rhs;
lhs = y(61);
rhs = 1e-18+0.003775000000*y(54)+0.992450000000*y(60)+0.003775000000*y(66);
residual(38) = lhs - rhs;
lhs = 1.042236122064*T(13);
rhs = y(1)*params(1)*(T(14)*0.004172634904+T(13)*1.034367239342+T(11)*0.003845481510)+0.000176382223;
residual(39) = lhs - rhs;
lhs = y(63);
rhs = 1.042236122064*T(13)-(y(64)-y(1)*y(65))*(-1.118399718654)*T(20);
residual(40) = lhs - rhs;
lhs = y(63);
rhs = T(16)+T(6)+y(1)*params(1)*(y(51)*0.003775000000+y(63)*0.992450000000+y(75)*0.003775000000);
residual(41) = lhs - rhs;
lhs = y(65);
rhs = 0.003775000000*y(58)+0.992450000000*y(64)+0.003775000000*y(70);
residual(42) = lhs - rhs;
lhs = y(68);
rhs = y(2)*0.882084236031+y(1)*y(67)-y(66)-y(14);
residual(43) = lhs - rhs;
lhs = y(67);
rhs = 1e-18+0.005662500000*y(72)+0.992450006827*y(78)+0.001887493173*y(84);
residual(44) = lhs - rhs;
lhs = 1.058263406894*T(12);
rhs = y(1)*params(1)*(T(14)*0.004172634904+T(13)*1.034367239342+T(11)*0.003845481510)-0.087261980631;
residual(45) = lhs - rhs;
lhs = y(69);
rhs = 1.058263406894*T(12)-(y(70)-y(1)*y(71))*(-1.171814211092)*T(21);
residual(46) = lhs - rhs;
lhs = y(69);
rhs = T(16)+T(6)+y(1)*params(1)*(y(51)*0.003775000000+y(63)*0.992450000000+y(75)*0.003775000000);
residual(47) = lhs - rhs;
lhs = y(71);
rhs = 0.005662500000*y(76)+0.992450006827*y(82)+0.001887493173*y(88);
residual(48) = lhs - rhs;
lhs = y(74);
rhs = y(2)*1.459493862568+y(1)*y(73)-y(72)-y(14);
residual(49) = lhs - rhs;
lhs = y(73);
rhs = 1e-18+0.003775000000*y(54)+0.992450000000*y(60)+0.003775000000*y(66);
residual(50) = lhs - rhs;
lhs = 1.018670598592*T(11);
rhs = 0.038258931254+y(1)*params(1)*(T(12)*0.005992416542+T(10)*1.022168332661+T(8)*0.001916231428);
residual(51) = lhs - rhs;
lhs = y(75);
rhs = 1.018670598592*T(11)-(y(76)-y(1)*y(77))*(-1.052456197938)*T(22);
residual(52) = lhs - rhs;
lhs = y(75);
rhs = T(16)+T(6)+y(1)*params(1)*(0.005662500000*y(69)+0.992450006827*y(81)+0.001887493173*y(93));
residual(53) = lhs - rhs;
lhs = y(77);
rhs = 0.003775000000*y(58)+0.992450000000*y(64)+0.003775000000*y(70);
residual(54) = lhs - rhs;
lhs = y(80);
rhs = y(2)*1.459493862568+y(1)*y(79)-y(78)-y(14);
residual(55) = lhs - rhs;
lhs = y(79);
rhs = 1e-18+0.005662500000*y(72)+0.992450006827*y(78)+0.001887493173*y(84);
residual(56) = lhs - rhs;
lhs = 1.029944405894*T(10);
rhs = y(1)*params(1)*(T(12)*0.005992416542+T(10)*1.022168332661+T(8)*0.001916231428)-0.000120269918;
residual(57) = lhs - rhs;
lhs = y(81);
rhs = 1.029944405894*T(10)-(y(82)-y(1)*y(83))*(-1.086909626500)*T(23);
residual(58) = lhs - rhs;
lhs = y(81);
rhs = T(16)+T(6)+y(1)*params(1)*(0.005662500000*y(69)+0.992450006827*y(81)+0.001887493173*y(93));
residual(59) = lhs - rhs;
lhs = y(83);
rhs = 0.005662500000*y(76)+0.992450006827*y(82)+0.001887493173*y(88);
residual(60) = lhs - rhs;
lhs = y(86);
rhs = y(2)*1.459493862568+y(1)*y(85)-y(84)-y(14);
residual(61) = lhs - rhs;
lhs = y(85);
rhs = 1e-18+0.007549972690*y(90)+0.992450027310*y(96);
residual(62) = lhs - rhs;
lhs = 1.039669027202*T(9);
rhs = y(1)*params(1)*(T(12)*0.005992416542+T(10)*1.022168332661+T(8)*0.001916231428)-0.052244522213;
residual(63) = lhs - rhs;
lhs = y(87);
rhs = 1.039669027202*T(9)-(y(88)-y(1)*y(89))*(-1.121502788713)*T(24);
residual(64) = lhs - rhs;
lhs = y(87);
rhs = T(16)+T(6)+y(1)*params(1)*(0.005662500000*y(69)+0.992450006827*y(81)+0.001887493173*y(93));
residual(65) = lhs - rhs;
lhs = y(89);
rhs = 0.007549972690*y(94)+0.992450027310*y(100);
residual(66) = lhs - rhs;
lhs = y(92);
rhs = y(2)*2.414874053818+y(1)*y(91)-y(90)-y(14);
residual(67) = lhs - rhs;
lhs = y(91);
rhs = 1e-18+0.005662500000*y(72)+0.992450006827*y(78)+0.001887493173*y(84);
residual(68) = lhs - rhs;
lhs = 1.015225620707*T(8);
rhs = 0.026809778013+y(1)*params(1)*(T(9)*0.007849472762+T(7)*1.014394042493);
residual(69) = lhs - rhs;
lhs = y(93);
rhs = 1.015225620707*T(8)-(y(94)-y(1)*y(95))*(-1.044154917184)*T(25);
residual(70) = lhs - rhs;
lhs = y(93);
rhs = T(16)+T(6)+y(1)*params(1)*(y(87)*0.007549972690+y(99)*0.992450027310);
residual(71) = lhs - rhs;
lhs = y(95);
rhs = 0.005662500000*y(76)+0.992450006827*y(82)+0.001887493173*y(88);
residual(72) = lhs - rhs;
lhs = y(98);
rhs = y(2)*2.414874053818+y(1)*y(97)-y(96)-y(14);
residual(73) = lhs - rhs;
lhs = y(97);
rhs = 1e-18+0.007549972690*y(90)+0.992450027310*y(96);
residual(74) = lhs - rhs;
lhs = 1.022110952269*T(7);
rhs = y(1)*params(1)*(T(9)*0.007849472762+T(7)*1.014394042493)-0.000206478253;
residual(75) = lhs - rhs;
lhs = y(99);
rhs = 1.022110952269*T(7)-(y(100)-y(1)*y(101))*(-1.066638944053)*T(26);
residual(76) = lhs - rhs;
lhs = y(99);
rhs = T(16)+T(6)+y(1)*params(1)*(y(87)*0.007549972690+y(99)*0.992450027310);
residual(77) = lhs - rhs;
lhs = y(101);
rhs = 0.007549972690*y(94)+0.992450027310*y(100);
residual(78) = lhs - rhs;
lhs = params(3)*y(14)^(params(3)-1);
rhs = 0.062028126707*y(99)+0.000471873293*y(93)+0.000471873293*y(87)+0.248112501707*y(81)+0.001415625000*y(75)+0.001415625000*y(69)+y(63)*0.372168750000+y(57)*0.001415625000+0.001415625000*y(51)+y(45)*0.248112501707+y(39)*0.000471873293+y(33)*0.000471873293+y(27)*0.062028126707;
residual(79) = lhs - rhs;
lhs = y(2);
rhs = T(27)*T(28);
residual(80) = lhs - rhs;
lhs = (y(1)-1+params(5))/(params(2)*y(8));
rhs = T(29);
residual(81) = lhs - rhs;
lhs = y(14)+y(13)+y(12);
rhs = y(11);
residual(82) = lhs - rhs;
lhs = y(4);
rhs = (1-params(2))*y(8)*T(30);
residual(83) = lhs - rhs;
lhs = y(3);
rhs = params(2)*y(8)*T(29);
residual(84) = lhs - rhs;
lhs = y(6);
rhs = T(29)*y(8)*params(2)*(1-params(2));
residual(85) = lhs - rhs;
lhs = y(5);
rhs = y(8)*params(2)*(params(2)-1)*T(31);
residual(86) = lhs - rhs;
lhs = y(8);
rhs = 1+y(7);
residual(87) = lhs - rhs;
lhs = y(7);
rhs = y(7)*params(7)+x(1);
residual(88) = lhs - rhs;
lhs = y(9);
rhs = y(10);
residual(89) = lhs - rhs;
lhs = y(13);
rhs = y(26)*0.062028126707+y(32)*0.000471873293+y(38)*0.000471873293+y(44)*0.248112501707+0.001415625000*y(50)+y(56)*0.001415625000+0.372168750000*y(62)+0.001415625000*y(68)+0.001415625000*y(74)+0.248112501707*y(80)+0.000471873293*y(86)+0.000471873293*y(92)+0.062028126707*y(98);
residual(90) = lhs - rhs;
lhs = y(12);
rhs = y(10)-y(10)*(1-params(5));
residual(91) = lhs - rhs;
lhs = y(9);
rhs = 0.062028126707*y(96)+0.000471873293*y(90)+0.000471873293*y(84)+0.248112501707*y(78)+0.001415625000*y(72)+0.001415625000*y(66)+0.372168750000*y(60)+0.001415625000*y(54)+y(48)*0.001415625000+y(42)*0.248112501707+0.3221476805593601;
residual(92) = lhs - rhs;
lhs = y(17);
rhs = y(7)*100;
residual(93) = lhs - rhs;
lhs = y(19);
rhs = 100*(y(10)/40.590193574133-1);
residual(94) = lhs - rhs;
lhs = y(22);
rhs = 100*(y(13)/2.474873519141-1);
residual(95) = lhs - rhs;
lhs = y(16);
rhs = 100*(y(2)/2.427792565820-1);
residual(96) = lhs - rhs;
lhs = y(15);
rhs = 100*(y(1)-1.008644415200);
residual(97) = lhs - rhs;
lhs = y(23);
rhs = 100*(y(14)/0.303797532989-1);
residual(98) = lhs - rhs;
lhs = y(18);
rhs = y(17);
residual(99) = lhs - rhs;
lhs = y(20);
rhs = 100*(y(11)/3.793425891477-1);
residual(100) = lhs - rhs;
lhs = y(21);
rhs = 100*(y(12)/1.014754839353-1);
residual(101) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end