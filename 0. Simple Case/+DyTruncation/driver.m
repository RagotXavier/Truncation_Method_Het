%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'DyTruncation';
M_.dynare_version = '4.6.2';
oo_.dynare_version = '4.6.2';
options_.dynare_version = '4.6.2';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('DyTruncation.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps'};
M_.exo_names_tex(1) = {'eps'};
M_.exo_names_long(1) = {'eps'};
M_.endo_names = cell(91,1);
M_.endo_names_tex = cell(91,1);
M_.endo_names_long = cell(91,1);
M_.endo_names(1) = {'r'};
M_.endo_names_tex(1) = {'r'};
M_.endo_names_long(1) = {'r'};
M_.endo_names(2) = {'w'};
M_.endo_names_tex(2) = {'w'};
M_.endo_names_long(2) = {'w'};
M_.endo_names(3) = {'Z'};
M_.endo_names_tex(3) = {'Z'};
M_.endo_names_long(3) = {'Z'};
M_.endo_names(4) = {'K'};
M_.endo_names_tex(4) = {'K'};
M_.endo_names_long(4) = {'K'};
M_.endo_names(5) = {'GDP'};
M_.endo_names_tex(5) = {'GDP'};
M_.endo_names_long(5) = {'GDP'};
M_.endo_names(6) = {'u'};
M_.endo_names_tex(6) = {'u'};
M_.endo_names_long(6) = {'u'};
M_.endo_names(7) = {'I'};
M_.endo_names_tex(7) = {'I'};
M_.endo_names_long(7) = {'I'};
M_.endo_names(8) = {'C'};
M_.endo_names_tex(8) = {'C'};
M_.endo_names_long(8) = {'C'};
M_.endo_names(9) = {'rt'};
M_.endo_names_tex(9) = {'rt'};
M_.endo_names_long(9) = {'rt'};
M_.endo_names(10) = {'wt'};
M_.endo_names_tex(10) = {'wt'};
M_.endo_names_long(10) = {'wt'};
M_.endo_names(11) = {'Zt'};
M_.endo_names_tex(11) = {'Zt'};
M_.endo_names_long(11) = {'Zt'};
M_.endo_names(12) = {'Kt'};
M_.endo_names_tex(12) = {'Kt'};
M_.endo_names_long(12) = {'Kt'};
M_.endo_names(13) = {'GDPt'};
M_.endo_names_tex(13) = {'GDPt'};
M_.endo_names_long(13) = {'GDPt'};
M_.endo_names(14) = {'ut'};
M_.endo_names_tex(14) = {'ut'};
M_.endo_names_long(14) = {'ut'};
M_.endo_names(15) = {'It'};
M_.endo_names_tex(15) = {'It'};
M_.endo_names_long(15) = {'It'};
M_.endo_names(16) = {'Ct'};
M_.endo_names_tex(16) = {'Ct'};
M_.endo_names_long(16) = {'Ct'};
M_.endo_names(17) = {'a1'};
M_.endo_names_tex(17) = {'a1'};
M_.endo_names_long(17) = {'a1'};
M_.endo_names(18) = {'at1'};
M_.endo_names_tex(18) = {'at1'};
M_.endo_names_long(18) = {'at1'};
M_.endo_names(19) = {'c1'};
M_.endo_names_tex(19) = {'c1'};
M_.endo_names_long(19) = {'c1'};
M_.endo_names(20) = {'a2'};
M_.endo_names_tex(20) = {'a2'};
M_.endo_names_long(20) = {'a2'};
M_.endo_names(21) = {'at2'};
M_.endo_names_tex(21) = {'at2'};
M_.endo_names_long(21) = {'at2'};
M_.endo_names(22) = {'c2'};
M_.endo_names_tex(22) = {'c2'};
M_.endo_names_long(22) = {'c2'};
M_.endo_names(23) = {'a3'};
M_.endo_names_tex(23) = {'a3'};
M_.endo_names_long(23) = {'a3'};
M_.endo_names(24) = {'at3'};
M_.endo_names_tex(24) = {'at3'};
M_.endo_names_long(24) = {'at3'};
M_.endo_names(25) = {'c3'};
M_.endo_names_tex(25) = {'c3'};
M_.endo_names_long(25) = {'c3'};
M_.endo_names(26) = {'a4'};
M_.endo_names_tex(26) = {'a4'};
M_.endo_names_long(26) = {'a4'};
M_.endo_names(27) = {'at4'};
M_.endo_names_tex(27) = {'at4'};
M_.endo_names_long(27) = {'at4'};
M_.endo_names(28) = {'c4'};
M_.endo_names_tex(28) = {'c4'};
M_.endo_names_long(28) = {'c4'};
M_.endo_names(29) = {'a5'};
M_.endo_names_tex(29) = {'a5'};
M_.endo_names_long(29) = {'a5'};
M_.endo_names(30) = {'at5'};
M_.endo_names_tex(30) = {'at5'};
M_.endo_names_long(30) = {'at5'};
M_.endo_names(31) = {'c5'};
M_.endo_names_tex(31) = {'c5'};
M_.endo_names_long(31) = {'c5'};
M_.endo_names(32) = {'a6'};
M_.endo_names_tex(32) = {'a6'};
M_.endo_names_long(32) = {'a6'};
M_.endo_names(33) = {'at6'};
M_.endo_names_tex(33) = {'at6'};
M_.endo_names_long(33) = {'at6'};
M_.endo_names(34) = {'c6'};
M_.endo_names_tex(34) = {'c6'};
M_.endo_names_long(34) = {'c6'};
M_.endo_names(35) = {'a7'};
M_.endo_names_tex(35) = {'a7'};
M_.endo_names_long(35) = {'a7'};
M_.endo_names(36) = {'at7'};
M_.endo_names_tex(36) = {'at7'};
M_.endo_names_long(36) = {'at7'};
M_.endo_names(37) = {'c7'};
M_.endo_names_tex(37) = {'c7'};
M_.endo_names_long(37) = {'c7'};
M_.endo_names(38) = {'a8'};
M_.endo_names_tex(38) = {'a8'};
M_.endo_names_long(38) = {'a8'};
M_.endo_names(39) = {'at8'};
M_.endo_names_tex(39) = {'at8'};
M_.endo_names_long(39) = {'at8'};
M_.endo_names(40) = {'c8'};
M_.endo_names_tex(40) = {'c8'};
M_.endo_names_long(40) = {'c8'};
M_.endo_names(41) = {'a9'};
M_.endo_names_tex(41) = {'a9'};
M_.endo_names_long(41) = {'a9'};
M_.endo_names(42) = {'at9'};
M_.endo_names_tex(42) = {'at9'};
M_.endo_names_long(42) = {'at9'};
M_.endo_names(43) = {'c9'};
M_.endo_names_tex(43) = {'c9'};
M_.endo_names_long(43) = {'c9'};
M_.endo_names(44) = {'a10'};
M_.endo_names_tex(44) = {'a10'};
M_.endo_names_long(44) = {'a10'};
M_.endo_names(45) = {'at10'};
M_.endo_names_tex(45) = {'at10'};
M_.endo_names_long(45) = {'at10'};
M_.endo_names(46) = {'c10'};
M_.endo_names_tex(46) = {'c10'};
M_.endo_names_long(46) = {'c10'};
M_.endo_names(47) = {'a11'};
M_.endo_names_tex(47) = {'a11'};
M_.endo_names_long(47) = {'a11'};
M_.endo_names(48) = {'at11'};
M_.endo_names_tex(48) = {'at11'};
M_.endo_names_long(48) = {'at11'};
M_.endo_names(49) = {'c11'};
M_.endo_names_tex(49) = {'c11'};
M_.endo_names_long(49) = {'c11'};
M_.endo_names(50) = {'a12'};
M_.endo_names_tex(50) = {'a12'};
M_.endo_names_long(50) = {'a12'};
M_.endo_names(51) = {'at12'};
M_.endo_names_tex(51) = {'at12'};
M_.endo_names_long(51) = {'at12'};
M_.endo_names(52) = {'c12'};
M_.endo_names_tex(52) = {'c12'};
M_.endo_names_long(52) = {'c12'};
M_.endo_names(53) = {'a13'};
M_.endo_names_tex(53) = {'a13'};
M_.endo_names_long(53) = {'a13'};
M_.endo_names(54) = {'at13'};
M_.endo_names_tex(54) = {'at13'};
M_.endo_names_long(54) = {'at13'};
M_.endo_names(55) = {'c13'};
M_.endo_names_tex(55) = {'c13'};
M_.endo_names_long(55) = {'c13'};
M_.endo_names(56) = {'a14'};
M_.endo_names_tex(56) = {'a14'};
M_.endo_names_long(56) = {'a14'};
M_.endo_names(57) = {'at14'};
M_.endo_names_tex(57) = {'at14'};
M_.endo_names_long(57) = {'at14'};
M_.endo_names(58) = {'c14'};
M_.endo_names_tex(58) = {'c14'};
M_.endo_names_long(58) = {'c14'};
M_.endo_names(59) = {'a15'};
M_.endo_names_tex(59) = {'a15'};
M_.endo_names_long(59) = {'a15'};
M_.endo_names(60) = {'at15'};
M_.endo_names_tex(60) = {'at15'};
M_.endo_names_long(60) = {'at15'};
M_.endo_names(61) = {'c15'};
M_.endo_names_tex(61) = {'c15'};
M_.endo_names_long(61) = {'c15'};
M_.endo_names(62) = {'a16'};
M_.endo_names_tex(62) = {'a16'};
M_.endo_names_long(62) = {'a16'};
M_.endo_names(63) = {'at16'};
M_.endo_names_tex(63) = {'at16'};
M_.endo_names_long(63) = {'at16'};
M_.endo_names(64) = {'c16'};
M_.endo_names_tex(64) = {'c16'};
M_.endo_names_long(64) = {'c16'};
M_.endo_names(65) = {'a17'};
M_.endo_names_tex(65) = {'a17'};
M_.endo_names_long(65) = {'a17'};
M_.endo_names(66) = {'at17'};
M_.endo_names_tex(66) = {'at17'};
M_.endo_names_long(66) = {'at17'};
M_.endo_names(67) = {'c17'};
M_.endo_names_tex(67) = {'c17'};
M_.endo_names_long(67) = {'c17'};
M_.endo_names(68) = {'a18'};
M_.endo_names_tex(68) = {'a18'};
M_.endo_names_long(68) = {'a18'};
M_.endo_names(69) = {'at18'};
M_.endo_names_tex(69) = {'at18'};
M_.endo_names_long(69) = {'at18'};
M_.endo_names(70) = {'c18'};
M_.endo_names_tex(70) = {'c18'};
M_.endo_names_long(70) = {'c18'};
M_.endo_names(71) = {'a19'};
M_.endo_names_tex(71) = {'a19'};
M_.endo_names_long(71) = {'a19'};
M_.endo_names(72) = {'at19'};
M_.endo_names_tex(72) = {'at19'};
M_.endo_names_long(72) = {'at19'};
M_.endo_names(73) = {'c19'};
M_.endo_names_tex(73) = {'c19'};
M_.endo_names_long(73) = {'c19'};
M_.endo_names(74) = {'a20'};
M_.endo_names_tex(74) = {'a20'};
M_.endo_names_long(74) = {'a20'};
M_.endo_names(75) = {'at20'};
M_.endo_names_tex(75) = {'at20'};
M_.endo_names_long(75) = {'at20'};
M_.endo_names(76) = {'c20'};
M_.endo_names_tex(76) = {'c20'};
M_.endo_names_long(76) = {'c20'};
M_.endo_names(77) = {'a21'};
M_.endo_names_tex(77) = {'a21'};
M_.endo_names_long(77) = {'a21'};
M_.endo_names(78) = {'at21'};
M_.endo_names_tex(78) = {'at21'};
M_.endo_names_long(78) = {'at21'};
M_.endo_names(79) = {'c21'};
M_.endo_names_tex(79) = {'c21'};
M_.endo_names_long(79) = {'c21'};
M_.endo_names(80) = {'a22'};
M_.endo_names_tex(80) = {'a22'};
M_.endo_names_long(80) = {'a22'};
M_.endo_names(81) = {'at22'};
M_.endo_names_tex(81) = {'at22'};
M_.endo_names_long(81) = {'at22'};
M_.endo_names(82) = {'c22'};
M_.endo_names_tex(82) = {'c22'};
M_.endo_names_long(82) = {'c22'};
M_.endo_names(83) = {'a23'};
M_.endo_names_tex(83) = {'a23'};
M_.endo_names_long(83) = {'a23'};
M_.endo_names(84) = {'at23'};
M_.endo_names_tex(84) = {'at23'};
M_.endo_names_long(84) = {'at23'};
M_.endo_names(85) = {'c23'};
M_.endo_names_tex(85) = {'c23'};
M_.endo_names_long(85) = {'c23'};
M_.endo_names(86) = {'a24'};
M_.endo_names_tex(86) = {'a24'};
M_.endo_names_long(86) = {'a24'};
M_.endo_names(87) = {'at24'};
M_.endo_names_tex(87) = {'at24'};
M_.endo_names_long(87) = {'at24'};
M_.endo_names(88) = {'c24'};
M_.endo_names_tex(88) = {'c24'};
M_.endo_names_long(88) = {'c24'};
M_.endo_names(89) = {'a25'};
M_.endo_names_tex(89) = {'a25'};
M_.endo_names_long(89) = {'a25'};
M_.endo_names(90) = {'at25'};
M_.endo_names_tex(90) = {'at25'};
M_.endo_names_long(90) = {'at25'};
M_.endo_names(91) = {'c25'};
M_.endo_names_tex(91) = {'c25'};
M_.endo_names_long(91) = {'c25'};
M_.endo_partitions = struct();
M_.param_names = cell(6,1);
M_.param_names_tex = cell(6,1);
M_.param_names_long = cell(6,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'alpha'};
M_.param_names_tex(2) = {'alpha'};
M_.param_names_long(2) = {'alpha'};
M_.param_names(3) = {'abar'};
M_.param_names_tex(3) = {'abar'};
M_.param_names_long(3) = {'abar'};
M_.param_names(4) = {'delta'};
M_.param_names_tex(4) = {'delta'};
M_.param_names_long(4) = {'delta'};
M_.param_names(5) = {'gamma'};
M_.param_names_tex(5) = {'gamma'};
M_.param_names_long(5) = {'gamma'};
M_.param_names(6) = {'rho_z'};
M_.param_names_tex(6) = {'rho\_z'};
M_.param_names_long(6) = {'rho_z'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 91;
M_.param_nbr = 6;
M_.orig_endo_nbr = 91;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.orig_eq_nbr = 91;
M_.eq_nbr = 91;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 16 107;
 0 17 0;
 0 18 0;
 1 19 0;
 0 20 0;
 2 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 0 25 0;
 0 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 3 32 0;
 0 33 0;
 0 34 108;
 4 35 0;
 0 36 0;
 0 37 109;
 0 38 0;
 0 39 0;
 0 40 0;
 0 41 0;
 0 42 0;
 0 43 0;
 0 44 0;
 0 45 0;
 0 46 0;
 5 47 0;
 0 48 0;
 0 49 110;
 6 50 0;
 0 51 0;
 0 52 111;
 7 53 0;
 0 54 0;
 0 55 112;
 0 56 0;
 0 57 0;
 0 58 0;
 0 59 0;
 0 60 0;
 0 61 0;
 0 62 0;
 0 63 0;
 0 64 0;
 8 65 0;
 0 66 0;
 0 67 113;
 9 68 0;
 0 69 0;
 0 70 114;
 10 71 0;
 0 72 0;
 0 73 115;
 0 74 0;
 0 75 0;
 0 76 0;
 0 77 0;
 0 78 0;
 0 79 0;
 0 80 0;
 0 81 0;
 0 82 0;
 11 83 0;
 0 84 0;
 0 85 116;
 12 86 0;
 0 87 0;
 0 88 117;
 13 89 0;
 0 90 0;
 0 91 118;
 0 92 0;
 0 93 0;
 0 94 0;
 0 95 0;
 0 96 0;
 0 97 0;
 0 98 0;
 0 99 0;
 0 100 0;
 14 101 0;
 0 102 0;
 0 103 119;
 15 104 0;
 0 105 0;
 0 106 120;]';
M_.nstatic = 62;
M_.nfwrd   = 14;
M_.npred   = 15;
M_.nboth   = 0;
M_.nsfwrd   = 14;
M_.nspred   = 15;
M_.ndynamic   = 29;
M_.dynamic_tmp_nbr = [13; 11; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'c1' ;
  2 , 'name' , 'at1' ;
  3 , 'name' , 'a1' ;
  4 , 'name' , 'c2' ;
  5 , 'name' , 'at2' ;
  6 , 'name' , '6' ;
  7 , 'name' , 'c3' ;
  8 , 'name' , 'at3' ;
  9 , 'name' , 'a3' ;
  10 , 'name' , 'c4' ;
  11 , 'name' , 'at4' ;
  12 , 'name' , 'a4' ;
  13 , 'name' , 'c5' ;
  14 , 'name' , 'at5' ;
  15 , 'name' , 'a5' ;
  16 , 'name' , 'c6' ;
  17 , 'name' , 'at6' ;
  18 , 'name' , 'a6' ;
  19 , 'name' , 'c7' ;
  20 , 'name' , 'at7' ;
  21 , 'name' , '21' ;
  22 , 'name' , 'c8' ;
  23 , 'name' , 'at8' ;
  24 , 'name' , '24' ;
  25 , 'name' , 'c9' ;
  26 , 'name' , 'at9' ;
  27 , 'name' , 'a9' ;
  28 , 'name' , 'c10' ;
  29 , 'name' , 'at10' ;
  30 , 'name' , 'a10' ;
  31 , 'name' , 'c11' ;
  32 , 'name' , 'at11' ;
  33 , 'name' , 'a11' ;
  34 , 'name' , 'c12' ;
  35 , 'name' , 'at12' ;
  36 , 'name' , '36' ;
  37 , 'name' , 'c13' ;
  38 , 'name' , 'at13' ;
  39 , 'name' , '39' ;
  40 , 'name' , 'c14' ;
  41 , 'name' , 'at14' ;
  42 , 'name' , '42' ;
  43 , 'name' , 'c15' ;
  44 , 'name' , 'at15' ;
  45 , 'name' , 'a15' ;
  46 , 'name' , 'c16' ;
  47 , 'name' , 'at16' ;
  48 , 'name' , 'a16' ;
  49 , 'name' , 'c17' ;
  50 , 'name' , 'at17' ;
  51 , 'name' , 'a17' ;
  52 , 'name' , 'c18' ;
  53 , 'name' , 'at18' ;
  54 , 'name' , '54' ;
  55 , 'name' , 'c19' ;
  56 , 'name' , 'at19' ;
  57 , 'name' , '57' ;
  58 , 'name' , 'c20' ;
  59 , 'name' , 'at20' ;
  60 , 'name' , '60' ;
  61 , 'name' , 'c21' ;
  62 , 'name' , 'at21' ;
  63 , 'name' , 'a21' ;
  64 , 'name' , 'c22' ;
  65 , 'name' , 'at22' ;
  66 , 'name' , 'a22' ;
  67 , 'name' , 'c23' ;
  68 , 'name' , 'at23' ;
  69 , 'name' , 'a23' ;
  70 , 'name' , 'c24' ;
  71 , 'name' , 'at24' ;
  72 , 'name' , '72' ;
  73 , 'name' , 'c25' ;
  74 , 'name' , 'at25' ;
  75 , 'name' , '75' ;
  76 , 'name' , 'w' ;
  77 , 'name' , 'r' ;
  78 , 'name' , 'Z' ;
  79 , 'name' , 'u' ;
  80 , 'name' , 'GDP' ;
  81 , 'name' , 'C' ;
  82 , 'name' , 'I' ;
  83 , 'name' , 'K' ;
  84 , 'name' , 'ut' ;
  85 , 'name' , 'Kt' ;
  86 , 'name' , 'Ct' ;
  87 , 'name' , 'wt' ;
  88 , 'name' , 'rt' ;
  89 , 'name' , 'Zt' ;
  90 , 'name' , 'GDPt' ;
  91 , 'name' , 'It' ;
};
M_.mapping.r.eqidx = [1 4 6 7 10 13 16 19 21 22 24 25 28 31 34 36 37 39 40 42 43 46 49 52 54 55 57 58 60 61 64 67 70 72 73 75 76 77 88 ];
M_.mapping.w.eqidx = [1 4 7 10 13 16 19 22 25 28 31 34 37 40 43 46 49 52 55 58 61 64 67 70 73 76 87 ];
M_.mapping.Z.eqidx = [76 77 78 80 ];
M_.mapping.K.eqidx = [77 80 82 83 85 ];
M_.mapping.GDP.eqidx = [80 90 ];
M_.mapping.u.eqidx = [78 79 84 ];
M_.mapping.I.eqidx = [82 91 ];
M_.mapping.C.eqidx = [81 86 ];
M_.mapping.rt.eqidx = [88 ];
M_.mapping.wt.eqidx = [87 ];
M_.mapping.Zt.eqidx = [89 ];
M_.mapping.Kt.eqidx = [85 ];
M_.mapping.GDPt.eqidx = [90 ];
M_.mapping.ut.eqidx = [84 89 ];
M_.mapping.It.eqidx = [91 ];
M_.mapping.Ct.eqidx = [86 ];
M_.mapping.a1.eqidx = [1 2 3 17 83 ];
M_.mapping.at1.eqidx = [1 2 ];
M_.mapping.c1.eqidx = [1 6 81 ];
M_.mapping.a2.eqidx = [2 4 17 83 ];
M_.mapping.at2.eqidx = [4 5 ];
M_.mapping.c2.eqidx = [4 6 21 24 81 ];
M_.mapping.a3.eqidx = [7 9 83 ];
M_.mapping.at3.eqidx = [8 ];
M_.mapping.c3.eqidx = [7 81 ];
M_.mapping.a4.eqidx = [10 12 83 ];
M_.mapping.at4.eqidx = [11 ];
M_.mapping.c4.eqidx = [10 81 ];
M_.mapping.a5.eqidx = [13 15 83 ];
M_.mapping.at5.eqidx = [14 ];
M_.mapping.c5.eqidx = [13 81 ];
M_.mapping.a6.eqidx = [5 16 18 20 35 83 ];
M_.mapping.at6.eqidx = [16 17 ];
M_.mapping.c6.eqidx = [6 16 81 ];
M_.mapping.a7.eqidx = [5 19 20 35 83 ];
M_.mapping.at7.eqidx = [19 20 ];
M_.mapping.c7.eqidx = [19 21 24 81 ];
M_.mapping.a8.eqidx = [5 20 22 35 83 ];
M_.mapping.at8.eqidx = [22 23 ];
M_.mapping.c8.eqidx = [22 24 36 39 42 81 ];
M_.mapping.a9.eqidx = [25 27 83 ];
M_.mapping.at9.eqidx = [26 ];
M_.mapping.c9.eqidx = [25 81 ];
M_.mapping.a10.eqidx = [28 30 83 ];
M_.mapping.at10.eqidx = [29 ];
M_.mapping.c10.eqidx = [28 81 ];
M_.mapping.a11.eqidx = [31 33 83 ];
M_.mapping.at11.eqidx = [32 ];
M_.mapping.c11.eqidx = [31 81 ];
M_.mapping.a12.eqidx = [23 34 38 53 83 ];
M_.mapping.at12.eqidx = [34 35 ];
M_.mapping.c12.eqidx = [21 24 34 36 81 ];
M_.mapping.a13.eqidx = [23 37 38 53 83 ];
M_.mapping.at13.eqidx = [37 38 ];
M_.mapping.c13.eqidx = [36 37 39 42 81 ];
M_.mapping.a14.eqidx = [23 38 40 53 83 ];
M_.mapping.at14.eqidx = [40 41 ];
M_.mapping.c14.eqidx = [40 42 54 57 60 81 ];
M_.mapping.a15.eqidx = [43 45 83 ];
M_.mapping.at15.eqidx = [44 ];
M_.mapping.c15.eqidx = [43 81 ];
M_.mapping.a16.eqidx = [46 48 83 ];
M_.mapping.at16.eqidx = [47 ];
M_.mapping.c16.eqidx = [46 81 ];
M_.mapping.a17.eqidx = [49 51 83 ];
M_.mapping.at17.eqidx = [50 ];
M_.mapping.c17.eqidx = [49 81 ];
M_.mapping.a18.eqidx = [41 52 56 71 83 ];
M_.mapping.at18.eqidx = [52 53 ];
M_.mapping.c18.eqidx = [36 39 42 52 54 81 ];
M_.mapping.a19.eqidx = [41 55 56 71 83 ];
M_.mapping.at19.eqidx = [55 56 ];
M_.mapping.c19.eqidx = [54 55 57 60 81 ];
M_.mapping.a20.eqidx = [41 56 58 71 83 ];
M_.mapping.at20.eqidx = [58 59 ];
M_.mapping.c20.eqidx = [58 60 72 75 81 ];
M_.mapping.a21.eqidx = [61 63 83 ];
M_.mapping.at21.eqidx = [62 ];
M_.mapping.c21.eqidx = [61 81 ];
M_.mapping.a22.eqidx = [64 66 83 ];
M_.mapping.at22.eqidx = [65 ];
M_.mapping.c22.eqidx = [64 81 ];
M_.mapping.a23.eqidx = [67 69 83 ];
M_.mapping.at23.eqidx = [68 ];
M_.mapping.c23.eqidx = [67 81 ];
M_.mapping.a24.eqidx = [59 70 74 83 ];
M_.mapping.at24.eqidx = [70 71 ];
M_.mapping.c24.eqidx = [54 57 60 70 72 81 ];
M_.mapping.a25.eqidx = [59 73 74 83 ];
M_.mapping.at25.eqidx = [73 74 ];
M_.mapping.c25.eqidx = [72 73 75 81 ];
M_.mapping.eps.eqidx = [79 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 6 17 20 32 35 38 50 53 56 68 71 74 86 89 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(91, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.endo_trends = struct('deflator', cell(91, 1), 'log_deflator', cell(91, 1), 'growth_factor', cell(91, 1), 'log_growth_factor', cell(91, 1));
M_.NNZDerivatives = [324; -1; -1; ];
M_.static_tmp_nbr = [19; 14; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(2) = 0.36;
alpha = M_.params(2);
M_.params(1) = 0.98;
beta = M_.params(1);
M_.params(3) = 0;
abar = M_.params(3);
M_.params(4) = 0.025;
delta = M_.params(4);
M_.params(5) = 1;
gamma = M_.params(5);
M_.params(6) = 0.95;
rho_z = M_.params(6);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(17) = 0.000000000000;
oo_.steady_state(18) = 0.167711210000;
oo_.steady_state(19) = 0.560187572721;
oo_.steady_state(20) = 8.689700000000;
oo_.steady_state(21) = 8.802382240000;
oo_.steady_state(22) = 0.648975637279;
oo_.steady_state(23) = 0.000000000000;
oo_.steady_state(24) = 0.000000000000;
oo_.steady_state(25) = 0.389683004549;
oo_.steady_state(26) = 0.000000000000;
oo_.steady_state(27) = 0.000000000000;
oo_.steady_state(28) = 0.389683004549;
oo_.steady_state(29) = 0.000000000000;
oo_.steady_state(30) = 0.000000000000;
oo_.steady_state(31) = 0.389683004549;
oo_.steady_state(32) = 0.000000000000;
oo_.steady_state(33) = 0.167711210000;
oo_.steady_state(34) = 0.975705914587;
oo_.steady_state(35) = 8.689900000000;
oo_.steady_state(36) = 8.802382240000;
oo_.steady_state(37) = 1.064293979146;
oo_.steady_state(38) = 19.397800000000;
oo_.steady_state(39) = 19.692585280000;
oo_.steady_state(40) = 1.427981660381;
oo_.steady_state(41) = 0.000000000000;
oo_.steady_state(42) = 0.000000000000;
oo_.steady_state(43) = 0.805201346415;
oo_.steady_state(44) = 0.000000000000;
oo_.steady_state(45) = 0.000000000000;
oo_.steady_state(46) = 0.805201346415;
oo_.steady_state(47) = 0.000000000000;
oo_.steady_state(48) = 0.000000000000;
oo_.steady_state(49) = 1.664226312243;
oo_.steady_state(50) = 8.761300000000;
oo_.steady_state(51) = 8.802382240000;
oo_.steady_state(52) = 1.851918944974;
oo_.steady_state(53) = 19.513600000000;
oo_.steady_state(54) = 19.692585280000;
oo_.steady_state(55) = 2.171206626209;
oo_.steady_state(56) = 48.910200000000;
oo_.steady_state(57) = 49.343458720000;
oo_.steady_state(58) = 2.919337987999;
oo_.steady_state(59) = 0.000000000000;
oo_.steady_state(60) = 0.000000000000;
oo_.steady_state(61) = 1.664226312243;
oo_.steady_state(62) = 0.000000000000;
oo_.steady_state(63) = 0.000000000000;
oo_.steady_state(64) = 3.439544574955;
oo_.steady_state(65) = 0.000000000000;
oo_.steady_state(66) = 0.000000000000;
oo_.steady_state(67) = 3.439544574955;
oo_.steady_state(68) = 20.035400000000;
oo_.steady_state(69) = 19.692585280000;
oo_.steady_state(70) = 3.424724888920;
oo_.steady_state(71) = 49.456000000000;
oo_.steady_state(72) = 49.343458720000;
oo_.steady_state(73) = 4.148856250710;
oo_.steady_state(74) = 114.271700000000;
oo_.steady_state(75) = 114.544305600000;
oo_.steady_state(76) = 5.619972986008;
oo_.steady_state(77) = 0.000000000000;
oo_.steady_state(78) = 0.000000000000;
oo_.steady_state(79) = 7.109669535474;
oo_.steady_state(80) = 0.000000000000;
oo_.steady_state(81) = 0.000000000000;
oo_.steady_state(82) = 7.109669535474;
oo_.steady_state(83) = 0.000000000000;
oo_.steady_state(84) = 0.000000000000;
oo_.steady_state(85) = 7.109669535474;
oo_.steady_state(86) = 50.987100000000;
oo_.steady_state(87) = 49.343458720000;
oo_.steady_state(88) = 6.287881211229;
oo_.steady_state(89) = 115.795100000000;
oo_.steady_state(90) = 114.544305600000;
oo_.steady_state(91) = 7.766697946527;
oo_.steady_state(1) = 0.016655763035;
oo_.steady_state(2) = 2.152944776511;
oo_.steady_state(4) = 29.072362347015;
oo_.steady_state(3) = 1;
oo_.steady_state(6) = 0;
oo_.steady_state(5) = 3.363976213299;
oo_.steady_state(7) = 0.726809058675;
oo_.steady_state(8) = 2.637169108041;
oo_.steady_state(14) = 0;
oo_.steady_state(12) = 0;
oo_.steady_state(16) = 0;
oo_.steady_state(10) = 0;
oo_.steady_state(9) = 0;
oo_.steady_state(11) = 0;
oo_.steady_state(13) = 0;
oo_.steady_state(15) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.003120000000)^2;
options_.irf = 200;
options_.order = 1;
options_.periods = 10000;
var_list_ = {'ut';'Zt';'GDPt';'Ct';'Kt';'It';'rt';'wt';'u';'Z';'GDP';'C';'K';'I';'r';'w'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('DyTruncation_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('DyTruncation_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('DyTruncation_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('DyTruncation_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('DyTruncation_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('DyTruncation_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('DyTruncation_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
