%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'DyRA';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('DyRA.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps'};
M_.exo_names_tex(1) = {'eps'};
M_.exo_names_long(1) = {'eps'};
M_.endo_names = cell(22,1);
M_.endo_names_tex = cell(22,1);
M_.endo_names_long = cell(22,1);
M_.endo_names(1) = {'R'};
M_.endo_names_tex(1) = {'R'};
M_.endo_names_long(1) = {'R'};
M_.endo_names(2) = {'w'};
M_.endo_names_tex(2) = {'w'};
M_.endo_names_long(2) = {'w'};
M_.endo_names(3) = {'Z'};
M_.endo_names_tex(3) = {'Z'};
M_.endo_names_long(3) = {'Z'};
M_.endo_names(4) = {'a'};
M_.endo_names_tex(4) = {'a'};
M_.endo_names_long(4) = {'a'};
M_.endo_names(5) = {'K'};
M_.endo_names_tex(5) = {'K'};
M_.endo_names_long(5) = {'K'};
M_.endo_names(6) = {'GDP'};
M_.endo_names_tex(6) = {'GDP'};
M_.endo_names_long(6) = {'GDP'};
M_.endo_names(7) = {'u'};
M_.endo_names_tex(7) = {'u'};
M_.endo_names_long(7) = {'u'};
M_.endo_names(8) = {'I'};
M_.endo_names_tex(8) = {'I'};
M_.endo_names_long(8) = {'I'};
M_.endo_names(9) = {'c'};
M_.endo_names_tex(9) = {'c'};
M_.endo_names_long(9) = {'c'};
M_.endo_names(10) = {'C'};
M_.endo_names_tex(10) = {'C'};
M_.endo_names_long(10) = {'C'};
M_.endo_names(11) = {'W'};
M_.endo_names_tex(11) = {'W'};
M_.endo_names_long(11) = {'W'};
M_.endo_names(12) = {'TT'};
M_.endo_names_tex(12) = {'TT'};
M_.endo_names_long(12) = {'TT'};
M_.endo_names(13) = {'Rt'};
M_.endo_names_tex(13) = {'Rt'};
M_.endo_names_long(13) = {'Rt'};
M_.endo_names(14) = {'wt'};
M_.endo_names_tex(14) = {'wt'};
M_.endo_names_long(14) = {'wt'};
M_.endo_names(15) = {'Zt'};
M_.endo_names_tex(15) = {'Zt'};
M_.endo_names_long(15) = {'Zt'};
M_.endo_names(16) = {'Kt'};
M_.endo_names_tex(16) = {'Kt'};
M_.endo_names_long(16) = {'Kt'};
M_.endo_names(17) = {'GDPt'};
M_.endo_names_tex(17) = {'GDPt'};
M_.endo_names_long(17) = {'GDPt'};
M_.endo_names(18) = {'ut'};
M_.endo_names_tex(18) = {'ut'};
M_.endo_names_long(18) = {'ut'};
M_.endo_names(19) = {'It'};
M_.endo_names_tex(19) = {'It'};
M_.endo_names_long(19) = {'It'};
M_.endo_names(20) = {'Ct'};
M_.endo_names_tex(20) = {'Ct'};
M_.endo_names_long(20) = {'Ct'};
M_.endo_names(21) = {'Wt'};
M_.endo_names_tex(21) = {'Wt'};
M_.endo_names_long(21) = {'Wt'};
M_.endo_names(22) = {'Tt'};
M_.endo_names_tex(22) = {'Tt'};
M_.endo_names_long(22) = {'Tt'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'alpha'};
M_.param_names_tex(2) = {'alpha'};
M_.param_names_long(2) = {'alpha'};
M_.param_names(3) = {'theta'};
M_.param_names_tex(3) = {'theta'};
M_.param_names_long(3) = {'theta'};
M_.param_names(4) = {'abar'};
M_.param_names_tex(4) = {'abar'};
M_.param_names_long(4) = {'abar'};
M_.param_names(5) = {'delta'};
M_.param_names_tex(5) = {'delta'};
M_.param_names_long(5) = {'delta'};
M_.param_names(6) = {'gamma'};
M_.param_names_tex(6) = {'gamma'};
M_.param_names_long(6) = {'gamma'};
M_.param_names(7) = {'rho_u'};
M_.param_names_tex(7) = {'rho\_u'};
M_.param_names_long(7) = {'rho_u'};
M_.param_names(8) = {'tau'};
M_.param_names_tex(8) = {'tau'};
M_.param_names_long(8) = {'tau'};
M_.param_names(9) = {'coef'};
M_.param_names_tex(9) = {'coef'};
M_.param_names_long(9) = {'coef'};
M_.param_names(10) = {'coef2'};
M_.param_names_tex(10) = {'coef2'};
M_.param_names_long(10) = {'coef2'};
M_.param_names(11) = {'coefz'};
M_.param_names_tex(11) = {'coefz'};
M_.param_names_long(11) = {'coefz'};
M_.param_names(12) = {'coefo'};
M_.param_names_tex(12) = {'coefo'};
M_.param_names_long(12) = {'coefo'};
M_.param_names(13) = {'coeft'};
M_.param_names_tex(13) = {'coeft'};
M_.param_names_long(13) = {'coeft'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 22;
M_.param_nbr = 13;
M_.orig_endo_nbr = 22;
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
M_.orig_eq_nbr = 22;
M_.eq_nbr = 22;
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
 0 4 26;
 0 5 0;
 0 6 0;
 1 7 0;
 2 8 27;
 0 9 0;
 3 10 0;
 0 11 0;
 0 12 28;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 0;
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 0 25 0;]';
M_.nstatic = 17;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [4; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'c' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'a' ;
  4 , 'name' , 'w' ;
  5 , 'name' , '5' ;
  6 , 'name' , 'u' ;
  7 , 'name' , 'Z' ;
  8 , 'name' , 'GDP' ;
  9 , 'name' , 'I' ;
  10 , 'name' , 'C' ;
  11 , 'name' , 'TT' ;
  12 , 'name' , 'W' ;
  13 , 'name' , 'ut' ;
  14 , 'name' , 'Kt' ;
  15 , 'name' , 'Ct' ;
  16 , 'name' , 'wt' ;
  17 , 'name' , 'Rt' ;
  18 , 'name' , 'Zt' ;
  19 , 'name' , 'GDPt' ;
  20 , 'name' , 'It' ;
  21 , 'name' , 'Wt' ;
  22 , 'name' , 'Tt' ;
};
M_.mapping.R.eqidx = [1 2 4 5 17 ];
M_.mapping.w.eqidx = [1 4 16 ];
M_.mapping.Z.eqidx = [4 5 7 8 ];
M_.mapping.a.eqidx = [1 3 ];
M_.mapping.K.eqidx = [3 5 8 9 11 14 ];
M_.mapping.GDP.eqidx = [8 10 11 19 ];
M_.mapping.u.eqidx = [6 7 13 ];
M_.mapping.I.eqidx = [9 10 20 ];
M_.mapping.c.eqidx = [1 2 12 15 ];
M_.mapping.C.eqidx = [10 ];
M_.mapping.W.eqidx = [12 21 ];
M_.mapping.TT.eqidx = [1 10 11 12 22 ];
M_.mapping.Rt.eqidx = [17 ];
M_.mapping.wt.eqidx = [16 ];
M_.mapping.Zt.eqidx = [18 ];
M_.mapping.Kt.eqidx = [14 ];
M_.mapping.GDPt.eqidx = [19 ];
M_.mapping.ut.eqidx = [13 18 ];
M_.mapping.It.eqidx = [20 ];
M_.mapping.Ct.eqidx = [15 ];
M_.mapping.Wt.eqidx = [21 ];
M_.mapping.Tt.eqidx = [22 ];
M_.mapping.eps.eqidx = [6 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 5 7 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(22, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(22, 1), 'log_deflator', cell(22, 1), 'growth_factor', cell(22, 1), 'log_growth_factor', cell(22, 1));
M_.NNZDerivatives = [58; -1; -1; ];
M_.static_tmp_nbr = [4; 1; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(2) = 0.36;
alpha = M_.params(2);
M_.params(1) = 0.99001;
beta = M_.params(1);
M_.params(3) = 0.23621;
theta = M_.params(3);
M_.params(8) = 0.080085;
tau = M_.params(8);
M_.params(4) = 6.2852e-11;
abar = M_.params(4);
M_.params(5) = 0.025;
delta = M_.params(5);
M_.params(6) = 1.0001;
gamma = M_.params(6);
M_.params(7) = 0.95;
rho_u = M_.params(7);
M_.params(11) = (-0.039898);
coefz = M_.params(11);
M_.params(12) = 0.068017;
coefo = M_.params(12);
M_.params(13) = 0.0020753;
coeft = M_.params(13);
M_.params(9) = 1.0;
coef = M_.params(9);
M_.params(10) = 1.;
coef2 = M_.params(10);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(4) = 40.59010827;
oo_.steady_state(9) = 2.47487278;
oo_.steady_state(1) = 1.00864446;
oo_.steady_state(2) = 2.42779073;
oo_.steady_state(5) = 40.59010827;
oo_.steady_state(3) = 1;
oo_.steady_state(7) = 0;
oo_.steady_state(6) = 3.79342302;
oo_.steady_state(8) = 1.01475271;
oo_.steady_state(10) = 2.47540948;
oo_.steady_state(12) = 0.30233731;
oo_.steady_state(11) = 11.46959010;
oo_.steady_state(18) = 0;
oo_.steady_state(16) = 0;
oo_.steady_state(20) = 0;
oo_.steady_state(14) = 0;
oo_.steady_state(13) = 0;
oo_.steady_state(15) = 0;
oo_.steady_state(17) = 0;
oo_.steady_state(19) = 0;
oo_.steady_state(21) = 0;
oo_.steady_state(22) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid;
options_.solve_tolf=10^-6;
options_.solve_algo = 3;
options_.steady.maxit = 100;
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.00312000)^2;
options_.TeX=1;
options_.irf = 200;
options_.order = 1;
options_.periods = 1000;
var_list_ = {'Zt';'GDPt';'Ct';'Kt';'It';'Wt';'Rt';'wt';'Tt';'Z';'GDP';'C';'K';'I';'W';'R';'w';'TT'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('DyRA_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('DyRA_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('DyRA_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('DyRA_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('DyRA_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('DyRA_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('DyRA_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
