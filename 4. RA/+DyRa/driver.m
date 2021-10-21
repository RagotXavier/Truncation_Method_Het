%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'DyRa';
M_.dynare_version = '4.6.2';
oo_.dynare_version = '4.6.2';
options_.dynare_version = '4.6.2';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('DyRa.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps'};
M_.exo_names_tex(1) = {'eps'};
M_.exo_names_long(1) = {'eps'};
M_.endo_names = cell(18,1);
M_.endo_names_tex = cell(18,1);
M_.endo_names_long = cell(18,1);
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
M_.endo_names(9) = {'W'};
M_.endo_names_tex(9) = {'W'};
M_.endo_names_long(9) = {'W'};
M_.endo_names(10) = {'rt'};
M_.endo_names_tex(10) = {'rt'};
M_.endo_names_long(10) = {'rt'};
M_.endo_names(11) = {'wt'};
M_.endo_names_tex(11) = {'wt'};
M_.endo_names_long(11) = {'wt'};
M_.endo_names(12) = {'Zt'};
M_.endo_names_tex(12) = {'Zt'};
M_.endo_names_long(12) = {'Zt'};
M_.endo_names(13) = {'Kt'};
M_.endo_names_tex(13) = {'Kt'};
M_.endo_names_long(13) = {'Kt'};
M_.endo_names(14) = {'GDPt'};
M_.endo_names_tex(14) = {'GDPt'};
M_.endo_names_long(14) = {'GDPt'};
M_.endo_names(15) = {'ut'};
M_.endo_names_tex(15) = {'ut'};
M_.endo_names_long(15) = {'ut'};
M_.endo_names(16) = {'It'};
M_.endo_names_tex(16) = {'It'};
M_.endo_names_long(16) = {'It'};
M_.endo_names(17) = {'Ct'};
M_.endo_names_tex(17) = {'Ct'};
M_.endo_names_long(17) = {'Ct'};
M_.endo_names(18) = {'Wt'};
M_.endo_names_tex(18) = {'Wt'};
M_.endo_names_long(18) = {'Wt'};
M_.endo_partitions = struct();
M_.param_names = cell(8,1);
M_.param_names_tex = cell(8,1);
M_.param_names_long = cell(8,1);
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
M_.param_names(6) = {'rho_u'};
M_.param_names_tex(6) = {'rho\_u'};
M_.param_names_long(6) = {'rho_u'};
M_.param_names(7) = {'coef'};
M_.param_names_tex(7) = {'coef'};
M_.param_names_long(7) = {'coef'};
M_.param_names(8) = {'coef2'};
M_.param_names_tex(8) = {'coef2'};
M_.param_names_long(8) = {'coef2'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 18;
M_.param_nbr = 8;
M_.orig_endo_nbr = 18;
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
M_.orig_eq_nbr = 18;
M_.eq_nbr = 18;
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
 0 3 21;
 0 4 0;
 0 5 0;
 1 6 22;
 0 7 0;
 2 8 0;
 0 9 0;
 0 10 23;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;
 0 19 0;
 0 20 0;]';
M_.nstatic = 14;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [4; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'w' ;
  4 , 'name' , '4' ;
  5 , 'name' , 'Z' ;
  6 , 'name' , 'u' ;
  7 , 'name' , 'GDP' ;
  8 , 'name' , 'I' ;
  9 , 'name' , 'W' ;
  10 , 'name' , 'ut' ;
  11 , 'name' , 'Kt' ;
  12 , 'name' , 'Ct' ;
  13 , 'name' , 'wt' ;
  14 , 'name' , 'rt' ;
  15 , 'name' , 'Zt' ;
  16 , 'name' , 'GDPt' ;
  17 , 'name' , 'It' ;
  18 , 'name' , 'Wt' ;
};
M_.mapping.r.eqidx = [1 2 4 14 ];
M_.mapping.w.eqidx = [1 3 13 ];
M_.mapping.Z.eqidx = [3 4 5 7 ];
M_.mapping.K.eqidx = [1 3 4 7 8 11 ];
M_.mapping.GDP.eqidx = [7 16 ];
M_.mapping.u.eqidx = [5 6 10 ];
M_.mapping.I.eqidx = [8 17 ];
M_.mapping.C.eqidx = [1 2 9 12 ];
M_.mapping.W.eqidx = [9 18 ];
M_.mapping.rt.eqidx = [14 ];
M_.mapping.wt.eqidx = [13 ];
M_.mapping.Zt.eqidx = [15 ];
M_.mapping.Kt.eqidx = [11 ];
M_.mapping.GDPt.eqidx = [16 ];
M_.mapping.ut.eqidx = [10 15 ];
M_.mapping.It.eqidx = [17 ];
M_.mapping.Ct.eqidx = [12 ];
M_.mapping.Wt.eqidx = [18 ];
M_.mapping.eps.eqidx = [6 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 6 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(18, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(8, 1);
M_.endo_trends = struct('deflator', cell(18, 1), 'log_deflator', cell(18, 1), 'growth_factor', cell(18, 1), 'log_growth_factor', cell(18, 1));
M_.NNZDerivatives = [45; -1; -1; ];
M_.static_tmp_nbr = [3; 2; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(2) = 0.36;
alpha = M_.params(2);
M_.params(1) = 0.98;
beta = M_.params(1);
M_.params(3) = 4.0945e-05;
abar = M_.params(3);
M_.params(4) = 0.025;
delta = M_.params(4);
M_.params(5) = 1.0001;
gamma = M_.params(5);
M_.params(6) = 0.95;
rho_u = M_.params(6);
M_.params(7) = 1.0;
coef = M_.params(7);
M_.params(8) = 1.;
coef2 = M_.params(8);
resid;
options_.steadystate.nocheck = true;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.00312000)^2;
options_.TeX=1;
options_.irf = 200;
options_.order = 1;
options_.periods = 10000;
var_list_ = {'ut';'Zt';'GDPt';'Ct';'Kt';'It';'Wt';'rt';'wt';'u';'Z';'GDP';'C';'K';'I';'W';'r';'w'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('DyRa_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('DyRa_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('DyRa_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('DyRa_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('DyRa_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('DyRa_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('DyRa_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
