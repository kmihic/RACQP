%
% RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
% Copyright (C) 2019 
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>
%
% This file is part of RACQP 
%
%

%TEST: randomly generated Markowitz min variance problem
%
% Output: rndMark_var_out_multi, rndMark_var_out_single
%

%construct a problem 
clear all
size = 500;
sparsity = .1;
mu = 1e-5;
rnd_seed = 123;
Q_construct.cond_num= 5;
Q_construct.zeta= 2;
Q_construct.eta= 0.5000;

m_mark_var = get_model_rnd_Markowitz(size, sparsity, mu, Q_construct, rnd_seed);

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = 5;
run_params.beta = 1;
run_params.epsilon = 1e-6;
%turn off verbose
%run_params.debug = 0;

%run the solver
disp('Running RAC multi-block')
rndMark_var_out_multi = RACQP(m_mark_var,run_params);

%change to single block run
disp('Running RAC single-block')
run_params.n_blocks = 1;
rndMark_var_out_single = RACQP(m_mark_var,run_params);


disp('Solution, multi-block run:')
rndMark_var_out_multi
disp('Solution, single-block run:')
rndMark_var_out_single


