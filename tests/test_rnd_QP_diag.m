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

%TEST: randomly generated diagonal QP problem
%  Using single-block implementation
%
% Output: rndQP_out_multi, rndQP_out_single
%

clear all
%construct a problem 
n_var= 500;
problem_construct.Q_sparsity= 0.9000;
problem_construct.n_var= n_var;
problem_construct.c_sparsity= 0.1000;
problem_construct.Aeq_sparsity= 0.9000;
problem_construct.Aineq_sparsity= 0.9000;
problem_construct.Aineq_n_row= 0;
problem_construct.Aeq_n_row= 30;
problem_construct.rnd_seed= 1234;

Q_construct.cond_num= 5;
Q_construct.zeta= 2;
Q_construct.eta= 0.5000;

%To use simple Q = sprand(.) call, no bounds
% m_qp = get_rnd_QP(problem_construct);
%
% to use a controlled way to construct Q call, with bounds
m_qp = get_model_rnd_QP(problem_construct,true);
%replace Q by a diagonal matrix
m_qp.Q=sparse(diag(rand(n_var,1)));
m_qp.c = sparse(n_var,1);

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = 5;
run_params.beta = 1;
run_params.epsilon = 1e-6;
run_params.n_blocks = 1;
%turn off verbose
run_params.debug = 0;

%change to single block run
disp('################')
disp('Running RAC single-block (LDL)')
run_params.n_blocks = 1;

%default Cholesky
rndQP_out_single_ldl = RACQP(m_qp,run_params);

% use diagonal sparse KKT implementation 
disp('################')
disp('Running RAC single-block (KKT)')
run_params.single_block_solver = 'diag_kkt';
rndQP_out_single_kkt = RACQP(m_qp,run_params);

% use diagonal sparse cholesky implementation 
run_params.single_block_solver = 'diag_cholesky';
rndQP_out_single_diagC = RACQP(m_qp,run_params);

disp('Solution, single-block run (LDL):')
rndQP_out_single_ldl
disp('Solution, single-block run (KKT):')
rndQP_out_single_kkt


