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

%TEST: randomly generated QP problem
%
% Output: rndQP_out_multi, rndQP_out_single
%

clear all
%construct a problem 
problem_construct.Q_sparsity= 0.8000;
problem_construct.n_var= 500;
problem_construct.c_sparsity= 0.8000;
problem_construct.Aeq_sparsity= 0.9000;
problem_construct.Aineq_sparsity= 0.9000;
problem_construct.Aineq_n_row= 4;
problem_construct.Aeq_n_row= 30;
problem_construct.rnd_seed= 1234;

Q_construct.cond_num= 5;
Q_construct.zeta= 2;
Q_construct.eta= 0.5000;

%To use simple Q = sprand(.) call, no bounds
% m_qp = get_rnd_QP(problem_construct);
%
% to use a controlled way to construct Q call, with bounds
m_qp = get_model_rnd_QP(problem_construct,true,Q_construct);

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = 5;
run_params.beta = 1;
run_params.epsilon = 1e-6;
%run_params.max_iter = 100;
%turn off verbose
run_params.debug = 0;

%run the solver
disp('################')
disp('Running RAC multi-block')
rndQP_out_multi = RACQP(m_qp,run_params);

%change to single block run
disp('################')
disp('Running RAC single-block')
run_params.n_blocks = 1;

% use default cholesky 
rndQP_out_single_ldl = RACQP(m_qp,run_params);

% use gurobi to handle bounds
%change to single block run
disp('################')
disp('Running RAC single-block by gurobi')
run_params.single_block_solver = 'bounds_by_gurobi';
default_gurobi_parameters
run_params.gurobi_params = gurobi_params;
rndQP_out_single_gurobi = RACQP(m_qp,run_params);


disp('################')
disp('Solution, multi-block run:')
rndQP_out_multi
disp('Solution, single-block run (LDL):')
rndQP_out_single_ldl
disp('Solution, single-block run (Gurobi):')
rndQP_out_single_gurobi

