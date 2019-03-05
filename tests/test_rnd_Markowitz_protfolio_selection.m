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

%TEST: randomly generated Markowitz portfolio selection problem (binary)
%
% Output: rnd_mark_bin

clear all
%construct a problem 
sparsity = .1;
kappa = 1e-5;
rnd_seed = 123;
Q_construct.cond_num= 5;
Q_construct.zeta= 2;
Q_construct.eta= 0.5000;
nsize = 500;
p_sel = 0.5;
nselected = min(nsize,ceil(p_sel*nsize));

model = get_model_rnd_Markowitz(nsize, sparsity, kappa, Q_construct, rnd_seed, nselected);

model.local_constraints.lb=model.lb;
model.local_constraints.ub=model.ub;
model.local_constraints.Aineq = sparse(0,nsize);
model.local_constraints.bineq = [];
model.local_constraints.Aeq = sparse(0,nsize);
model.local_constraints.beq = [];

%load default runtime parameters to be used when solving each iteration of MIP
default_run_params
run_params.n_blocks = 5;
run_params.beta = 1;
run_params.epsilon = 1e-7;
run_params.sub_solver_type='gurobi';
%turn off verbose
run_params.debug = 0;


%load default gurobi parameters
default_gurobi_parameters
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 5;
 
%gurobi_params.presolve=1;
%gurobi_params.outputflag=1;
run_params.gurobi_params = gurobi_params;

%set MIP run parameters

run_mip.permute_type = 'permute';
run_mip.permute_dist = 'exponential'; 
run_mip.permute_min = 2;
run_mip.permute_max = nsize;
run_mip.permute_mu = 0.4*nsize; 
run_mip.max_iter = 10000000; %time limited 
run_mip.max_nperturb = 10000000; %time limited
run_mip.max_rtime = 10; 
run_mip.rnd_seed = 123; 
run_mip.n_perturb_trial = max(2,round(0.1*nsize));
run_mip.run_sub = run_params;
run_mip.debug = 1;



%run the solver
disp('Running RAC MIP')
rnd_mark_bin = RACQP(model,run_mip);

disp("SOLUTION")
disp(rnd_mark_bin)
disp("SUM OF SOL")
disp(nnz(round(rnd_mark_bin.sol_x(1:nsize-1))))






