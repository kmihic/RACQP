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
%TEST: Markowitz min-variance problem
%      Data source Center for Research in Security Price (CRSP)
%      Time period Jan-Dec 2018
%
% 
%

function [sol, model]=verify_Markowitz_min_variance_data(filename, low_rank_s,solver)
addpath('../racqp');
addpath('../utils');

if(strcmpi("true",low_rank_s))
  low_rank = true;
else
  low_rank = false;
end
%setup problem data
if(low_rank)
  n_blocks = 1;
else
  n_blocks = 50;
end
beta = 1;
min_iter = 50;
epsilon = 1e-4;
kappa = 1e-5;

%get the problem 
disp("LOADING THE MODEL...")
model = load_Markowitz_model(filename, -1, kappa, low_rank, false);

if(nargin > 2 && strcmp(solver,'gurobi'))
  disp('Running Gurobi')
  sol = GUROBI(model);
  disp("Done")
  return;
else if(nargin > 2 && strcmp(solver,'mosek'))
  disp('Running Mosek')
  sol = MOSEK(model);
  disp("Done")
  return;
end

nsize=model.size;
model.local_constraints.lb=model.lb;
model.local_constraints.ub=model.ub;
model.local_constraints.Aineq = sparse(0,nsize);
model.local_constraints.bineq = [];
model.local_constraints.Aeq = sparse(0,nsize);
model.local_constraints.beq = [];

%load default runtime parameters to be used when solving each iteration of MIP
run_params = default_run_params;
%change some parameters
run_params.beta = beta;
run_params.epsilon = epsilon;
run_params.rnd_seed = 123;
%turn off verbose
run_params.debug = 0;
if(low_rank)
  %if using lowrank model, we use single block, kkt mode
  %with localized constraints. Q is diagonal, so model is sparse
  run_params.n_blocks = n_blocks;
  run_params.single_block_solver = 'diag_kkt';
  run_params.use_sparse = true;
else
  %regular model has dense Q+A'A
  run_params.use_sparse = false;
end
run_params.n_blocks = n_blocks;
%when using single block, primal feasibility can be found, by a luck
%in the first couple of iterations. However, the solution is far away
%from the optimal. Making sure we do at least a couple of rounds
run_params.min_iter = min_iter;

%print run parameters
%disp("RUN PARAMETERS")
%disp("run_params")
%disp(run_params)

%run the solver
disp('Running RACQP')
sol = RACQP(model,run_params, true);
disp("Done")
%format long
%disp("SOLUTION")
%disp(sol)

%quit
end







