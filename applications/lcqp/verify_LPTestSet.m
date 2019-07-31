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

% Verify results: LPTestSet benchmark problems + diag Q
%      http://plato.asu.edu/ftp/lptestset/

function [sol, model]=verify_LPTestSet(filename, beta, n_blocks, pL,solver) 
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_MPS(filename,true);
if(nargin > 4 && strcmp(solver,'gurobi'))
  disp('Running Gurobi')
  sol = GUROBI(model);
  disp("Done")
  return;
else if(nargin > 4 && strcmp(solver,'mosek'))
  disp('Running Mosek')
  sol = MOSEK(model);
  disp("Done")
  return;
end
if(strcmpi("l_bounds",pL))
  solver = 'gurobi';
  model.local_constraints.lb=model.lb;
  model.local_constraints.ub=model.ub;
  model.lb = -inf*model.lb;
  model.ub = inf*model.ub;
  model.local_constraints.Aineq = sparse(0,model.size);
  model.local_constraints.bineq = [];
  model.local_constraints.Aeq = sparse(0,model.size);
  model.local_constraints.beq = [];
elseif(strcmpi("none",pL))
  solver = 'cholesky';
else
  error("Wrong part. Lagrange mode")
end

%load default runtime parameters
run_params = default_run_params;
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
run_params.epsilon = 1e-4;
run_params.max_iter = 4000;
run_params.max_rtime = 1800;
run_params.sub_solver_type = solver;

gurobi_params = default_gurobi_parameters;
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 10;
run_params.gurobi_params = gurobi_params;

run_params.use_sparse = true; %this is the default setings!
%turn off verbose
run_params.debug = 0;

%print run parameters
%disp("RUN PARAMETERS")
%disp("run_params")
%disp(run_params)

disp('Running RACQP')
sol = RACQP(model,run_params,true);
disp("Done");
%disp('Solution:')
%disp(sol)

%quit
end
