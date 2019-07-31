%
% RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
% Copyright (C) 2019 
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>
%
% This file is part of RACQP 
%

% Verify results: CUTEr benchmark problems
%      Loading data from OSQP solver dataset (mat files)
%      Original mps files:
%      http://www.cuter.rl.ac.uk/Problems/marmes.html
%

function [sol, model]=verify_CUTEr(filename, beta, pL_mode, scale, solver) 
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_OSQP(filename);

if(strcmpi('true',scale))
  [model, scale_Q] = scale_model(model);
else
  scale_Q = inf;
end

%verify the model (some simple checks)
if(~verify_model(model))
  error('Error. Model has problems')
end
%clean the model from embedded, unbounded slacks
model = remove_embedded_slacks(model);
if(nargin > 4 && strcmp(solver,'gurobi'))
  disp('Running Gurobi')
  sol = GUROBI(model);
  disp("Done")
  if(~isinf(scale_Q))
    sol.sol_obj_val = sol.sol_obj_val * scale_Q;
  end
  return;
elseif(nargin > 4 && strcmp(solver,'mosek'))
  disp('Running Mosek')
  sol = MOSEK(model);
  disp("Done")
  if(~isinf(scale_Q))
    sol.sol_obj_val = sol.sol_obj_val * scale_Q;
  end
  return;
end

%load default runtime parameters
run_params = default_run_params;
%change some parameters
run_params.n_blocks = 1;
m = size(model.Aeq,1) + size(model.Aineq,1);
if(strcmpi('m',beta))
  run_params.beta = m;
elseif(strcmpi('m3',beta))
  run_params.beta = m^3;
else
  error('beta can be m or m3')
end
run_params.epsilon = 1e-4;
run_params.max_iter = 4000;
run_params.n_blocks = 1;
run_params.max_rtime = 100;
if(strcmpi('l_eq', pL_mode))
  run_params.single_block_solver = 'diag_kkt';
elseif(strcmpi('l_bounds',pL_mode))
  run_params.single_block_solver = 'bounds_by_gurobi';  
  % no need for extra iterations
  gurobi_params = default_gurobi_parameters;
  %change max time for gurobi for each subproblem
  gurobi_params.TimeLimit = 200;
  run_params.gurobi_params = gurobi_params;
  run_params.min_iter = 0;
elseif(strcmpi('none', pL_mode))
  run_params.single_block_solver = 'cholesky';
end

run_params.use_sparse = true; %this is the default setings!
%turn off verbose
run_params.debug = 0;

%print run parameters
%disp("RUN PARAMETERS")
%disp("run_params")
%disp(run_params)

disp('Running RACQP single-block')
sol = RACQP(model,run_params,true);
disp("Done");

%if scalled, get the proper objVal
 if(~isinf(scale_Q))
   sol.sol_obj_val = cuter_out.sol_obj_val * scale_Q;
 end

%disp('Solution:')
%disp(cuter_out)

%quit
end

function [m, max_Q] = scale_model(m)

max_Q = max(max(abs(m.Q)));
m.Q = m.Q/max_Q;
m.c = m.c/max_Q;
if(size(m.Aeq,1)>0)
  max_A=max(max(abs(m.Aeq)));%,[],2);
  m.Aeq = m.Aeq./max_A;
  m.beq = m.beq./max_A;
end
if(size(m.Aineq,1)>0)
  max_A=max(max(abs(m.Aineq)));%m.Aineq,[],2);
  m.Aineq = m.Aineq./max_A;
  m.bineq = m.bineq./max_A;
end

end















