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

% DEMO: randomly generated QP problem
%

function [solution, model] = demo_rnd_min_portfolio(model,n_size, sparsity_Q, kappa, ...
                   rnd_seed, n_blocks, beta, solver_mode, max_rtime, verbose, p_trial)

%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%construct a problem 
problem_construct.n_var= n_size;
problem_construct.Q_sparsity= sparsity_Q;
problem_construct.c_sparsity= 0;
problem_construct.Aeq_sparsity= 0;
problem_construct.Aineq_sparsity= 0;
problem_construct.Aineq_n_row= 0;
problem_construct.Aeq_n_row= 0;
problem_construct.rnd_seed= rnd_seed;

%To use simple Q = sprand(.) call, no bounds
if(length(model)==0)
  model = get_model_rnd_QP(problem_construct,false);
  model.Q = model.Q + kappa*speye(n_size);
  model.c = model.c + 1;
  model.Aeq = sparse(ones(1,n_size)); 
  model.lb = zeros(n_size,1);  
  model.ub = ones(n_size,1);
  model.beq = ceil(n_size/2);
  model.binary = 1:n_size;
  % partial Lagrange on the bounds
  model.local_constraints.lb=model.lb;
  model.local_constraints.ub=model.ub;
  model.local_constraints.Aineq = sparse(0,n_size);
  model.local_constraints.bineq = [];
  model.local_constraints.Aeq = sparse(0,n_size);
  model.local_constraints.beq = [];
end

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
%turn off verbose
run_params.debug = 1;

%load default gurobi parameters
default_gurobi_parameters
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 2; 
%gurobi_params.presolve=1;
%gurobi_params.outputflag=0;
run_params.gurobi_params = gurobi_params;

%set MIP run parameters
lambda = 0.4;
run_mip.permute_type = 'swap';
run_mip.permute_dist = 'exponential'; 
run_mip.permute_min = 2; 
run_mip.permute_max = n_size; 
run_mip.permute_mu = lambda*n_size; 
run_mip.max_iter = 10000000; %time limited 
run_mip.max_nperturb = 10000000; %time limited
run_mip.max_rtime = max_rtime; 
run_mip.rnd_seed = rnd_seed; 
run_mip.n_perturb_trial = max(2,round(p_trial*n_size));
run_mip.run_sub = run_params;
run_mip.debug = verbose;
run_mip.mip_epsilon = 0;



if(strcmpi('RP_ADMM', solver_mode) || strcmpi('CYCLIC_ADMM', solver_mode))
  rng(run_params.rnd_seed);
  groups={};
  g_size = model.size/run_params.n_blocks;
  ix_x = randperm(model.size);
  for ii = 0:(run_params.n_blocks-1)
    st = (ii * g_size);
    en = st + g_size;
    ix = ix_x((st+1):en);
    groups((ii+1),:) = {ix, []};
  end  
  model.groups = groups;    
  if(strcmpi('RP_ADMM', solver_mode))
    model.group_mode = 'RP';
    msg = " RP-ADMM ";
  else
    model.group_mode = 'CADMM';
    msg = " CYCLIC-ADMM ";
  end
elseif(strcmpi('RAC', solver_mode))
   msg = " RAC-ADMM ";
else
  error("Solver mode incorrect. Accepted: RAC, RP_ADMM, CYCLIC_ADMM")
end

if(n_blocks == 1)
  msg = "Solver mode ignored. Running RACQP single-block";
else
  msg = "Running RACQP multi-block:" + msg +" mode.";
end

%run the solver
fprintf("\n##############\n\n%s\n\n",msg)
solution = RACQP(model,run_mip);

end
