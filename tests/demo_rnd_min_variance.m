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

function [solution, model] = demo_rnd_min_variance(model,n_size, sparsity_Q, kappa, ...
                   rnd_seed, n_blocks, beta, solver_mode, epsilon, max_iter, max_rtime, verbose)

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
  model.beq = 1;
end

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
run_params.epsilon = epsilon;
run_params.max_iter = max_iter;
run_params.max_rtime = max_rtime;
%turn off verbose
run_params.debug = verbose;


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
solution = RACQP(model,run_params);

end
