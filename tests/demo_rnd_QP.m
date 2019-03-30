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

function [solution, model] = demo_rnd_QP(model,n_size, m_eq, m_ineq, sparsity_Q, sparsity_A, ...
                   rnd_seed, n_blocks, beta, solver_mode, epsilon, max_iter, max_rtime, verbose)

addpath('../racqp');
addpath('../utils');
%construct a problem 
problem_construct.n_var= n_size;
problem_construct.Q_sparsity= sparsity_Q;
problem_construct.c_sparsity= 0;
problem_construct.Aeq_sparsity= sparsity_A;
problem_construct.Aineq_sparsity= sparsity_A;
problem_construct.Aineq_n_row= m_ineq;
problem_construct.Aeq_n_row= m_eq;
problem_construct.rnd_seed= rnd_seed;

%To use simple Q = sprand(.) call, no bounds
if(length(model)==0)
  model = get_model_rnd_QP(problem_construct,true);
end
%
% to use a controlled way to construct Q call, with bounds
% can be too slow for large n_size
% Q_construct.cond_num= 5;
% Q_construct.zeta= 2;
% Q_construct.eta= 0.5000;
% model = get_model_rnd_QP(problem_construct,true,Q_construct);

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
