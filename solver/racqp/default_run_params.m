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
%Define default runtime parameters

function run_params = default_run_params(n_blocks,beta,max_time,...
                       epsilon,max_iter)
  
  if(nargin <= 2)
    max_time = Inf;
  end
  if(nargin <= 3) 
    epsilon = 1e-5;
  end
  if(nargin <= 4)
    max_iter = 4000;
  end


  run_params.n_threads =  1;
  run_params.n_blocks =  n_blocks;
  run_params.beta =  beta;
  run_params.sub_solver_type =  'cholesky';
  run_params.single_block_solver = 'cholesky';
  run_params.max_iter = max_iter;
  run_params.min_iter = 0;
  run_params.max_rtime =  max_time;
  run_params.rnd_seed =  123;
  run_params.use_sparse = true;
  run_params.debug = 0;
  run_params.single_block_embed_slacks = true;
  % recent updates - dual and relative residuals  
  run_params.use_dual_res = true;
  run_params.do_rel_tol = true;  
  run_params.epsilon_prim =  epsilon;
  run_params.epsilon_dual =  epsilon;
  run_params.calc_dual_res = true;
  % experimental  
  run_params.clean_embedd = false;  %experimental
  run_params.stepZ = 1;
  run_params.stepY = 1;
  % if using single block, kkt mode, in some
  % cases (e.g. LP problems where Q=[] or actual QP but..), 
  % we may get semidefinite or even small indefinite
  % Fixed by adding a small diagonal values
  run_params.single_block_kkd_diag_fact = 0; 
end
