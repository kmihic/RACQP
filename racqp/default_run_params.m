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

run_params.n_threads =  1;
run_params.n_blocks =  10;
run_params.beta =  50;
run_params.sub_solver_type =  'cholesky';
run_params.max_iter =  500;
run_params.epsilon =  1e-4;
run_params.max_rtime =  600;
run_params.rnd_seed =  1234;
run_params.use_sparse = false;
run_params.debug = 1;