%
% RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
% Copyright (C) 2019 
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>
%
% This file is part of RACQP 
%


% Verify results: QAP relaxed (x\in[0,1])
%      Using QAPLIB benchmark problems 
%      (http://anjos.mgi.polymtl.ca/qaplib/)
%
% Output: qapcont_out


function [qapcont_out, model]=verify_QAPLIB_cont(filename) 
addpath('../racqp');
addpath('../utils');


%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_QAPLIB(filename);
nsize = model.size;
%num blocks = r, block size = r => RP mode
model.group_mode = 'RP';
%This is the relaxed QAP, keep x within bounds, but remove binary constraint
model.binary=[];

%load default runtime parameters and change them
default_run_params
r = floor(sqrt(model.size));
run_params.n_blocks = r;
run_params.beta = r;

run_params.epsilon_d = -Inf;
run_params.max_iter = 1000;
%Check density
density_Q=nnz(model.Q)/size(model.Q,1)^2;
%model.Aeq is sparse, no need to check
if(density_Q > 0.5)
  run_params.use_sparse = false;
else
  run_params.use_sparse = true;
end

%turn off verbose
run_params.debug = 0;

%print run parameters
disp("RUN PARAMETERS")
disp("run_params")
disp(run_params)

%run the solver
disp('Running RAC')
qapcont_out = RACQP(model, run_params);

disp('Solution, binary QAP:')
disp(qapcont_out)

end






