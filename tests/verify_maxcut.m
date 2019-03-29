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

% Verify results: Max-cut problem
%      Using GSET instances 
%      (http://web.stanford.edu/ yyye/yyye/Gset)
%


function [gset_maxcut, model] = verify_maxcut(filename, max_time)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%hardcode experimental run-parameters
n_blocks = 4;
p_trial = 0.005;
lambda = 0.4;

%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_GSET(filename);

nsize=model.size;
%load default runtime parameters to be used when solving each iteration of MIP
default_run_params
%change some parameters
run_params.n_blocks = n_blocks;
run_params.rnd_seed = 123;
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
run_mip.permute_mu = lambda*nsize; 
run_mip.max_iter = 10000000; %time limited 
run_mip.max_nperturb = 10000000; %time limited
run_mip.max_rtime = max_time; 
run_mip.rnd_seed = 123; 
run_mip.n_perturb_trial = max(2,round(p_trial*nsize));
run_mip.run_sub = run_params;
run_mip.debug = 1;
run_mip.mip_epsilon = 0;

%print run parameters
disp("RUN PARAMETERS")
disp("run_mip")
disp(run_mip)
disp("run_sub")
disp(run_mip.run_sub)
disp("gurobi_parameters")
disp(run_mip.run_sub.gurobi_params)

%run the solver
disp('Running RAC MIP')
gset_maxcut = RACQP(model,run_mip);

disp("SOLUTION")
disp(gset_maxcut)

%quit

end





