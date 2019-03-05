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

%TEST: Max-cut problem
%      Using GSET instances 
%      (http://web.stanford.edu/ yyye/yyye/Gset)
%
% Output: gset_maxcut
%

clear all
%load and convert the problem to rac style
datapath="../data/data_gset/";
filename=datapath+"G55";

max_time = 30;
n_blocks = 4;
beta = 0.005;
p_trial = 0.001;
lambda = 0.2;

m_gset_mc = get_GSET(filename);

nsize=m_gset_mc.size;
%load default runtime parameters to be used when solving each iteration of MIP
default_run_params
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
run_params.epsilon = 1e-2;
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
gset_maxcut = RACQP(m_gset_mc,run_mip);

disp("SOLUTION")
disp(gset_maxcut)





