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
%TEST: Markowitz portfolio selection problem
%      Data source Center for Research in Security Price (CRSP)
%      Time period Jan-Dec 2018
%
% Output: dataMark_bin
%

function [dataMark_bin, model]=verify_Markowitz_portfolio_selection_data(filename, low_rank_s, max_time)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../tests');

if(strcmpi("true",low_rank_s))
  low_rank = true;
else
  low_rank = false;
end

%setup problem data
kappa = 1e-5;
p_cut = 0.5;
if(low_rank)
  n_blocks = 50;
  beta = 0.5;
else
  n_blocks = 100;
  beta = 0.05;
end
lambda=0.4;
p_trial = 0.005;

%get the problem 
disp("LOADING THE MODEL...")
model = load_Markowitz_model(filename, p_cut, kappa, low_rank, true);
nsize=model.size;
model.local_constraints.lb=model.lb;
model.local_constraints.ub=model.ub;
model.local_constraints.Aineq = sparse(0,nsize);
model.local_constraints.bineq = [];
model.local_constraints.Aeq = sparse(0,nsize);
model.local_constraints.beq = [];

%load default runtime parameters to be used when solving each iteration of MIP
default_run_params
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
run_params.rnd_seed = 123;
%turn off verbose
run_params.debug = 1;
%Check density - no point, time is driven by Gurobi which uses sparse Q and A


%load default gurobi parameters
default_gurobi_parameters
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 2;
 
%gurobi_params.presolve=1;
%gurobi_params.outputflag=0;
run_params.gurobi_params = gurobi_params;

%set MIP run parameters
run_mip.permute_type = 'swap';
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
dataMark_bin = RACQP(model,run_mip);

format long
disp("SOLUTION")
disp(dataMark_bin)
if(p_cut > 0)
  disp("SUM OF X SOL ")
  disp(nnz(round(dataMark_bin.sol_x(model.binary))))
end

%quit
end







