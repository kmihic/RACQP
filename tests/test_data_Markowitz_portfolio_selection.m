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

clear all
%Load the problem
%regular
datapath="../data/data_markowitz/";
filename=datapath+"regular_daily";
low_rank = false;

%low rank
%datapath="../data/data_markowitz/";
%filename=datapath+"lowrank_daily";
%low_rank = true;

%setup problem data
kappa = 1e-5;
p_cut = 0.5;
if(low_rank)
  n_blocks = 10;
else
  n_blocks = 50;
end
beta = 1;
lambda=0.2;
max_time = 60;
max_iter = 10;
p_trial = 0.001;

%get the problem 
%load_Markowitz_model(filename, p_cut, kappa, low_rank, mip)
model = load_Markowitz_model(filename, p_cut, kappa, low_rank, true);
N=model.size;
model.local_constraints.lb=model.lb;
model.local_constraints.ub=model.ub;
model.local_constraints.Aineq = sparse(0,N);
model.local_constraints.bineq = [];
model.local_constraints.Aeq = sparse(0,N);
model.local_constraints.beq = [];

%load default runtime parameters to be used when solving each iteration of MIP
default_run_params
run_params.n_blocks = n_blocks;
run_params.beta = beta;
%must meet p_cut, (for low rank others are actually min error on x, wheat we minimize for)
run_params.epsilon = 1;

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

%run_mip.permute_type = 'permute';
run_mip.permute_type = 'swap';
run_mip.permute_dist = 'exponential'; 
run_mip.permute_min = 2;
run_mip.permute_max = N;
run_mip.permute_mu = lambda*N;
if(low_rank) %n_uter limited
  run_mip.max_iter = max_iter; 
  run_mip.max_rtime = 100000;
else %time limited
  run_mip.max_iter = 1000000; 
  run_mip.max_rtime = max_time;
end

run_mip.max_nperturb = 10000000; %time limited
run_mip.rnd_seed = 123; 
run_mip.n_perturb_trial = max(2,round(p_trial*N));
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
dataMark_bin = RACQP(model,run_mip);

disp("SOLUTION")
disp(dataMark_bin)
if(p_cut > 0)
  disp("SUM OF X SOL ")
  disp(nnz(round(dataMark_bin.sol_x(model.binary))))
end







