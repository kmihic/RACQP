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

%TEST: Markowitz min variance problem
%      Data source Center for Research in Security Price (CRSP)
%      Time period Jan-Dec 2018
%
% Output: dataMark_var_out_regular, dataMark_var_out_lowrank
%

clear all
%Load the problem
datapath="../data/data_markowitz/";
filename=datapath+"regular_monthly";
%load_Markowitz_model(filename, p_cut, kappa, low_rank, mip)
model = load_Markowitz_model(filename, -1, 1e-5, false,false);

%test low rank model
%filename=datapath+"lowrank_monthly";
%model = load_Markowitz_model(filename, -1, 1e-5, true,false);

%using partial Lagrange to handle bounds
model.local_constraints.lb=model.lb;
model.local_constraints.ub=model.ub;
model.lb = -inf*ones(model.size,1);
model.ub = inf*ones(model.size,1);
model.local_constraints.Aineq = sparse(0,model.size);
model.local_constraints.bineq = [];
model.local_constraints.Aeq = sparse(0,model.size);
model.local_constraints.beq = [];

%load default runtime parameters
default_run_params
%change some parameters
run_params.n_blocks = 50;
run_params.epsilon_d = 1e-2;
run_params.beta = 1;
run_params.sub_solver_type='gurobi';
%turn off verbose
%run_params.debug = 0;

%load gurobi parameters
default_gurobi_parameters
run_params.gurobi_params = gurobi_params;

%run the solver
disp('Running RAC multi-block')
dataMark_var_out_regular = RACQP(model,run_params);


%single block for lowrank
%disp('Running RAC single-block')
%run_params.n_blocks = 1;
%run_params.beta = .1;
%rndMark_var_out_lowrank = RACQP(model,run_params);

disp('Solution, regular model:')
dataMark_var_out_regular

%disp('Solution, low rank model:')
%dataMark_var_out_lowrank


