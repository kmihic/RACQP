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

% Verify results: Max-bisection problems
%      Using GSET instances 
%      (http://web.stanford.edu/ yyye/yyye/Gset)
%
%

function [sol, model]= verify_maxcut(filename, max_time, solver)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%hardcode experimental run-parameters
n_blocks = 4;
p_trial = 0.005;
lambda = 0.4;
beta = 0.005;

%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_GSET(filename,123,0.5);
if(nargin > 2 && strcmp(solver,'gurobi'))
  disp('Running Gurobi')
  sol = GUROBI(model,max_time);
  disp("Done")
  return;
end
nsize=model.size;
%load default runtime parameters to be used when solving each iteration of MIP
run_params = default_run_params;
%change some parameters
run_params.n_blocks = n_blocks;
run_params.beta = beta;
run_params.epsilon = 0;
run_params.rnd_seed = 123;
run_params.sub_solver_type='user_defined';
run_params.sub_solver_f = @solve_subproblem_gset;
%turn off verbose
run_params.debug = 0;

%load default gurobi parameters
gurobi_params = default_gurobi_parameters;
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 5;
 
%gurobi_params.presolve=1;
%gurobi_params.outputflag=1;
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
run_mip.debug = 0;
run_mip.mip_epsilon = 0;

%print run parameters
% disp("RUN PARAMETERS")
% disp("run_mip")
% disp(run_mip)
% disp("run_sub")
% disp(run_mip.run_sub)
% disp("gurobi_parameters")
% disp(run_mip.run_sub.gurobi_params)

%run the solver
disp('Running RACQP MIP')
sol = RACQP(model,run_mip, true);
disp("Done")

% disp("SOLUTION")
% disp(sol)
% disp("SUM OF SOL")
% disp(nnz(round(gset_maxbisect.sol_x(1:nsize))))

%quit
end




