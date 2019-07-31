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
% TEST: randomly generated problems using procedure described in Section 4.2.1 
%       Mihic, K., Zhu, M. and Ye, Y., 2019. 
%       "Managing Randomization in the Multi-Block Alternating Direction 
%        Method of Multipliers for Quadratic Optimization"
%        arXiv preprint arXiv:1903.01786.
%      
%
% 
%

function [sol, model]=verify_Markowitz_min_variance_rnd(filename, solver)
addpath('../racqp');
addpath('../utils');

beta = 10;
epsilon = 1e-6;

%get the problem 
disp("LOADING THE MODEL...")
model = read_rac(filename);
if(nargin > 1 && strcmp(solver,'gurobi'))
  disp('Running Gurobi')
  sol = GUROBI(model);
  disp("Done")
  return;
else if(nargin > 1 && strcmp(solver,'mosek'))
  disp('Running Mosek')
  sol = MOSEK(model);
  disp("Done")
  return;
end
nsize=model.size;
if(nsize == 6000)
  n_blocks = 100;
elseif(nsize == 9000)
  n_blocks = 150;
else
  error('Coded for n = 6000 or 9000');
end

%load default runtime parameters to be used when solving each iteration of MIP
run_params = default_run_params;
%change some parameters
run_params.beta = beta;
run_params.epsilon = epsilon;
run_params.rnd_seed = 123;
%turn off verbose
run_params.debug = 0;
%get density
density = nnz(model.Q)/model.size^2;
if(density > 0.5)
  run_params.use_sparse = false;
else
  run_params.use_sparse = true;
end
run_params.n_blocks = n_blocks;

%print run parameters
%disp("RUN PARAMETERS")
%disp("run_params")
%disp(run_params)

%run the solver
disp("Running RACQP")
sol = RACQP(model,run_params, true);
disp("Done")

%format long
%disp("SOLUTION")
%disp(sol)

%quit
end







