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

%Default parameters when using gurobi as the sub-problem solver
function gurobi_params = default_gurobi_parameters()

  gurobi_params.presolve = 0;
  %All constraints must be satisfied to a tolerance. default: 1e-6
  %gurobi_params.FeasibilityTol=1e-4;
  %the relative difference between the primal and dual objective values
  %default: 1e-8
  %gurobi_params.BarConvTol=1e-4;
  %the improving direction in order for a model to be declared optimal. 
  %default: 1e-6
  %gurobi_params.OptimalityTol=1e-4;
  gurobi_params.outputflag = 0;
  gurobi_params.threads = 1;
  % do not change time limit!!
  % change it at application level, if needed
  gurobi_params.TimeLimit = 5;

end
