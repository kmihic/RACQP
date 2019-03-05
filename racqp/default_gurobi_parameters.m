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

gurobi_params.presolve = 0;
gurobi_params.FeasibilityTol=1e-4;
gurobi_params.BarConvTol=1e-4;
gurobi_params.outputflag = 0;
gurobi_params.threads = 1;
gurobi_params.TimeLimit = 100;