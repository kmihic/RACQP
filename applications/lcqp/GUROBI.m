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

function sol = GUROBI(model,r_time)

  gp = default_gurobi_parameters; 
  if(length(model.binary) > 0 || length(model.integers) > 0)
    gp.presolve = 1;
  else
    gp.presolve = 0;
  end
  gp.outputflag = 1;
  if(nargin == 2)
   gp.TimeLimit = r_time;
  else
   gp.TimeLimit = inf;
  end

  g_out = call_gurobi(model,gp);
  % prepare the output to be read by run*test scripts
  if(isfield(g_out,'objval'))
    sol.sol_obj_val = g_out.objval;
  else
    sol.sol_obj_val = inf;
  end
  sol.rac_time = g_out.runtime;
  sol.n_iter = -inf;
  sol.sol_residue_p = -inf;
  sol.sol_residue_d = -inf;

end


