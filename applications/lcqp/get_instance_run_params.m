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


function inst_param = get_instance_run_params(model, solver, max_time, binary, ...
                         epsilon, max_iter)

  inst_param = struct();
  inst_param.model = model;
  inst_param.max_time = max_time;
  inst_param.solver = solver;
  inst_param.binary = binary;
  if(nargin >=5)
    inst_param.epsilon = epsilon;
  end
  if(nargin >=6)
    inst_param.max_iter = max_iter;
  end
  inst_param.racqp_run_p = struct();
end
