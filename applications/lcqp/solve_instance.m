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



function sol = solve_instance(ip)

  % call the solver
  if(strcmp(ip.solver,'gurobi'))
    disp('Running Gurobi')
    if(ip.binary) %no residuals, no epsilon
      sol = use_gurobi(ip.model,ip.max_time, -1, false);
    else
      sol = use_gurobi(ip.model,ip.max_time,ip.epsilon,true);
    end
  elseif(strcmp(ip.solver,'mosek'))
    disp('Running Mosek')
    if(ip.binary)
     disp("Mosek not set to run binary")
     quit
    end
    sol = use_mosek(ip.model,ip.max_time,ip.epsilon);
  elseif(strcmp(ip.solver,'osqp'))
    disp('Running OSQP')
    if(ip.binary)
      disp("OSQP can not solve binary problems")
      quit
    end
    sol = use_osqp(ip.model,ip.max_time, ip.epsilon, ip.max_iter);
  elseif(strcmp(ip.solver,'racqp'))
    disp('Running RACQP')
    sol = RACQP(ip.model, ip.racqp_run_p, true);
  end
  disp("Done");
end
