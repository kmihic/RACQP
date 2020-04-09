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

function sol = use_gurobi(model,r_time,epsilon, do_residuals)

  if(nargin <= 2)
   epsilon = -1;
  end
  if(nargin <= 3)
   do_residuals = false;
  end

  % gurobi does not return dual variable for variable 
  % lower/upper bounds. Converting bounds to constraints 
  % to get duals for dual residual.
  if(do_residuals)
    n = model.size;
    B=[-speye(n);speye(n)];
    b=[-model.lb;model.ub];
    ix = ~isfinite(b);
    B(ix,:)=[];
    b(ix)=[];
    model.Aineq = [model.Aineq;B];
    model.bineq = [model.bineq;b];
    model.lb = ones(n,1)*(-Inf);  
    model.ub = ones(n,1)*(Inf);
  end

  gp = struct();
  gp.threads = 1;
  gp.presolve = 1;
  gp.outputflag = 1;
  if(r_time>0)
   gp.TimeLimit = r_time;
  end
  if(epsilon > 0)  
    %primal feasibility
    gp.FeasibilityTol=epsilon;
    %primal/dual gap 
    gp.BarConvTol=epsilon;
    % dual feasibility 
    gp.OptimalityTol=epsilon;
  end

  g_out = call_gurobi(model,gp);
  % prepare the output to be read by run*test scripts
  if(strcmpi(g_out.status, 'OPTIMAL'))
    sol.sol_x = g_out.x;
    % get residuals
    if(do_residuals)
      sol.sol_y = g_out.pi;  
      [res_p res_d] = get_residuals(model, g_out.x, g_out.pi, true, ...
                          false, [], false);
      sol.sol_res_p = res_p;
      sol.sol_res_d = res_d;
    else
      sol.sol_y = [];
      sol.sol_res_p.abs = 0;
      sol.sol_res_d.abs = 0;
      sol.sol_res_p.rel = 0;
      sol.sol_res_d.rel = 0;
      sol.sol_res_p.L2 = 0;
      sol.sol_res_d.L2 = 0;
   end
  else
    sol.sol_x = [];
    sol.sol_y = [];
    sol.sol_obj_val = inf;
    sol.sol_res_p.abs = -inf;
    sol.sol_res_d.abs = -inf;
    sol.sol_res_p.rel = -inf;
    sol.sol_res_d.rel = -inf;
    sol.sol_res_p.L2 = -inf;
    sol.sol_res_d.L2 = -inf;
  end
  sol.n_iter = max(g_out.itercount, g_out.baritercount);
  sol.runtime = g_out.runtime;
  if(isfield(g_out,'objval'))
    sol.sol_obj_val = g_out.objval;
  else
    sol.sol_obj_val = -inf;
  end

end


