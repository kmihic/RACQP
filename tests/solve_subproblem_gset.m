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

function x = solve_subproblem_gset(Q_current, c_current, x_ix, rac_model, run_p, ...
             x_curr, y_curr,beta_h,RP)
  %Implements special split for max-bisection problem 

  %\min  0.5\x^TQ\x + \c^T\x + (\beta/2)r^2
  n_vars = length(x_ix)+1;
  c_current(n_vars)=-y_curr;
  Q_current(n_vars,n_vars)=beta_h;
  model.Q = Q_current;
  model.obj = full(c_current);
  
  %s.t.   \e^T\x - r=0
  x_not=x_curr;
  x_not(x_ix)=0;
  A=rac_model.Aeq(:,x_ix);
  A(1,n_vars)=-1;
  model.A=A;
  Ax=rac_model.Aeq*x_not;
  model.rhs = rac_model.beq - Ax;
  model.sense = repmat('=', 1, 1); 

  %s.t. \x, r: binary
  lc = rac_model.local_constraints;
  model.lb = zeros(n_vars,1);
  model.ub = ones(n_vars,1);

  %define var type
  if(isfield(rac_model,'vtype'))
     model.vtype = rac_model.vtype(x_ix);
  end
  model.vtype(n_vars)='B';


  if(~isfield(run_p.gurobi_params,'TimeLimit'))
    error('Error. Time limit for gurobi sub-solver must be set')
  end
  g_out = gurobi(model,run_p.gurobi_params);
  if(~strcmpi(g_out.status, 'OPTIMAL') && ~isfield(run_p.gurobi_params,'TimeLimit')) 
    error('Error. Gurobi did not return optimal result')
  end

  if(isfield(g_out,'x'))
    x = g_out.x;
  elseif(strcmpi(g_out.status, 'TIME_LIMIT') && size(model.A,1) == 0)
    error('Time limit too short')
    g_out
    error('no solution returned by gurobi')
  end
x(n_vars)=[];
end