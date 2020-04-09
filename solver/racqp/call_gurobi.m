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
function g_out = call_gurobi(m, gurobi_params)
% CALL_GUROBI Solves QP model using Gurobi
%
% Note: Gurobi solves QP models of the following format
%         min  x'H'x + c'x
%         s.t. Aeq x = beq
%              Aineq x <= bineq
%              lb <= x <= ub
%              x integer, continuous
  
   
  % setup gurobi model
  if(isfield(m,'Q') && size(m.Q,1) > 0)
    model.Q = m.Q;
  end
  model.obj = full(m.c);
  model.lb = m.lb;
  model.ub = m.ub;
  
  model.A = [m.Aeq;m.Aineq];
  model.rhs = full([m.beq;m.bineq]);
  model.sense = [repmat('=', length(m.beq), 1);repmat('<', length(m.bineq), 1)]; 

  if(isfield(m,'x_ws') && isfield(m,'y_ws'))
    model.pstart = m.x_ws;
    model.dstart = m.y_ws;
  end

  %set var type to be used by Gurobi
  vtype = repmat('C', m.size, 1);
  vtype(m.integers) = 'I';
  vtype(m.binary) = 'B';
  model.vtype = vtype;

  %add special gurobi general constraints
  if(isfield(m,'gurobi') & isfield(m.gurobi,'genconmax'))
    model.genconmax = m.gurobi.genconmax;
  end
 
  %solve
  g_out = gurobi(model,gurobi_params);
  
  if(~strcmpi(g_out.status, 'OPTIMAL') && ~isfield(gurobi_params,'TimeLimit')) 
    error('Error. Gurobi did not return optimal result')
  end
  %prepare return data
  if(isfield(g_out,'objval'))
    g_out.objval = g_out.objval + m.const;
    if(isfield(m,'norm_coeff'))
      g_out.objval = g_out.objval*m.norm_coeff;
    end
  else
    g_out.objval = Inf;
  end

end
