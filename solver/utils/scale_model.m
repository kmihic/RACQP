

function [m,max_Obj, max_A_eq, max_A_ineq] = scale_model(m, scale_rows, scale_objective, use_max_Q_c)

  if(nargin <= 3)
   use_max_Q_c = false;
  end
  max_Obj = 1;
  max_A_eq = [];
  max_A_ineq = [];

  if(scale_objective)
    % max returns a scalar as sparse one cell field ?!?
    max_Q = full(max(abs(m.Q),[],'all'));
    if(max_Q == 0 || use_max_Q_c)
      max_c = full(max(abs(m.c)));
    else
      max_c = 0;
    end
    max_Qc = max(max_Q, max_c);
    if(max_Qc > 0)
      max_Obj = max_Qc;
      m.Q = m.Q./max_Obj;
      m.c = m.c./max_Obj;
    end
  end

  if(~scale_rows)
    return;
  end
  
  if(length(m.beq>0))
    max_A_eq = max(abs(m.Aeq),[],2);
    m.Aeq = m.Aeq./max_A_eq;
    m.beq = m.beq./max_A_eq;
  end

  if(length(m.bineq>0))
    max_A_ineq = max(abs(m.Aineq),[],2);
    m.Aineq = m.Aineq./max_A_ineq;
    m.bineq = m.bineq./max_A_ineq;
  end
end
