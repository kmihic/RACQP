
function [res_p res_d] = get_residuals(model, x, y, do_dual, Q_handle,...
                               z_current, do_z)

    %
    % bounds
    res_lb = max(0, model.lb-x);
    res_ub = max(0, x-model.ub);
    
    abs_lb = norm(res_lb,Inf);
    L2_lb = norm(res_lb)/max(1,norm(model.lb));
    max_x = norm(x,Inf);
    rel_lb = abs_lb/(1+max(norm(model.lb,Inf), max_x));

    abs_ub = norm(res_ub,Inf);
    L2_ub = norm(res_ub)/max(1,norm(model.ub));
    rel_ub = abs_ub/(1+max(norm(model.ub,Inf), max_x));

    res_p.abs_bounds = max(abs_lb,abs_ub);
    res_p.L2_bounds = max(L2_lb,L2_ub);
    res_p.rel_bounds = max(rel_lb,rel_ub);

    

%     res_b_lbub = [res_lb; res_ub];
%     b_lbub = [model.lb; model.ub];
%     b_lbub(~isfinite(b_lbub)) = 0;
%     res_p.abs_bounds = norm(res_b_lbub, Inf);
%     res_p.L2_bounds = norm(res_b_lbub)/max(1,norm(b_lbub));
%     res_p.rel_bounds = res_p.abs_bounds/(1+norm(b_lbub,Inf));

    % primal residual
    Ax_ineq = model.Aineq*x;
    Ax_eq = model.Aeq*x;
    res_eq= Ax_eq - model.beq;       
    res_ineq= max(0,Ax_ineq - model.bineq);     
    if(length(res_eq) == 0)
      res_p.abs_eq = 0;
      res_p.L2_eq = 0;
      res_p.rel_eq = 0;
      res_p.abs_ineq = norm(res_ineq, Inf);
      res_p.L2_ineq = norm(res_ineq)/max(1,norm(model.bineq));
      res_p.rel_ineq = res_p.abs_ineq/(1+max(norm(model.bineq,Inf), norm(Ax_ineq,Inf)));
    elseif(length(res_ineq) == 0)
      res_p.abs_eq = norm(res_eq, Inf);
      res_p.L2_eq = norm(res_eq)/max(1,norm(model.beq));
      res_p.rel_eq = res_p.abs_eq/(1+max(norm(model.beq,Inf), norm(Ax_eq,Inf)));
      res_p.abs_ineq = 0;
      res_p.L2_ineq = 0;
      res_p.rel_ineq = 0;
    else
      res_p.abs_eq = norm(res_eq, Inf);
      res_p.L2_eq = norm(res_eq)/max(1,norm(model.beq));
      res_p.rel_eq = res_p.abs_eq/(1+max(norm(model.beq,Inf), norm(Ax_eq,Inf)));
      res_p.abs_ineq = norm(res_ineq, Inf);
      res_p.L2_ineq = norm(res_ineq)/max(1,norm(model.bineq));    
      res_p.rel_ineq = res_p.abs_ineq/(1+max(norm(model.bineq,Inf), norm(Ax_ineq,Inf)));  
    end
    res_p.abs = max(max(res_p.abs_eq, res_p.abs_ineq), res_p.abs_bounds);
    res_p.L2 = max(max(res_p.L2_eq, res_p.L2_ineq), res_p.L2_bounds);
    res_p.rel = max(max(res_p.rel_eq, res_p.rel_ineq), res_p.rel_bounds);
    
    % dual residual
    if(~do_dual)
      res_d.abs = -inf;
      res_d.L2 = -inf;
      res_d.rel = -inf;
    else
      if(Q_handle)
        Qx = model.Q(':',':',x);
      else
        Qx = model.Q*x;
      end
      A = [model.Aeq;model.Aineq];
      % NOTE: 2Qx is bcs racqp is solving x'Qx+c'x
      %       if changed to 1/2 x'Qx (as it is in the paper), remove 2!
      Ay = A'*y;
      curr_res_d = 2*Qx + model.c - Ay ;
      if(do_z) %racqp splits x for bounds: x-x_k = 0; x_k >= 0
        curr_res_d = curr_res_d - z_current;
        res_d.abs = norm(curr_res_d, Inf);
        res_d.L2 = norm(curr_res_d)/max(1,norm(model.c));
        res_d.rel = res_d.abs/(1+max(max(norm(Qx,Inf), norm(model.c,Inf)), ...
                                     max(norm(Ay,Inf), norm(z_current,Inf))));
      else
        res_d.abs = norm(curr_res_d, Inf);
        res_d.L2 = norm(curr_res_d)/max(1,norm(model.c, Inf));
        res_d.rel = res_d.abs/(1+max(max(norm(Qx,Inf), norm(model.c,Inf)), ...
                                     norm(Ay,Inf)));
      end
    end
    
end
