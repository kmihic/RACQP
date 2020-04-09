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



function rac_out = rac_single_block(model, run_p, time_start)
% RAC_SINGLE_BLOCK Solves QP using RAC single x block approach 
  
  %start the clock
  stime_start = time_start;   

  % Remove inequalities
  n_var_x_orig = model.size;
  if(length(model.bineq > 0) && run_p.single_block_embed_slacks && strcmpi('diag_kkt',run_p.single_block_solver))
    m1 = size(model.Aeq,1);
    m2 = size(model.Aineq,1);  
    model.Q = [model.Q, sparse(n_var_x_orig,m2); sparse(m2,n_var_x_orig+m2)];
    model.c = [model.c; sparse(m2,1)];
    model.Aeq = [model.Aeq, sparse(m1,m2); model.Aineq, speye(m2,m2)];
    model.Aineq = sparse(0,n_var_x_orig+m2);
    model.beq = [model.beq;model.bineq];
    model.bineq = [];
    model.x0(n_var_x_orig+m2) = 0;
    model.lb(n_var_x_orig+m2) = 0;
    model.ub(n_var_x_orig+1:n_var_x_orig+m2) = inf;
    model.size = n_var_x_orig+m2;
  end

  
  n_var_x = model.size;
  m1 = size(model.Aeq,1);
  m2 = size(model.Aineq,1);  
  m = m1+m2;
  beta = run_p.beta;
  beta_h = beta/2;  
  beta_i = 1/beta;
  beta_s = sqrt(beta);
  beta_si = 1/beta_s;

   % constants and defaults for residuals
  b_ineq_max = norm(model.bineq,Inf);
  b_eq_max = norm(model.beq,Inf);
  c_max = norm(model.c,Inf);
  d_res_z = 0;
  d_res_A = 0;
  p_res_eq = 0;
  p_res_ineq = 0;      
  p_res_z = 0;  
  use_dual_res = run_p.use_dual_res;
  do_rel_tol = run_p.do_rel_tol;

  % check if the problem is bounded
  n_inf = sum(isinf(model.lb)) + sum(isinf(model.ub));
  if(n_inf == 2*n_var_x)
    bounds_present = false;
  else
    bounds_present = true;
  end
  %Check what numerical approach we are using
  use_ldl = false;
  use_kkt = false;
  use_gurobi = false;
  do_dual_bounds = false;
  if(isfield(run_p,'single_block_solver'))
    if(strcmpi('cholesky',run_p.single_block_solver))
       use_ldl = true;
    elseif(strcmpi('diag_kkt',run_p.single_block_solver))
       if(length(model.bineq) > 0)
         error("diag_kkt mode can not be used with ineq constraints")
       end
       if(~isdiag(model.Q))
         error("diag_kkt mode can not be used with non diagonal H")
       end
       use_kkt = true;
    elseif(strcmpi('bounds_by_gurobi',run_p.single_block_solver))
       % bounds handled by gurobi, other constraints normal way
       use_gurobi = true;
       if(bounds_present)
         do_dual_bounds = true;
         % gurobi does not return dual variable for variable 
         % lower/upper bounds. Converting bounds to constraints 
         % to get duals for dual residual.
         
         lb_ix = isfinite(model.lb);
         n_lb = length(find(lb_ix));
         ub_ix = isfinite(model.ub);
         n_ub = length(find(ub_ix));

         if(n_lb+n_ub > 0)
           %model.g_A = [-speye(n_lb);speye(n_ub)];
           Ab = sparse(n_lb+n_ub,model.size);
           if(n_lb > 0)
             a_ix = sub2ind(size(Ab),1:n_lb, find(lb_ix)');
             Ab(a_ix) = -1;
           end
           if(n_ub > 0)
             a_ix = sub2ind(size(Ab),(1+n_lb):(n_lb+n_ub), find(ub_ix)');
             Ab(a_ix) = 1;
           end
           model.g_A = Ab;
           model.g_rhs = [-model.lb(lb_ix);model.ub(ub_ix)];         
           model.g_sense = [repmat('<', length(model.g_rhs), 1)]; 
         else
           model.g_A = sparse(0,model.size);
           model.g_rhs = [];
           model.g_sense = [];
         end
         model.g_lb = ones(model.size,1)*(-Inf);  
         model.g_ub = ones(model.size,1)*(Inf);
         model.n_lb = n_lb;
         model.n_ub = n_ub;
         model.lb_ix = lb_ix;
         model.ub_ix = ub_ix;
       end
     end
 else
    use_ldl = true;
 end

  % use split technique for bounds if defined and not handled by partial Lagrange
  if(use_gurobi || ~bounds_present)
    do_xk = 0;
    Q = model.Q;
    z_current = [];
  else
    do_xk = 1;
    z_current = zeros(n_var_x,1);  
    xk = model.x0;    
    if(use_kkt) % testing 1/2xQx model
      Q = model.Q+beta_h*speye(n_var_x);
    else
      Q = model.Q+beta_h*speye(n_var_x);
    end
  end

  %Define auxiliary vectors and matrice
  if(~run_p.use_sparse)
    Q=full(Q);
    Aeq = full(model.Aeq);
    Aineq = full(model.Aineq);
    mc = full(model.c);
    mQ = full(model.Q);
    sp_m1 = zeros(m1,1);
  else
    Aeq = model.Aeq;
    Aineq = model.Aineq;
    mc = model.c;
    mQ = model.Q;
    sp_m1 = sparse(m1,1);
  end
  A = [Aeq; Aineq];
  A_T = A';
  b = [model.beq;model.bineq];

  
  x_current = model.x0;
  y_current = zeros(m,1);
  s_current = max(model.bineq-Aineq*model.x0,0);
  Ares = [Aeq;Aineq];
  

  %initialize counters and vars
  curr_iter = 0;
  sub_model_time = 0;
  sub_solver_time = 0;
  res_iter=[];
  res_iter_dual=[];
  curr_res.prim = Inf;
  curr_res.dual = Inf;

  %prepare data arrays 
  if(size(A,1) == 0 || use_kkt)
    c_stat = mc;
  else
    c_stat =  mc-beta*A'*b;
  end

  if(use_gurobi) %gurobi used to handle bounds
    Q_current = Q+beta_h*A'*A;
  elseif(use_kkt) %equality constrained problems only
    Qinv = inv(Q);
    Qinv_m = -Qinv;
    AQinv = Aeq * Qinv;
    AQinv_m = -AQinv;
    AQinv_mt = AQinv_m';
    Q_fact = AQinv*Aeq';
    % testing: Q_fact = AQinv*Aeq' + run_p.single_block_kkd_diag_fact*speye(m);

    R = chol(Q_fact);

  else %default ldl factorization,     
    sbA=beta_s*A;
   % using speye, for some reasons it wors better than is using 
   % eye(m) regardless of Q,A sparsity
     [L,D,P] = ldl([2*Q, sbA';sbA, -speye(m)]);
  end

  init_time = toc(stime_start);
  %main loop
  while(~terminate(run_p, curr_iter, curr_res, time_start,false))

    %prepare rhs
    stime_start = tic; 
    %have bounds (and Gurobi not used)
    if(do_xk == 1) 
      c_current = c_stat- beta*xk - z_current;
    else
      c_current = c_stat;
    end
    %dual vars
    if(~use_kkt && size(A,1)>0) 
      c_current = c_current - A'*y_current;
    end
    %add slacks if used (will be skipped if kkt, as no ineq)
    if(length(model.bineq) > 0) 
      c_current = c_current + beta*Aineq'*s_current;
    end     
    sub_model_time = toc(stime_start) + sub_model_time;
    
    %solve sub-model
    stime_start = tic;
    if(use_gurobi)
      [result_x, z_current] = solve_subproblem_gurobi(Q_current, c_current,...
                             model, run_p, do_dual_bounds);
    elseif(use_kkt)
      q = AQinv_m * c_current -2*model.beq;   
 %     y_kkt =  P*(L'\(D\(L\(P'*(q)))));

      y_kkt = R\(R'\q);
      result_x = (Qinv_m*c_current + AQinv_mt*y_kkt)/2; 
      if(~do_xk)
        %no bounds, nothing else to do
        y_current = -y_kkt;
        x_current=result_x(1:n_var_x);        
        curr_iter = curr_iter + 1;
        break;
      end
    else %default ldl factorization        
      q = [-c_current;zeros(m,1)];
      result_x = P*(L'\(D\(L\(P'*(q)))));
    end
    sub_solver_time = toc(stime_start) + sub_solver_time;
    %update x  
    x_current=result_x(1:n_var_x);
    
    if(~use_kkt) 
      %update Ax (y variables)
      if(use_gurobi)
        Ax = A*x_current;
      else
        Ax = (beta_si)*result_x(n_var_x+1:n_var_x+m);
      end
      %update s if needed
      if(length(model.bineq) > 0)
        s_current = max(0,beta_i*y_current(m1+1:m)+model.bineq-Ax(m1+1:m));
        res_y = (Ax + [sp_m1;s_current] - b);
      elseif(length(model.beq) > 0)
        res_y = (Ax - b);
      else 
        res_y = 0;
      end
     % curr_res = max(abs(res_y)); 
      %update dual y
      y_current = y_current - run_p.stepY*beta*res_y;  
    else
      % kkt mode has no ineq constraints
      y_current = -y_kkt;
      Ax = A*x_current;
      if(length(model.beq) > 0)
        res_y = (Ax - b);
      else 
        res_y = 0;
      end      
   %   curr_res = full(max(abs(res_y))); %matlab returns a sparse scalar !?! 
    end
    
   % get residuals -primal    
    if(length(model.bineq) > 0)
      res_y_ineq = norm(res_y(m1+1:m),Inf); 
      if(do_rel_tol)
        p_res_ineq = res_y_ineq / (1+max(b_ineq_max, ...
                               norm(Ax(m1+1:m),Inf)));
      else
        p_res_ineq = res_y_ineq;
      end
    end 
    if(length(model.beq) > 0)
      res_y_eq = norm(res_y(1:m1),Inf);
      if(do_rel_tol)
        p_res_eq = res_y_eq / (1+max(b_eq_max, ...
                               norm(Ax(1:m1),Inf)));
      else
        p_res_eq = res_y_eq;
      end
    end

    %update auxiliary x and dual z if needed
    if(do_xk)
      xk = min(max(model.lb,x_current-beta_i*z_current),model.ub);
      diff_z = (x_current-xk);
      z_current = z_current - run_p.stepZ*beta*diff_z;
     % res_z = norm(diff_z,Inf);
      if(do_rel_tol)
       res_lb = max(0, model.lb-x_current(1:n_var_x));
       res_ub = max(0, x_current(1:n_var_x)-model.ub);
       abs_lb = norm(res_lb,Inf);
       max_x = norm(x_current(1:n_var_x),Inf);
       rel_lb = abs_lb/(1+max(norm(model.lb,Inf), max_x));
       abs_ub = norm(res_ub,Inf);
       rel_ub = abs_ub/(1+max(norm(model.ub,Inf), max_x));
       p_res_z = max(rel_lb,rel_ub);
%        p_res_z = res_z / (1+max(norm(x_current(1:n_var_x), Inf),...
%                             norm(xk, Inf)));  
      else
        p_res_z = norm(diff_z,Inf);%res_z;
      end
    end
    if(use_dual_res && (do_xk || (use_gurobi && do_dual_bounds)))
      d_res_z = norm(z_current,Inf);
    end
    % calc the overall residuals
    % primal
    curr_res.prim = max(p_res_z,max(p_res_eq, p_res_ineq));
    % dual
    if(use_dual_res)
      Ay = A_T*y_current;
      Qx = model.Q*x_current;
      curr_res_d = 2*Qx(1:n_var_x) + model.c - Ay ;
      if(do_xk || (use_gurobi && do_dual_bounds)) 
        curr_res_d = curr_res_d - z_current;
      end
      d_res_abs = norm(curr_res_d, Inf);
      if(do_rel_tol)
        curr_res.dual = d_res_abs/(1+max(max(norm(Qx,Inf), c_max), ...
                                     max(norm(Ay,Inf), d_res_z)));
      else
         curr_res.dual = d_res_abs;
      end
     else
       curr_res.dual = 0;
     end
     
    %finalize the iteration    
    curr_iter = curr_iter + 1;

    if(isfield(run_p,'debug') && run_p.debug > 0)
      x = x_current(1:n_var_x);
      obj_val = x'*mQ*x+mc'*x+model.const;  
      if(run_p.debug > 0)
        s=sprintf("single_block (%d): %.3e %.3e %.3e  (residual prim, dual, obj val)",...
            curr_iter, curr_res.prim, curr_res.dual, obj_val);
        disp(s)
      end
      %store for output  
      res_iter(curr_iter) = curr_res.prim;
      res_iter_dual(curr_iter) = curr_res.dual;
    end
%%% END WHILE
  end %main loop
  rac_time = toc(time_start);

  %prepare the output struct  
  x = x_current(1:n_var_x);
  rac_out.sol_x = x_current(1:n_var_x_orig);
  rac_out.sol_y = y_current;
  %gurobi used to handle bounds
  if(do_dual_bounds)
     do_xk = true;
  end 
  if(do_xk)
    rac_out.sol_z = z_current(1:n_var_x_orig);
  else
    rac_out.sol_z = [];
  end
  obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
  rac_out.sol_obj_val = obj_val;
  rac_out.n_iter = curr_iter;

  % get residuals
 [res_p res_d ] = get_residuals(model, x, y_current, true, ...
                         false, z_current, do_xk);
 rac_out.sol_res_p = res_p;
 rac_out.sol_res_d = res_d;

 % finalize the output  
 rac_out.runtime = rac_time;  
 rac_out.init_time = init_time;
 rac_out.solver_time = sub_solver_time;
 rac_out.model_time = sub_model_time;
 if(isfield(run_p,'debug') && run_p.debug > 0)
   rac_out.res_iter = res_iter;
   rac_out.res_iter_dual = res_iter_dual;
 end
end %function solve

function [x,z] = solve_subproblem_gurobi(Q_current, c_current, ...
    model_rac, run_p, do_dual_bounds)

  model.Q = Q_current;
  model.obj = full(c_current);
 
  if(do_dual_bounds && length(model_rac.g_rhs > 0))    
    model.A = model_rac.g_A;
    model.rhs = model_rac.g_rhs;
    model.sense = model_rac.g_sense;  
    model.lb = model_rac.g_lb;  
    model.ub = model_rac.g_ub;
    run_p.gurobi_params.presolve = 1;
  else
    % no local constraints
    model.A= sparse(0,length(c_current));
    model.rhs = [];
    model.sense = repmat('=', length(model.rhs), 1);
    %add local bounds
    model.lb = model_rac.lb;
    model.ub = model_rac.ub;
  end


%  run_p.gurobi_params.presolve = 1;
 % run_p.gurobi_params.outputflag = 1;
  % solve the problem
  if(~isfield(run_p.gurobi_params,'TimeLimit'))
    error('Error. Time limit for gurobi sub-solver must be set')
  end
  g_out = gurobi(model,run_p.gurobi_params);

  %g_out = gurobi(model,gp);

  % get the output
  if(~strcmpi(g_out.status, 'OPTIMAL') && ~isfield(run_p.gurobi_params,'TimeLimit')) 
    error('Error. Gurobi did not return optimal result')
  end

  if(isfield(g_out,'x'))
    x = g_out.x;
    if(do_dual_bounds)
      z_lb = zeros(model_rac.size,1);
      z_ub = zeros(model_rac.size,1);
      n_lb = model_rac.n_lb;
      n_ub = model_rac.n_ub;
      lb_ix = model_rac.lb_ix;
      ub_ix = model_rac.ub_ix;
      if(n_lb > 0 && n_ub > 0)
        z_lb(lb_ix) = g_out.pi(1:n_lb);
        z_ub(ub_ix) = g_out.pi(1+n_lb:n_lb+n_ub);
      elseif(n_lb > 0)        
        z_lb(lb_ix) = g_out.pi(1:n_lb);
      elseif(n_ub > 0)
        z_ub(ub_ix) = g_out.pi(1:n_ub);
      end
      z = -z_lb+z_ub;
    else
      z = [];
    end     
  elseif(strcmpi(g_out.status, 'TIME_LIMIT') && size(model.A,1) == 0)
    error('Time limit too short')
  else
    g_out
    error('No solution returned by gurobi')
  end
  
end


