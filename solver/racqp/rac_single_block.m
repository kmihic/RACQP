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
  if(length(model.bineq > 0) && run_p.single_block_embed_slacks)
    n_var_x = model.size;
    m1 = size(model.Aeq,1);
    m2 = size(model.Aineq,1);  
    model.Q = [model.Q, sparse(n_var_x,m2); sparse(m2,n_var_x+m2)];
    model.c = [model.c; sparse(m2,1)];
    model.Aeq = [model.Aeq, sparse(m1,m2); model.Aineq, speye(m2,m2)];
    model.Aineq = sparse(0,n_var_x+m2);
    model.beq = [model.beq;model.bineq];
    model.bineq = [];
    model.x0(n_var_x+m2) = 0;
    model.lb(n_var_x+m2) = 0;
    model.ub(n_var_x+1:n_var_x+m2) = inf;
    model.size = n_var_x+m2;
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

  


  %Check what numerical approach we are using
  use_ldl = false;
  use_kkt = false;
  use_gurobi = false;
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
    end
 else
    use_ldl = true;
 end

  % use split technique for bounds if defined and not handled by partial Lagrange
  n_inf = sum(isinf(model.lb)) + sum(isinf(model.ub));
  if(use_gurobi || n_inf == 2*n_var_x)
    do_xk = 0;
    Q = model.Q;
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
  curr_res = Inf;

  %prepare data arrays 
  if(size(A,1) == 0 || use_kkt)
    c_stat = mc;
  else
    c_stat =  mc-beta*A'*b;
  end

  if(use_gurobi) %gurobi used to handle bounds
    Q_current = Q+beta_h*A'*A;
    lb = model.lb;    
    ub = model.ub;
  elseif(use_kkt) %equality constrained problems only
    Qinv = inv(Q);
    Qinv_m = -Qinv;
    AQinv = Aeq * Qinv;
    AQinv_m = -AQinv;
    AQinv_mt = AQinv_m';
    Q_fact = AQinv*Aeq';
    
    R = chol(Q_fact);

  else %default ldl factorization,     
    sbA=beta_s*A;
   % using speye, for some reasons it wors better than is using 
   % eye(m) regardless of Q,A sparsity
     [L,D,P] = ldl([2*Q, sbA';sbA, -speye(m)]);
  end

  init_time = toc(stime_start);
  %main loop
  while(~terminate(run_p, curr_iter, curr_res, time_start))

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
      %q = [c_current;zeros(m,1)];
      result_x = solve_subproblem_gurobi(Q_current, c_current, model.lb, model.ub, run_p);
    elseif(use_kkt)
      q = AQinv_m * c_current -2*model.beq;   
      y_kkt = R\(R'\q);
      result_x = (Qinv_m*c_current + AQinv_mt*y_kkt)/2; 
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
      curr_res = max(abs(res_y)); 
      %update dual y
      y_current = y_current - beta*res_y;  
    else
      % kkt mode has no ineq constraints and eq constraints are met
      curr_res = 0;
    end
    
    %update auxiliary x and dual z if needed
    if(do_xk)
      xk = min(max(model.lb,x_current-beta_i*z_current),model.ub);
      res_z = (x_current-xk);
      z_current = z_current - beta*res_z;
      curr_res = max(curr_res, max(abs(res_z)));
    end
     
    %finalize the iteration    
    curr_iter = curr_iter + 1;

    if(isfield(run_p,'debug') && run_p.debug > 0)
      x = x_current(1:n_var_x);
      obj_val = x'*mQ*x+mc'*x+model.const;  
      if(run_p.debug > 1)
        s=sprintf("single_block (%d): %e %e (residual prim, residual obj val)",...
            curr_iter, curr_res,obj_val);
        disp(s)
      end
      %store for output  
      res_iter(curr_iter) = curr_res;
    end

  end %main loop
  rac_time = toc(time_start);

  %prepare the output struct  
  x = x_current(1:n_var_x);
  obj_val = x'*mQ*x+mc'*x+model.const;
  %dual residual
  %bounds
  curr_res_d = inf;  
  %if kkt, no ineq, and eq are local, thus no dual vars avaiable
  if(~use_kkt)
    if(do_xk)
      curr_res_d = - z_current;
    end
    curr_res_d = curr_res_d + 2*mQ*x + mc - Ares'*y_current;
    curr_res_d = max(abs(curr_res_d));
  end  
  
  rac_out.sol_x = x;
  rac_out.sol_obj_val = obj_val;  
  rac_out.sol_residue_p = curr_res;
  %if(~isinf(curr_res_d))
      rac_out.sol_residue_d = curr_res_d;
  %end
  rac_out.sol_y = y_current;
  rac_out.n_iter = curr_iter;
  rac_out.rac_time = rac_time;  
  rac_out.init_time = init_time;
  rac_out.solver_time = sub_solver_time;
  rac_out.model_time = sub_model_time;
  if(isfield(run_p,'debug') && run_p.debug > 0)
    rac_out.res_iter = res_iter;
  end
end %function solve

function x = solve_subproblem_gurobi(Q_current, c_current, lb,ub, run_p)

  if(issparse(Q_current))
    model.Q = Q_current;
  else
    model.Q = sparse(Q_current);
  end
  model.obj = full(c_current);
  % no local constraints
  model.A= sparse(0,length(c_current));
  model.rhs = [];
  model.sense = repmat('=', length(model.rhs), 1);
  %add local bounds
  model.lb = lb;
  model.ub = ub;
 

  if(~isfield(run_p.gurobi_params,'TimeLimit'))
    error('Error. Time limit for gurobi sub-solver must be set')
  end
  %gurobi_write(model, 'qp.lp');
  g_out = gurobi(model,run_p.gurobi_params);
  if(~strcmpi(g_out.status, 'OPTIMAL') && ~isfield(run_p.gurobi_params,'TimeLimit')) 
    error('Error. Gurobi did not return optimal result')
  end

  if(isfield(g_out,'x'))
    x = g_out.x;
  elseif(strcmpi(g_out.status, 'TIME_LIMIT') && size(model.A,1) == 0)
    error('Time limit too short')
  else
    g_out
    error('No solution returned by gurobi')
  end
end
