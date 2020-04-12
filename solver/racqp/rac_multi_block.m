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


function rac_out = rac_multi_block(model, run_p, time_start)
% RAC_MULTI_BLOCK Solves QP using RAC multi-block approach
  
  %define global constants;
  global CHOLESKY_SOL GUROBI_SOL USER_SOL;
  CHOLESKY_SOL = 1;
  GUROBI_SOL = 2;
  USER_SOL = 3;


  %start the clock
  stime_start = time_start;   
  
  %check if using user defined sub-problem solver
  if(strcmpi('user_defined',run_p.sub_solver_type))
    user_sp = true;
  else
    user_sp = false;
  end

  n_var_x = model.size;
  m1 = size(model.Aeq,1);
  m2 = size(model.Aineq,1);  
  m = m1+m2;  
  beta_i = 1/run_p.beta;
  beta_h = run_p.beta/2;
  beta = run_p.beta;


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
  do_dual_bounds_gurobi = false;

  %prepare for mip, if used
  mip = run_p.mip;
  if(mip)
    x_best = model.x0;
    tol_best_mip = inf;
    obj_best = inf;
  end

  %Define auxiliary vectors and matrice  
   mQ = model.Q;
   mc = model.c;
   if(isa(model.Q,'function_handle'))
     Q_handle = true;
     if(m2 > 0) 
      error('Not coded for Q as a function handle and inequality');
     else
      Q = model.Q;
     end
   else
     Q_handle = false;
     Q = [model.Q, sparse(n_var_x,m2); sparse(m2,n_var_x+m2)];
   end

   c = [model.c; sparse(m2,1)];
   c_x = model.c(1:n_var_x);
   A = [model.Aeq, sparse(m1,m2); model.Aineq, speye(m2,m2)];   
   A_T = [model.Aeq; model.Aineq]';
   b = [model.beq;model.bineq];
   s0 = max(model.bineq-model.Aineq*model.x0,0);
   x_current = [model.x0;s0];
   y_current = zeros(m,1);
   Ares = [model.Aeq;model.Aineq];
   zero_m2=sparse(m2,1);
   % force dense mode
   if(~run_p.use_sparse)
    % mQ and mC are used for opj val calc and user defined sub-problem solver
    % leave it as defined by the model
    %mQ = sparse(mQ);
    %mc = sparse(mc);
    if(~Q_handle)
      Q = full(Q);
    end
    A = full(A);
    c = full(c);
    c_x = full(c_x);
    s0 = full(s0);
    Ares = full(Ares);
    zero_m2 = zeros(m2,1);
   end

  % sub problem solvers
  use_gurobi = false;
  if(strcmpi('cholesky',run_p.sub_solver_type))
    solver = CHOLESKY_SOL;
  elseif(strcmpi('gurobi',run_p.sub_solver_type))
    solver = GUROBI_SOL;
    use_gurobi = true;
  elseif(strcmpi('user_defined',run_p.sub_solver_type))
    solver = USER_SOL;
  end

  bounds_present_lc = false;
  if(isfield(model,'local_constraints'))
    lc = model.local_constraints;  
    n_inf = sum(isinf(lc.lb)) + sum(isinf(lc.ub));
    if(n_inf ~= 2*n_var_x)
      bounds_present_lc = true;
    end
  end

  n_inf = sum(isinf(model.lb)) + sum(isinf(model.ub));
  if(n_inf == 2*n_var_x)
    bounds_present = false;
  else
    bounds_present = true;
  end
  % use split technique for bounds if defined and not handled by partial Lagrange
  if(mip || ~bounds_present)
    do_xk = 0;
    z_current = [];
  elseif(use_gurobi && bounds_present_lc) 
    do_xk = 0;
    if(use_dual_res)
      z_current = zeros(model.size,1);      
      do_dual_bounds_gurobi = true;
    else
      z_current = [];      
      do_dual_bounds_gurobi = false;
    end
  else
    do_xk = 1;
    z_current = zeros(n_var_x,1);  
    xk = model.x0;
  end
  
  RP_all = cell(0);
  RP = cell(0);
  if(isfield(model,'group_mode')&& strcmpi('RAC_NO_SPLIT',model.group_mode))
    no_split = true;
  else
    no_split = false;
  end
  % optimize if using RP mode, thus no overlap groups, num blocks =num groups  
  use_cyclic_admm = false;  
  use_distrib_admm = false;
  f_block={};
  l_block={};
  if((isfield(model,'group_mode') && (strcmpi('RP',model.group_mode) ...
      || strcmpi('CADMM',model.group_mode)...
      || strcmpi('DADMM',model.group_mode))))
    b_size = [];
    if(strcmpi('cholesky',run_p.sub_solver_type))      
      RP_all = prepare_RP_matrices(Q,A,model.groups,beta,true, run_p.use_sparse, do_xk);
    else
      RP_all = prepare_RP_matrices(Q,A,model.groups,beta,false, run_p.use_sparse, do_xk);
    end
    use_rp = true;
    n_blocks = size(RP_all,1);
    if(strcmpi('CADMM',model.group_mode))
      use_cyclic_admm = true;
      %no needed to do rnd; keep user's preference
%      cadmm_order = randperm(n_blocks);
      cadmm_order = 1:n_blocks;
    elseif(strcmpi('DADMM',model.group_mode))
      use_distrib_admm = true;
      dadmm_order = 1:n_blocks;
    end
  else
    % block sizes: minimize difference in sizes
    % if a group is designated as the first or the last, then it is put in 
    % a stand alone block
    % first -> update done after the dual, but before any other primal block
    % last -> update done after all other promal blocks, but before the dual update 

    if(isfield(model,'group_mode') && isfield(model,'groups') &&...
        size(model.groups,2) >= 3 && ...
           (strcmpi(model.group_mode, 'RAC_NO_SPLIT') || ...
            strcmpi(model.group_mode, 'RAC')))
      n_grp = size(model.groups,1);
      n_first = 0;
      n_last = 0;
      cnt_f = 0;
      cnt_l = 0;
      rm_grp = [];
      for grp = 1:n_grp
        pos = model.groups{grp,3};
        if(pos == 0)
          g_vars = model.groups{grp,1};
          n_first = n_first + length(g_vars);
          cnt_f = cnt_f + 1;
          f_block{cnt_f} = g_vars;
          rm_grp = [rm_grp; grp];
        elseif(pos == 2)
          g_vars = model.groups{grp,1};
          n_last = n_last + length(g_vars);
          cnt_l = cnt_l + 1;
          l_block{cnt_l} = g_vars;
          rm_grp = [rm_grp; grp];
        end
      end
      % remove first/last vars from sizing
      n_vars_grp = n_var_x - n_first - n_last;
      %remove the groups
      model.groups(rm_grp,:)=[]; 
    else
     n_vars_grp = n_var_x;
    end

    block_size = floor(n_vars_grp/run_p.n_blocks);
    o_blocks = n_vars_grp-block_size*run_p.n_blocks;
    n_blocks = run_p.n_blocks;
    b_size = ones(n_blocks,1) * block_size;
    b_size(1:o_blocks) = b_size(1:o_blocks)+1; 
    use_rp = false;
  end
    
  %initialize counters and vars
  curr_iter = 0;
  sub_model_time = 0;
  sub_solver_time = 0;
  x_iter=[];
  y_iter=[];
  z_iter=[];
  res_iter=[];
  res_iter_dual=[];
  obj_iter = [];
  time_iter = [];
  curr_res.prim = Inf;
  if(use_dual_res)
    curr_res.dual = Inf;
  else 
    curr_res.dual = -Inf;
  end
  cnt = 0;
  
  %prepare data arrays
  Ax = A*x_current;
  if(use_distrib_admm)
    l_sum = (sum(Ax)-model.beq)/n_blocks;
    l_current = get_lambda_current(RP_all, n_blocks,l_sum,A,x_current);
  end
  if(size(A,1) == 0)
    c_stat = c;
  elseif(use_distrib_admm)
    c_stat =  c;
  else
    c_stat =  c-beta*A'*b;
  end
  
  
  if(length(find(x_current)) == 0)
    Qx = zeros(length(x_current),1);
  else
    if(Q_handle)
      Qx = Q(':',':',x_current);
      if(isfield(run_p,'debug') && run_p.debug > 0)
       disp("Initial point not zero, using Q as a function. This may be an (memory) expensive operation")
      end
    else
      Qx = Q*x_current;
    end
  end
  init_time = toc(stime_start);
  %main loop
  while(~terminate(run_p, curr_iter, curr_res, time_start, mip))
    c_res = c_stat - A'*y_current;
    if(do_xk == 1)
      c_res = c_res  - [z_current;zero_m2];
    end
    
    if(use_cyclic_admm)      
      grp_order = cadmm_order;
    elseif(use_distrib_admm)
      grp_order = dadmm_order;
    elseif(use_rp)
      grp_order = randperm(n_blocks);
    else
      mblocks = get_blocks(model,b_size, no_split, n_vars_grp,run_p); 
      x_blocks = [f_block, mblocks.blocks, l_block];
      grp_order = [];
      %if using no split, num_blocks can change
      n_blocks = length(x_blocks);
    end
    %update primal vars
    for block_ix = 1:n_blocks
      if(use_rp)
        RP = RP_all(grp_order(block_ix),:);
        x_ix = RP{3};
      else
        x_ix = x_blocks{block_ix};
      end
      %prepare sub-model
      stime_start = tic;
      if(user_sp)  
        Q_current = mQ(x_ix,x_ix);
        c_current = mc(x_ix) + 2*(Qx(x_ix)-Q_current*x_current(x_ix));
      else
        A_sub = A(:,x_ix); 
        if(use_rp)   
          if(Q_handle) 
            % get the row here we'll need it a couple of rows later anyway
            Qh = Q(':',x_ix);
          end
          Q_current = RP{4};
        else
          if(Q_handle)
            Qh = Q(':',x_ix);
            Q_current = Qh(x_ix,:)+beta_h*A_sub'*A_sub;
          else
            Q_current = Q(x_ix,x_ix)+beta_h*A_sub'*A_sub;
          end
        end
        if(use_distrib_admm)
          c_current = c_res(x_ix) - beta * l_current(block_ix);
        else
          c_current_sub = 2*(Qx(x_ix)+(beta_h)*A_sub'*Ax-Q_current*x_current(x_ix));
          c_current = c_current_sub + c_res(x_ix);
        end
        if(do_xk)
         c_current = c_current -beta*xk(x_ix); 
         if(~use_rp)
           Q_current = Q_current+beta_h*speye(length(x_ix));
         end
        end
      end
      
      sub_model_time = toc(stime_start) + sub_model_time;    
      stime_start = tic;

      %call solver
      [x_sub z_sub] = solve_subproblem(Q_current, c_current, x_ix, model, run_p,x_current,...
               y_current,beta_h,RP, solver);
      if(do_dual_bounds_gurobi)        
        z_current(x_ix) = z_sub;
      end

      sub_solver_time = toc(stime_start) + sub_solver_time;
      stime_start = tic;

      %update Qx, Ax  
      diff_x=x_current(x_ix) - x_sub; 
      if(Q_handle)
          Qx = Qx - Qh*diff_x;
      else
        Qx = Qx - Q(:,x_ix)*diff_x;  
      end
      Ax = Ax - A(:,x_ix)*diff_x;
      sub_model_time = toc(stime_start) + sub_model_time;   
      %update x
      x_current(x_ix) = x_sub;
      %update auxiliary x, if needed
      if(do_xk == 1)
        xk(x_ix) = min(max(model.lb(x_ix),x_current(x_ix)-beta_i*z_current(x_ix)),model.ub(x_ix));      
      end
    end %x block

    stime_start = tic;
    %update s if needed
    if(length(model.bineq) > 0)
      A_ineqx1=(Ax(m1+1:m)-x_current(n_var_x+1:n_var_x+m2));
      s_current = max(0,beta_i*y_current(m1+1:m)+model.bineq-A_ineqx1);
      %update Ax
      Ax(m1+1:m) = Ax(m1+1:m) - x_current(n_var_x+1:n_var_x+m2) + s_current;
      %update s (last block)
      x_current(n_var_x+1:n_var_x+m2)=s_current; 
    end

    % update lambda and c_stat if using distributed admm
    if(use_distrib_admm)
      l_sum = (sum(Ax)-model.beq)/n_blocks;
      l_current = get_lambda_current(RP_all, n_blocks,l_sum,A,x_current); 
    end
    %update dual y
    res_y = Ax - b;   
    if(length(res_y) > 0)
      if(use_distrib_admm)
        y_current = y_current - beta*l_sum;
      else
        y_current = y_current - run_p.stepY*beta*res_y;
      end
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
    %update dual z if needed and get z residual
    if(do_xk)
      diff_z = (x_current(1:n_var_x)-xk);
      z_current = z_current - run_p.stepZ*beta*diff_z;
      res_z = norm(diff_z,Inf);
      if(do_rel_tol)
        p_res_z = res_z / (1+max(norm(x_current(1:n_var_x), Inf),...
                             norm(xk, Inf)));  
      else
        p_res_z = res_z;
      end
    end    
    if(use_dual_res && (do_xk || do_dual_bounds_gurobi))
      d_res_z = norm(z_current,Inf);
    end

    % calc the overall residuals
    % primal
    curr_res.prim = max(p_res_z,max(p_res_eq, p_res_ineq));
    % dual
    if(use_dual_res)
      Ay = A_T*y_current;
      curr_res_d = 2*Qx(1:n_var_x) + model.c - Ay ;
      if(do_xk || do_dual_bounds_gurobi)
        curr_res_d = curr_res_d - z_current;
      end
      d_res_abs = norm(curr_res_d, Inf);
      if(do_rel_tol)
        curr_res.dual = d_res_abs/(1+max(max(norm(Qx,Inf), c_max), ...
                                     max(norm(Ay,Inf), d_res_z)));
      else
         curr_res.dual = d_res_abs;
      end
%      else
%        curr_res.dual = 0;
     end



    %finalize the iteration for MIP
    if(mip)
      x = x_current(1:n_var_x);
      if(Q_handle)
        obj_val = mQ(':',':',x,x)+mc'*x+model.const;
      else
        obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
      end
      if( curr_res.prim < tol_best_mip)
        tol_best_mip = curr_res.prim;
        obj_best = obj_val;
        x_best = x_current(1:n_var_x);
        cnt=0;
      elseif(curr_res.prim == tol_best_mip && obj_val < obj_best)
        obj_best = obj_val;
        x_best = x_current(1:n_var_x);    
        cnt = 0;
      else
        cnt = cnt+1;
      end
      %time to perturb
      if(cnt == run_p.n_perturb_trial)
        break;
      end  
    end    

    curr_iter = curr_iter + 1;
    if(isfield(run_p,'debug') && run_p.debug > 0)
      x = x_current(1:n_var_x);
      % no need for obj val if not mip. Using dual residual
    %  if(~mip)
    %    if(Q_handle)
    %      obj_val = mQ(':',':',x,x)+mc'*x+model.const;
    %    else
    %      obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
    %    end
    %  end
      if(mip)
         s=sprintf("multi_block (%d): %e %e %e (residual prim, obj val, best obj val)",...
            curr_iter, curr_res.prim, obj_val, obj_best);
         obj_iter(curr_iter) = obj_val;

      else
        s=sprintf("multi_block (%d): %e %.3e (res prim, dual) %3f [s]",...
            curr_iter, curr_res.prim, curr_res.dual, toc(time_start));
      end
      disp(s)   

      %store for output  
      x_iter(curr_iter,:) = x_current;
      y_iter(curr_iter,:) = y_current;
      if(do_xk)
        z_iter(curr_iter,:) = z_current;
      end
      res_iter(curr_iter) = curr_res.prim;
      res_iter_dual(curr_iter) = curr_res.dual;
      time_iter(curr_iter) = toc(time_start);
    end
  % END WHILE

  end %main loop
  rac_time = toc(time_start);  

  %prepare the output struct  
  if(mip)
    rac_out.sol_x = x_best;
    rac_out.sol_obj_val = obj_best;  
    rac_out.sol_res_p = tol_best_mip;
    rac_out.sol_res_d = -inf;
    rac_out.sol_y = [];
    rac_out.sol_z = [];
  else    
    x = x_current(1:n_var_x);
    rac_out.sol_x = x;
    rac_out.sol_y = y_current;
    rac_out.sol_z = z_current;
    if(Q_handle)
      obj_val = mQ(':',':',x,x)+mc'*x+model.const;
    else
      obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
    end
    rac_out.sol_obj_val = obj_val;  
    

    % get residuals
    if(do_dual_bounds_gurobi)
      do_xk = true;
    end 
    
    [res_p res_d] = get_residuals(model, x, y_current, run_p.calc_dual_res, ...
                          Q_handle, z_current, do_xk);
    rac_out.sol_res_p = res_p;
    rac_out.sol_res_d = res_d;
  end

  rac_out.n_iter = curr_iter;
  rac_out.runtime = rac_time;  
  rac_out.init_time = init_time;
  rac_out.solver_time = sub_solver_time;
  rac_out.model_time = sub_model_time;
  if(isfield(run_p,'debug') && run_p.debug > 0)
    rac_out.x_iter = x_iter;
    rac_out.y_iter = y_iter;
    rac_out.z_iter = z_iter;
    rac_out.res_iter = res_iter;
    rac_out.res_iter_dual = res_iter_dual;
    if(length(obj_iter)>0)
      rac_out.obj_iter = obj_iter;
    end
    rac_out.time_iter = time_iter;
  end
  

end %function solve

function mblocks = get_blocks(model,b_size,no_split, n_vars, run_p)

  if(isfield(run_p, 'make_blocks_f'))
    mblocks = run_p.make_blocks_f(n_vars, model.groups, b_size, no_split);
  elseif(isfield(model,'groups') && size(model.groups,1) > 0) 
    if(no_split)
      mblocks = make_blocks_no_split(n_vars, model.groups, b_size); 
    else
      mblocks = make_blocks(n_vars, model.groups, b_size, no_split); 
    end
  else
    if(length(b_size) == 0)
      error('Error. Blocks not defined')
    end
    mblocks.blocks = get_x_blocks(model.size,b_size);
    mblocks.grp_order = [];
  end
end

%simple grouping, no super-vars
function blocks = get_x_blocks(n_var,b_size)
  n_blocks = length(b_size);
  x_ix_perm = randperm(n_var);
  blocks = cell(n_blocks,1);
  bb = 0;
  for block_ix = 1:n_blocks
    block_size = b_size(block_ix);
    be = bb+block_size;
    bnew= x_ix_perm((bb+1):be);
    blocks{block_ix,1} = bnew;
    bb = be;
    % BUG!!!!
    % variables skipped on overlap to regular size blocks
%     bnew= x_ix_perm(((block_ix-1)*block_size+1):block_ix*block_size);
%     blocks{block_ix,1} = bnew;

  end
end

function [x z_gurobi] = solve_subproblem(Q_current, c_current, x_ix, model, run_p, x_curr, ...
                y_curr,beta_h, RP,  solver)
  
  global CHOLESKY_SOL GUROBI_SOL USER_SOL;
  z_gurobi = [];
  if(solver == CHOLESKY_SOL)
    x = solve_cholesky(Q_current, c_current,RP);
  elseif(solver == GUROBI_SOL)
    [x z_gurobi] = solve_subproblem_gurobi(Q_current, c_current, x_ix, model, ...
        run_p, x_curr);
  elseif(solver == USER_SOL)
    x = run_p.sub_solver_f(Q_current, c_current, x_ix, model, run_p, ...
        x_curr, y_curr,beta_h,RP);
  else
    error('Error. Sub-solver type not recognized. Allowed: cholesky and gurobi');
  end
end

function x = solve_cholesky(Q_current, c_current, RP)
  if(size(RP,1)>0)
     R = RP{1};
     S = RP{2};
     if(size(S,1)>0)
       x = S*(R\(R'\(S'*(-c_current))));
     else
       x = R\(R'\(-c_current));
     end
  else
    x = (2*Q_current)\(-c_current);
  end
end



function [x z] = solve_subproblem_gurobi(Q_current, c_current, x_ix, ...
                  rac_model, run_p, x_curr)
    
  n_vars = length(x_ix);

  % take care of the objective
  if(~issparse(Q_current))
    model.Q = sparse(Q_current);
  else
    model.Q = Q_current;
  end
  model.obj = full(c_current);

  %add local bounds, if defined 
  A=[];
  dual_bounds = false;
  if(isfield(rac_model,'local_constraints'))
    lc = rac_model.local_constraints;  

    if(run_p.use_dual_res)    
      % gurobi does not return dual variable for variable 
      % lower/upper bounds. Converting bounds to constraints 
      % to get duals for dual residual.  
      lb = lc.lb(x_ix);
      lb_ix = isfinite(lb);
      n_lb = length(find(lb_ix));    
      ub = lc.ub(x_ix);
      ub_ix = isfinite(ub);
      n_ub = length(find(ub_ix));

      if(n_lb+n_ub > 0)
        %A = [-speye(n_lb);speye(n_ub)];      
        A = sparse(n_lb+n_ub,n_vars);
        if(n_lb > 0)
          a_ix = sub2ind(size(A),1:n_lb, find(lb_ix)');
          A(a_ix) = -1;
        end
        if(n_ub > 0)
          a_ix = sub2ind(size(A),(1+n_lb):(n_lb+n_ub), find(ub_ix)');
          A(a_ix) = 1;
        end
        rhs = [-lb(lb_ix);ub(ub_ix)];
        sense = [repmat('<', length(rhs), 1)];  
      else
        A = sparse(0,n_vars);
        rhs = [];
        sense = [];
      end
      model.lb = ones(n_vars,1)*(-Inf);  
      model.ub = ones(n_vars,1)*(Inf);
      run_p.gurobi_params.presolve = 1;
      dual_bounds = true;
    else %no dual residual needed
      A = sparse(0,n_vars);
      rhs = [];
      sense = [];
      if(length(lc.lb) > 0)
        model.lb = lc.lb(x_ix);
      else
        model.lb = -Inf * ones(length(x_ix),1);
      end
      if(length(lc.ub) > 0)
        model.ub = lc.ub(x_ix);
      else
        model.ub = Inf * ones(length(x_ix),1);
      end
    end

    %add local constraints, if defined
    if(length(lc.beq) > 0)
      A = [A;lc.Aeq(:,x_ix)];
      rhs = [rhs;lc.beq];
      sense = [sense;repmat('=', length(lc.beq), 1)];
    end
    if(length(lc.bineq) > 0)
      A = [A;lc.Aineq(:,x_ix)];
      rhs = [rhs;lc.bineq];
      sense = [sense; repmat('<', length(lc.bineq), 1)];
    end
  end

  %add to model
  [row,col] = find(A);
  row = unique(row);
  if(length(row) > 0)
    model.A = A(row,:);
    model.rhs = rhs(row);
    model.sense = sense(row);
  else
    model.A= sparse(0,n_vars);
    model.rhs = [];
    model.sense = repmat('=', length(model.rhs), 1);
 end
 
  %if MIP, define var type
  if(isfield(rac_model,'vtype'))
     model.vtype = rac_model.vtype(x_ix);
  end

   
%run_p.gurobi_params.outputflag=1;

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
    if(dual_bounds)
      % first n_lb + n_ub dual vars are bounds  
      z_lb = zeros(n_vars,1);
      z_ub = zeros(n_vars,1);
      if(n_lb > 0 && n_ub > 0)
        z_lb(lb_ix) = g_out.pi(1:n_lb);
        z_ub(ub_ix) = g_out.pi(1+n_lb:n_lb+n_ub);
      elseif(n_lb > 0)        
        z_lb(lb_ix) = g_out.pi(1:n_lb);
      elseif(n_ub > 0)
        z_ub(ub_ix) = g_out.pi(1:n_ub);
      end
      z = -z_lb+z_ub;
      g_out.pi = g_out.pi(n_lb+n_ub+1:end);
    else
      z = [];
    end %dual bounds
  elseif(strcmpi(g_out.status, 'TIME_LIMIT') && size(model.A,1) == 0)
    error('Time limit too short')
  else
    g_out
    error('No solution returned by gurobi')
  end
end


function l_current = get_lambda_current(RP_all, n_blocks,l_sum,A,x_current)

  l_current = [];
  for block_ix=1:n_blocks
    RP = RP_all(block_ix,:);
    x_ix = RP{3};
    l_i = A(:,x_ix)*x_current(x_ix)-l_sum;
    l_current =  [l_current,l_i];
  end
end















