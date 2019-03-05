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

function rac_out = rac_multi_block(model, run_p, tau)
% RAC_MULTI_BLOCK Solves QP using RAC multi-block approach

  %start the clock
  stime_start = tic;   
  time_start = tic;

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

  %prepare for mip, if used
  mip = run_p.mip;
  if(mip)
    x_best = model.x0;
    tol_best = inf;
    obj_best = inf;
  end

  %Define auxiliary vectors and matrice  
  Q = [model.Q, sparse(n_var_x,m2); sparse(m2,n_var_x+m2)];
  c = [model.c; sparse(m2,1)];
  c_x = model.c(1:n_var_x);

  A = [model.Aeq, sparse(m1,m2); model.Aineq, speye(m2,m2)];
  b = [model.beq;model.bineq];
  s0 = max(model.bineq-model.Aineq*model.x0,0);
  x_current = [model.x0;s0];
  y_current = zeros(m,1);
  Ares = [model.Aeq;model.Aineq];

  % use split technique for bounds if defined and not handled by partial Lagrange
  n_inf = sum(isinf(model.lb)) + sum(isinf(model.ub));
  if(mip || n_inf == 2*n_var_x)
    do_xk = 0;
  else
    do_xk = 1;
    z_current = sparse(n_var_x,1);  
    xk = model.x0;
  end

  % block sizes: minimize difference in sizes
  block_size = floor(n_var_x/run_p.n_blocks);
  o_blocks = n_var_x-block_size*run_p.n_blocks;
  RP_all = cell(0);
  RP = cell(0);
  % optimize if using RP mode, thus no overlap groups, num blocks =num groups
  if((isfield(model,'group_mode') && strcmpi('RP',model.group_mode)))
    b_size = [];
    if(strcmpi('cholesky',run_p.sub_solver_type))      
      RP_all = prepare_RP_matrices(Q,A,model.groups,beta,true, run_p.use_sparse, do_xk);
    else
      RP_all = prepare_RP_matrices(Q,A,model.groups,beta,false, run_p.use_sparse, do_xk);
    end
    use_rp = true;
    n_blocks = size(RP_all,1);
  else
    n_blocks = run_p.n_blocks;
    b_size = ones(n_blocks,1) * block_size;
    b_size(1:o_blocks) = b_size(1:o_blocks)+1; 
    use_rp = false;
  end
  %initialize counters and vars
  curr_iter = 0;
  sub_model_time = 0;
  sub_solver_time = 0;
  res_iter=[];
  curr_res = Inf;
  
  %prepare data arrays
  if(size(A,1) == 0)
    c_stat = c;
  else
    c_stat =  c-beta*A'*b;
  end
  
  Ax = A*x_current;
  Qx = Q*x_current;
 
  init_time = toc(stime_start);
  %main loop
  while(~terminate(run_p, curr_iter, curr_res, time_start))
    c_res = c_stat - A'*y_current;
    if(do_xk == 1)
      c_res = c_res  - [z_current;sparse(m2,1)];
    end

    if(use_rp)
      grp_order = randperm(n_blocks);
    else
      mblocks = get_blocks(model,b_size);
      x_blocks = mblocks.blocks;
      grp_order = [];
    end
    %update primal vars
    for block_ix = 1:n_blocks
      if(use_rp)
        RP = RP_all(grp_order(block_ix),:);
        x_ix = RP{4};
      else
        x_ix = x_blocks{block_ix};
      end
      %prepare sub-model
      stime_start = tic;
      if(user_sp)
        Q_current = model.Q(x_ix,x_ix);
        c_current = model.c(x_ix) + 2*(Qx(x_ix)-Q_current*x_current(x_ix));
      else
        A_sub = A(:,x_ix); 
        if(use_rp)          
          Q_current = RP{1};
        else
          Q_current = Q(x_ix,x_ix)+beta_h*A_sub'*A_sub;
        end
        c_current_sub = 2*(Qx(x_ix)+(beta_h)*A_sub'*Ax-Q_current*x_current(x_ix));
        c_current = c_current_sub + c_res(x_ix);
        if(do_xk)
         c_current = c_current -beta*xk(x_ix); 
         if(~use_rp)
           Q_current = Q_current+beta_h*speye(length(x_ix));
         end
        end
      end
      
      sub_model_time = toc(stime_start) + sub_model_time;    

      %call solver
      stime_start = tic;
      x_sub = solve_subproblem(Q_current, c_current, x_ix, model, run_p,x_current,...
               y_current,beta_h,RP);
      sub_solver_time = toc(stime_start) + sub_solver_time;

      %update Qx, Ax  
      diff_x=x_current(x_ix) - x_sub;      
      Qx = Qx - Q(:,x_ix)*diff_x;      
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

    %update dual y
    res_y = Ax - b;   
    if(length(res_y) > 0)
      y_current = y_current - beta*res_y;
      curr_res = max(abs(res_y));
    else
      curr_res = 0;
    end
    %update dual z if needed
    if(do_xk)
      res_z = (x_current(1:n_var_x)-xk);
      z_current = z_current - beta*res_z;
      curr_res = max(curr_res, max(abs(res_z)));
    end

    %finalize the iteration    
    if(mip)
      x = x_current(1:n_var_x);
      obj_val = x'*model.Q*x+model.c'*x+model.const;
      if(curr_res < tol_best)
        tol_best = curr_res;
        obj_best = obj_val;
        x_best = x_current(1:n_var_x);
      elseif(curr_res == tol_best && obj_val < obj_best)
        obj_best = obj_val;
        x_best = x_current(1:n_var_x);        
      end
    end
    

    curr_iter = curr_iter + 1;
    if(isfield(run_p,'debug') && run_p.debug > 0)
      x = x_current(1:n_var_x);
      if(~mip)
        obj_val = x'*model.Q*x+model.c'*x+model.const;
      end
      s=sprintf("multi_block (%d): %e %e(residual prim, residual dual, obj val)",...
            curr_iter, curr_res, obj_val);
      disp(s)      
      %store for output  
      res_iter(curr_iter) = curr_res;
    end


  end %main loop
  rac_time = toc(time_start);
  %Can not do dual residual if Gurobi was used -- missing dual vars
  if(strcmpi('gurobi',run_p.sub_solver_type))
    curr_res_d = -inf;
  else
    curr_res_d = 2*Qx(1:n_var_x) + c_x - Ares'*y_current;
    if(do_xk)
      curr_res_d = curr_res_d - z_current;
    end
    curr_res_d = max(abs(curr_res_d));
  end

  %prepare the output struct  
  if(mip)
    rac_out.sol_x = x_best;
    rac_out.sol_obj_val = obj_best;  
    rac_out.sol_residue_p = tol_best;
    rac_out.sol_dual = [];
  else
    x = x_current(1:n_var_x);
    obj_val = x'*model.Q*x+model.c'*x+model.const;
    rac_out.sol_x = x;
    rac_out.sol_obj_val = obj_val;  
    rac_out.sol_residue_p = curr_res;
    if(~isinf(curr_res_d))
      rac_out.sol_residue_d = curr_res_d;
    end
    rac_out.sol_y = y_current;
  end
  rac_out.n_iter = curr_iter;
  rac_out.rac_time = rac_time;  
  rac_out.init_time = init_time;
  rac_out.solver_time = sub_solver_time;
  rac_out.model_time = sub_model_time;
  if(isfield(run_p,'debug') && run_p.debug > 0)
    rac_out.res_iter = res_iter;
  end
end %function solve

function mblocks = get_blocks(model,b_size)

  if(isfield(model,'groups') && size(model.groups,1) > 0)
    mblocks = make_blocks(model.size, model.groups, b_size);    
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
  for block_ix = 1:n_blocks
    block_size = b_size(block_ix);
    bnew= x_ix_perm(((block_ix-1)*block_size+1):block_ix*block_size);
      blocks{block_ix} = bnew;
  end
end

function x = solve_subproblem(Q_current, c_current, x_ix, model, run_p, x_curr, ...
                y_curr,beta_h, RP)

  if(strcmpi('cholesky',run_p.sub_solver_type))
    x = solve_cholesky(Q_current, c_current,RP);
  elseif(strcmpi('gurobi',run_p.sub_solver_type))
    x = solve_subproblem_gurobi(Q_current, c_current, x_ix, model, ...
        run_p, x_curr);
  elseif(strcmpi('user_defined',run_p.sub_solver_type))
    x = run_p.sub_solver_f(Q_current, c_current, x_ix, model, run_p, ...
        x_curr, y_curr,beta_h,RP);
  else
    error('Error. Sub-solver type not recognized. Allowed: cholesky and gurobi');
  end
end

function x = solve_cholesky(Q_current, c_current, RP)
  if(size(RP,1)>0)
     R = RP{2};
     S = RP{3};
     if(size(S,1)>0)
       x = S*(R\(R'\(S'*(-c_current))));
     else
       x = R\(R'\(-c_current));
     end
  else
    x = (2*Q_current)\(-c_current);
  end
end


function x = solve_subproblem_user(Q_current, c_current, x_ix, rac_model, ...
          run_p, x_curr, y_curr,R)
  x=zeros(x_ix,1);
end


function x = solve_subproblem_gurobi(Q_current, c_current, x_ix, rac_model, run_p, x_curr)
  %disp('gurobi');
  n_vars = length(x_ix);
  model.Q = Q_current;
  model.obj = full(c_current);
  %add local bounds, if defined    

  lc = rac_model.local_constraints;
  
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

  %add local constraints, if defined
  A = sparse(0,n_vars);
  rhs = [];
  sense = [];
  N=length(c_current);
  x_not=x_curr;
  x_not(x_ix)=0;

   if(length(lc.beq) > 0)
     A = lc.Aeq(:,x_ix);
     Ax=lc.Aeq*x_not;
     rhs = lc.beq - Ax;
     sense = repmat('=', length(lc.beq), 1); 
   end
 
  if(length(lc.bineq) > 0)
    A = [A;lc.Aineq(:,x_ix)];
    Ax=lc.Aineq*x_not;
    rhs = [rhs;lc.bineq-Ax];
    sense = [sense; repmat('=', length(lc.beq), 1)]; 
  end
  %add to model
  [row,col] = find(A);
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













