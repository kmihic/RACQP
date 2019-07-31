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

  %prepare for mip, if used
  mip = run_p.mip;
  if(mip)
    x_best = model.x0;
    tol_best = inf;
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


  % use split technique for bounds if defined and not handled by partial Lagrange
  n_inf = sum(isinf(model.lb)) + sum(isinf(model.ub));
  if(mip || n_inf == 2*n_var_x)
    do_xk = 0;
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
  f_block={};
  l_block={};
  if((isfield(model,'group_mode') && (strcmpi('RP',model.group_mode) || strcmpi('CADMM',model.group_mode))))
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
      cadmm_order = randperm(n_blocks);
    end
  else
    % block sizes: minimize difference in sizes
    % if a group is designated as the first or the last, then it is put in 
    % a stand alone block
    % first -> update done after the dual, but before any other primal block
    % last -> update done after all other promal blocks, but before the dual update 

    if(isfield(model,'group_mode') && isfield(model,'groups') &&...
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
  res_iter=[];
  obj_iter = [];
  time_iter = [];
  curr_res = Inf;
  cnt = 0;
  
  %prepare data arrays
  if(size(A,1) == 0)
    c_stat = c;
  else
    c_stat =  c-beta*A'*b;
  end
  
  Ax = A*x_current;
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
  while(~terminate(run_p, curr_iter, curr_res, time_start))
    c_res = c_stat - A'*y_current;
    if(do_xk == 1)
      c_res = c_res  - [z_current;zero_m2];
    end
    
    if(use_cyclic_admm)      
      grp_order = cadmm_order;
    elseif(use_rp)
      grp_order = randperm(n_blocks);
    else
      mblocks = get_blocks(model,b_size, no_split, n_vars_grp); 
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
      if(Q_handle)
        obj_val = mQ(':',':',x,x)+mc'*x+model.const;
      else
        obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
      end
      if( curr_res < tol_best)
        tol_best = curr_res;
        obj_best = obj_val;
        x_best = x_current(1:n_var_x);
        cnt=0;
      elseif(curr_res == tol_best && obj_val < obj_best)
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
      if(~mip)
        if(Q_handle)
          obj_val = mQ(':',':',x,x)+mc'*x+model.const;
        else
          obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
        end
      end
      if(mip)
         s=sprintf("multi_block (%d): %e %e %e (residual prim, obj val, best obj val)",...
            curr_iter, curr_res, obj_val, obj_best);

      else
        s=sprintf("multi_block (%d): %.15e %e(residual prim,  obj val) %f",...
            curr_iter, curr_res, obj_val, toc(time_start));
      end
      disp(s)   
      %store for output  
      res_iter(curr_iter) = curr_res;
      obj_iter(curr_iter) = obj_val;
      time_iter(curr_iter) = toc(time_start);
    end


  end %main loop
  rac_time = toc(time_start);
  %Can not do dual residual if Gurobi was used -- missing dual vars
  if(strcmpi('gurobi',run_p.sub_solver_type) || run_p.no_calc_dual_res)
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
    rac_out.sol_residue_d = -inf;
    rac_out.sol_y = [];
  else
    x = x_current(1:n_var_x);
     if(Q_handle)
      obj_val = mQ(':',':',x,x)+mc'*x+model.const;
     else
      obj_val = x'*mQ(:,:)*x+mc'*x+model.const;
    end
    rac_out.sol_x = x;
    rac_out.sol_obj_val = obj_val;  
    rac_out.sol_residue_p = curr_res;
    %if(~isinf(curr_res_d))
      rac_out.sol_residue_d = curr_res_d;
    %end
    rac_out.sol_y = y_current;
  end
  rac_out.n_iter = curr_iter;
  rac_out.rac_time = rac_time;  
  rac_out.init_time = init_time;
  rac_out.solver_time = sub_solver_time;
  rac_out.model_time = sub_model_time;
  if(isfield(run_p,'debug') && run_p.debug > 0)
    rac_out.res_iter = res_iter;
    rac_out.obj_iter = obj_iter;
    rac_out.time_iter = time_iter;
  end
end %function solve

function mblocks = get_blocks(model,b_size,no_split, n_vars)

  if(isfield(model,'groups') && size(model.groups,1) > 0) 
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


function x = solve_subproblem_user(Q_current, c_current, x_ix, rac_model, ...
          run_p, x_curr, y_curr,R)
  x=zeros(x_ix,1);
end


function x = solve_subproblem_gurobi(Q_current, c_current, x_ix, rac_model, run_p, x_curr)
  %disp('gurobi');
  n_vars = length(x_ix);
  if(~issparse(Q_current))
    model.Q = sparse(Q_current);
  else
    model.Q = Q_current;
  end
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

  if(length(lc.beq) > 0)
    A = lc.Aeq(:,x_ix);
    rhs = lc.beq;
    sense = repmat('=', length(lc.beq), 1);
  end
  if(length(lc.bineq) > 0)
    A = [A;lc.Aineq(:,x_ix)];
    rhs = [rhs;lc.bineq];
    sense = [sense; repmat('<', length(lc.bineq), 1)];
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
  elseif(strcmpi(g_out.status, 'TIME_LIMIT') && size(model.A,1) == 0)
    error('Time limit too short')
  else
    g_out
    error('No solution returned by gurobi')
  end
end


















