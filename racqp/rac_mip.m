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

function rac_out = rac_mip(model, run_mip, time_start)


  run_sub = run_mip.run_sub; % subproblem run params
  run_sub.max_rtime = run_mip.max_rtime;
  run_sub.n_perturb_trial = run_mip.n_perturb_trial;

  %if not external,make sure we use gurobi
  if(strcmpi('cholesky',run_sub.sub_solver_type))
    run_sub.sub_solver_type = 'gurobi'; %make sure we use Gurobi
  end
  %set var type to be used by Gurobi
  vtype = repmat('C', model.size, 1);
  vtype(model.integers) = 'I';
  vtype(model.binary) = 'B';
  model.vtype = vtype;

  %initialize counters and vars
  curr_iter = 0;
  time_start = tic;
  obj_val_iter=[];
  res_best = Inf;
  obj_best = Inf;
  x_best = model.x0; 
  x0_curr = model.x0;
  cnt = 0;
  cnt_perturb = 0;
  obj_iter = [];
  res_iter = [];
  time_iter = [];

  %main loop
  while(~terminate(run_mip, curr_iter, cnt_perturb, time_start))

    %set the initial point
    model.x0 = x0_curr;
    %call the solver
    rac_current = rac_multi_block(model,run_sub,time_start);

    obj_curr = rac_current.sol_obj_val;
    res_curr = rac_current.sol_residue_p;

    %better solution found?
    %closer to feasibility - take the result
    if(res_curr < res_best) 
      res_best = res_curr;
      obj_best = obj_curr;
      x_best = rac_current.sol_x;
    %same tolerance, get if objective better
    elseif(res_curr == res_best && obj_curr < obj_best)
      obj_best = obj_curr;
      x_best = rac_current.sol_x;
    end

    %time to perturb
    x0_curr = perturb(x_best, model, run_mip, res_curr);
    cnt_perturb = cnt_perturb + 1;
    

    curr_iter = curr_iter +rac_current.n_iter;
    if(isfield(run_mip,'debug') && run_mip.debug > 0)      
      s=sprintf("MIP (%d): %.3e %.3e : %.3e %.3e %.3f (curr residual, obj val: best res,obj, time)",...
               curr_iter,res_curr,obj_curr,res_best,obj_best,toc(time_start));
      disp(s)  
      %store for output 
      if(isfield(rac_current,'obj_iter'))
        obj_iter = [obj_iter,rac_current.obj_iter] ;
      end
      if(isfield(rac_current,'res_iter'))
        res_iter = [res_iter,rac_current.res_iter] ;
      end
      if(isfield(rac_current,'time_iter'))
        time_iter = [time_iter,rac_current.time_iter] ;
      end
    end

  end %main loop
  rac_time = toc(time_start);

  rac_out.sol_x = x_best;
  rac_out.sol_obj_val = obj_best;  
  rac_out.sol_residue = res_best;
  rac_out.n_iter = curr_iter;
  rac_out.n_perturb = cnt_perturb;
  rac_out.rac_time = rac_time;  
  if(isfield(run_mip,'debug') && run_mip.debug > 0)
    rac_out.obj_iter = obj_iter;
    rac_out.res_iter = res_iter;
    rac_out.time_iter = time_iter;
  end

end

function term = terminate(r, curr_iter, cnt_perturb, time_start)

  term = false;
  if(curr_iter >= r.max_iter ...
     || cnt_perturb >= r.max_nperturb ...
     || toc(time_start) >= r.max_rtime)
     term = true;
  end

end


function x0 = perturb(x, model, run_mip, res_p)

  %find vars to permute
  n_var = length(x);
  n_ch = get_num_vars(run_mip);
  
  %decide on perm type to do
  if(strcmpi('user_defined',run_mip.permute_type))
    x0 = run_mip.permute_f(x, model,n_ch, res_p);
  elseif(strcmpi('swap', run_mip.permute_type))
    x0 = perturb_swap(x,  model, n_ch);
  elseif(strcmpi('permute', run_mip.permute_type))
    v_perm = datasample([1:n_var],n_ch,'Replace',false);
    x0 = x;
    for ii = v_perm
      is_int = length(find(model.integers == ii) + find(model.binary == ii)) > 0;
      x0(ii) = get_rnd_val_within_bounds(model.lb(ii), model.ub(ii),x(ii), is_int );
    end
  else
    error('Perturbation type not recognized')
  end
  
end

function n = get_num_vars(r)

  if(strcmpi('Uniform', r.permute_dist))
    n = random('Uniform', r.permute_min, r.permute_max);
  elseif(strcmpi('Normal', r.permute_dist))    
    n = random('Normal', r.permute_mu, r.permute_std);
    n = max(r.permute_min, min(r.permute_max, n));
  elseif(strcmpi('Exponential', r.permute_dist))    
    n = random('Exponential', r.permute_mu);
    n = max(r.permute_min, min(r.permute_max, n));
  else
    error('Error. Distribution not supported (only U, N and Exp)')
  end
  n = round(n);
  
end

function x = get_rnd_val_within_bounds(lb, ub, x_curr, is_int)  

   if(lb == 0 && ub == 1) %binary var, do swap
     x = 1 - x_curr;
     return;
  elseif( isinf(lb) && isinf(ub) ) %no limits
    x = randn(1);
  elseif( isinf(ub) )  %have a lower bound
    x = lb+abs(randn(1));
  elseif( isinf(ub) ) %have an upper bound
    x = ub - abs(randn(1));
  else %have both bounds
    x = rand(1)*(ub-lb)+lb;
  end
  if(is_int)
    x=round(x);
    %do not want the same value as it was
    if(x == x_curr)
      x = x+1;    
      if(x>ub)
        x = ub-1;
      end
    end
  end
end














