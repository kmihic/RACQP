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

function [sol res] = use_mosek(m, max_time, epsilon)
  
  % mosek does not return dual variable for variable 
  % lower/upper bounds. Converting bounds to constraints 
  % to get duals for dual residual to match Gurobi.
  % Note: Mosek returns the duals as slx and sux output values
  n = m.size;
  B=[-speye(n);speye(n)];
  b=[-m.lb;m.ub];
  ix = ~isfinite(b);
  B(ix,:)=[];
  b(ix)=[];
  m.Aineq = [m.Aineq;B];
  m.bineq = [m.bineq;b];
  m.lb = ones(n,1)*(-Inf);  
  m.ub = ones(n,1)*(Inf);

  Q = 2*m.Q;
  c = m.c;
  A = [m.Aeq;m.Aineq];
  clb = [m.beq;-inf(length(m.bineq),1)];
  cub = [m.beq;m.bineq];
  lb = m.lb;
  ub = m.ub;

  param.MSK_IPAR_INTPNT_MULTI_THREAD= 'MSK_OFF' ;
  param.MSK_IPAR_NUM_THREADS=1;
  if(max_time > 0)
    param.MSK_DPAR_OPTIMIZER_MAX_TIME = max_time;
  end
  if(epsilon > 0) 
    %Dual feasibility tolerance 
    param.MSK_DPAR_INTPNT_QO_TOL_DFEAS=epsilon;
    %Primal feasibility tolerance
    param.MSK_DPAR_INTPNT_QO_TOL_PFEAS=epsilon;
    %Relative gap termination tolerance used by the interior-point optimizer for quadratic problems.
    param.MSK_DPAR_INTPNT_QO_TOL_REL_GAP=epsilon;
  end

  s_time = tic;
  res = mskqpopt(Q,c,A,clb,cub,lb,ub,param,'minimize info');

  % prepare the output to be read by run*test scripts
  g_out = res.sol.itr;
  if(strcmpi(g_out.solsta, 'OPTIMAL'))
    sol.sol_x = g_out.xx;
    sol.sol_y = g_out.y;
    sol.sol_obj_val = g_out.pobjval;    
    sol.n_iter = res.info.MSK_IINF_INTPNT_ITER;

    % y (double[]) – Identical to sol.slc-sol.suc (not in integer solution).
    % slc (double[]) – Dual variable for constraint lower bounds (not in integer solution).
    % suc (double[]) – Dual variable for constraint upper bounds (not in integer solution).
    % slx (double[]) – Dual variable for variable lower bounds (not in integer solution).
    % sux (double[]) – Dual variable for variable upper bounds (not in integer solution).
    % get residuals
    % turning off dual for bounds; having them as inequalities
    do_bounds = false;
    [res_p res_d] = get_residuals(m, g_out.xx, g_out.y, true, ...
                          false, (g_out.slx-g_out.sux), do_bounds);    
    sol.sol_res_p = res_p;
    sol.sol_res_d = res_d;
  else
    disp(res.rcodestr)
    sol.sol_x = [];
    sol.sol_y = [];
    sol.sol_obj_val = inf;
    sol.n_iter = -inf;
    sol.sol_res_p.abs = -inf;
    sol.sol_res_d.abs = -inf;
    sol.sol_res_p.rel = -inf;
    sol.sol_res_d.rel = -inf;
    sol.sol_res_p.L2 = -inf;
    sol.sol_res_d.L2 = -inf;
  end
  sol.runtime = toc(s_time);

end









