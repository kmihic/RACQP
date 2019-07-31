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

function sol = MOSEK(m)

  Q = 2*m.Q;
  c = m.c;
  A = [m.Aeq;m.Aineq];
  clb = [m.beq;-inf(length(m.bineq),1)];
  cub = [m.beq;m.bineq];
  lb = m.lb;
  ub = m.ub;

  param.MSK_IPAR_INTPNT_MULTI_THREAD= 'MSK_OFF' ;
  param.MSK_IPAR_NUM_THREADS=1;
  param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
  
  s_time = tic;
  res = mskqpopt(Q,c,A,clb,cub,lb,ub,param);

  if(isfield(res,'sol') && isfield(res.sol,'itr') && isfield(res.sol.itr,'pobjval'))
    sol.sol_obj_val = res.sol.itr.pobjval;
  else
    sol.sol_obj_val = inf;
  end
  sol.rac_time = toc(s_time);
  sol.n_iter = -inf;
  sol.sol_residue_p = -inf;
  sol.sol_residue_d = -inf;

end









