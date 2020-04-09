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

function sol = get_stats(solutions)

  N = length(solutions);
  run_time = zeros(N,1);
  n_iter = zeros(N,1);
  p_res = zeros(N,1);
  p_res_rel = zeros(N,1);
  p_res_L2 = zeros(N,1);
  d_res = zeros(N,1); 
  d_res_rel = zeros(N,1);
  d_res_L2 = zeros(N,1);
  objval = zeros(N,1);
  i_time = zeros(N,1);  
  s_time = zeros(N,1);
  m_time = zeros(N,1);
  for ii = 1: N
    run_time(ii) = solutions(ii).runtime;
    n_iter(ii) = solutions(ii).n_iter;
    p_res( ii) = solutions(ii).sol_res_p.abs;
    p_res_rel(ii) = solutions(ii).sol_res_p.rel;
    p_res_L2(ii) = solutions(ii).sol_res_p.L2;
    d_res(ii) = solutions(ii).sol_res_d.abs;
    d_res_rel(ii) = solutions(ii).sol_res_d.rel;
    d_res_L2(ii) = solutions(ii).sol_res_d.L2;
    objval(ii) = solutions(ii).sol_obj_val;
    % RAC runtime detailed info
    if(isfield(solutions(ii),'init_time'))
      i_time(ii) = solutions(ii).init_time;
    end
    if(isfield(solutions(ii),'solver_time'))
      s_time(ii) = solutions(ii).solver_time;
    end
    if(isfield(solutions(ii),'model_time'))
      m_time(ii) = solutions(ii).model_time;
    end
  end
  %mean
  s.run_time = mean(run_time);
  s.n_iter = mean(n_iter);
  s.p_res = mean(p_res);
  s.p_res_rel = mean(p_res_rel);
  s.d_res = mean(d_res);
  s.d_res_rel = mean(d_res_rel);
  s.p_res_L2 = mean(p_res_L2);
  s.d_res_L2 = mean(d_res_L2);
  s.objval = mean(objval);
  s.init_time = mean(i_time);
  s.solver_time = mean(s_time);
  s.model_time = mean(m_time);
  sol.mean = s;
  %min
  s.run_time = min(run_time);
  s.n_iter = min(n_iter);
  s.p_res = min(p_res);
  s.p_res_rel = min(p_res_rel);
  s.d_res = min(d_res);
  s.d_res_rel = min(d_res_rel);
  s.p_res_L2 = min(p_res_L2);
  s.d_res_L2 = min(d_res_L2);
  s.objval = min(objval);
  s.init_time = min(i_time);
  s.solver_time = min(s_time);
  s.model_time = min(m_time);
  sol.min = s;
  %max
  s.run_time = max(run_time);
  s.n_iter = max(n_iter);
  s.p_res = max(p_res);
  s.p_res_rel = max(p_res_rel);
  s.d_res = max(d_res);
  s.d_res_rel = max(d_res_rel);
  s.p_res_L2 = max(p_res_L2);
  s.d_res_L2 = max(d_res_L2);
  s.objval = max(objval);
  s.init_time = max(i_time);
  s.solver_time = max(s_time);
  s.model_time = max(m_time);
  sol.max = s;
  %std
  s.run_time = std(run_time);
  s.n_iter = std(n_iter);
  s.p_res = std(p_res);
  s.p_res_rel = std(p_res_rel);
  s.d_res = std(d_res);
  s.d_res_rel = std(d_res_rel);
  s.p_res_L2 = std(p_res_L2);
  s.d_res_L2 = std(d_res_L2);
  s.objval = std(objval);
  s.init_time = std(i_time);
  s.solver_time = std(s_time);
  s.model_time = std(m_time);
  sol.std = s;
end