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

function print_stats(s, show_time, show_all, first_col_title)

  if(nargin <= 1)
    show_time = false;
  end
  if(nargin <= 2)
    show_all = false;
  end
  if(nargin <= 3)
   first_col_title = 'Instances';
  end

  disp(" ")
  disp("#####################")
  N = length(s);

  disp(" ")  
  disp(' MIN ')
  name = [];
  run_time = zeros(N,1);
  n_iter = zeros(N,1);
  p_res = zeros(N,1);
  p_res_rel = zeros(N,1);
  d_res = zeros(N,1); 
  d_res_rel = zeros(N,1);
  p_res_L2 = zeros(N,1);
  d_res_L2 = zeros(N,1);
  objval = zeros(N,1);
  if(show_time)
    i_time = zeros(N,1);  
    s_time = zeros(N,1);
    m_time = zeros(N,1);
  end
  for ii = 1: N
    name = [name;s(ii).name];
    run_time(ii) = s(ii).min.run_time;
    n_iter(ii) = s(ii).min.n_iter;

    p_res( ii) = s(ii).min.p_res; %abs
    p_res_rel( ii) = s(ii).min.p_res_rel;
    d_res(ii) = s(ii).min.d_res; %abs
    d_res_rel(ii) = s(ii).min.d_res_rel;     
    p_res_L2( ii) = s(ii).min.p_res_L2;
    d_res_L2(ii) = s(ii).min.d_res_L2;    
    objval(ii) = s(ii).min.objval;
    if(show_time)      
      i_time(ii) = s(ii).min.init_time;
      s_time(ii) = s(ii).min.solver_time;
      m_time(ii) = s(ii).min.model_time;
    end
  end
  if(show_time)
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all, ...
                    i_time, s_time, m_time);
  else
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all);
  end
  disp(T);

  disp(" ")  
  disp(' MAX ')
  name = [];
  run_time = zeros(N,1);
  n_iter = zeros(N,1);
  p_res = zeros(N,1);
  p_res_rel = zeros(N,1);
  d_res = zeros(N,1); 
  d_res_rel = zeros(N,1);
  p_res_L2 = zeros(N,1);
  d_res_L2 = zeros(N,1);
  objval = zeros(N,1);  
  if(show_time)
    i_time = zeros(N,1);  
    s_time = zeros(N,1);
    m_time = zeros(N,1);
  end
  for ii = 1: N
    name = [name;s(ii).name];
    run_time(ii) = s(ii).max.run_time;
    n_iter(ii) = s(ii).max.n_iter;
    p_res( ii) = s(ii).max.p_res;
    p_res_rel( ii) = s(ii).max.p_res_rel;
    d_res(ii) = s(ii).max.d_res;
    d_res_rel(ii) = s(ii).max.d_res_rel;     
    p_res_L2( ii) = s(ii).max.p_res_L2;
    d_res_L2(ii) = s(ii).max.d_res_L2;      
    objval(ii) = s(ii).max.objval;
    if(show_time)      
      i_time(ii) = s(ii).max.init_time;
      s_time(ii) = s(ii).max.solver_time;
      m_time(ii) = s(ii).max.model_time;
    end
  end
  if(show_time)
    T = make_table(first_col_title, name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all, ...
                    i_time, s_time, m_time);
  else
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all);
  end
  disp(T);

  disp(" ")  
  disp(' STDEV ')
  name = [];
  run_time = zeros(N,1);
  n_iter = zeros(N,1);
  p_res = zeros(N,1);
  p_res_rel = zeros(N,1);
  d_res = zeros(N,1); 
  d_res_rel = zeros(N,1);
  p_res_L2 = zeros(N,1);
  d_res_L2 = zeros(N,1);
  objval = zeros(N,1);
  if(show_time)
    i_time = zeros(N,1);  
    s_time = zeros(N,1);
    m_time = zeros(N,1);
  end
  for ii = 1: N
    name = [name;s(ii).name];
    run_time(ii) = s(ii).std.run_time;
    n_iter(ii) = s(ii).std.n_iter;
    p_res( ii) = s(ii).std.p_res;
    p_res_rel( ii) = s(ii).std.p_res_rel;
    d_res(ii) = s(ii).std.d_res;
    d_res_rel(ii) = s(ii).std.d_res_rel;    
    p_res_L2( ii) = s(ii).std.p_res_L2;
    d_res_L2(ii) = s(ii).std.d_res_L2;       
    objval(ii) = s(ii).std.objval;
    if(show_time)      
      i_time(ii) = s(ii).std.init_time;
      s_time(ii) = s(ii).std.solver_time;
      m_time(ii) = s(ii).std.model_time;
    end
  end
  if(show_time)
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all, ...
                    i_time, s_time, m_time);
  else
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all);
  end
  disp(T);

  disp(" ")  
  disp(' MEAN ')
  name = [];
  run_time = zeros(N,1);
  n_iter = zeros(N,1);
  p_res = zeros(N,1);
  p_res_rel = zeros(N,1);
  d_res = zeros(N,1); 
  d_res_rel = zeros(N,1);
  p_res_L2 = zeros(N,1);
  d_res_L2 = zeros(N,1);
  objval = zeros(N,1);
  if(show_time)
    i_time = zeros(N,1);  
    s_time = zeros(N,1);
    m_time = zeros(N,1);
  end
  for ii = 1: N
    name = [name;s(ii).name];
    run_time(ii) = s(ii).mean.run_time;
    n_iter(ii) = s(ii).mean.n_iter;
    p_res( ii) = s(ii).mean.p_res;
    p_res_rel( ii) = s(ii).mean.p_res_rel;
    d_res(ii) = s(ii).mean.d_res;
    d_res_rel(ii) = s(ii).mean.d_res_rel;       
    p_res_L2( ii) = s(ii).mean.p_res_L2;
    d_res_L2(ii) = s(ii).mean.d_res_L2;    
    objval(ii) = s(ii).mean.objval;
    if(show_time)      
      i_time(ii) = s(ii).mean.init_time;
      s_time(ii) = s(ii).mean.solver_time;
      m_time(ii) = s(ii).mean.model_time;
    end
  end
  if(show_time)
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all, ...
                    i_time, s_time, m_time);
  else
    T = make_table(first_col_title,name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all);
  end
  disp(T);
end


