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

function print_solutions(solutions, show_time, show_all, first_col_title, msg)

  if(nargin <= 1)
    show_time = false;
  end
  if(nargin <= 2)
    show_all = 0;
  end
  if(nargin <= 3)
   first_col_title = 'Instance_name';
  end
  if(nargin <= 4)
    msg = "";
  end

  name = [];
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
  if(show_time)
    i_time = zeros(N,1);  
    s_time = zeros(N,1);
    m_time = zeros(N,1);
  end
  for ii = 1: N
    name = [name;solutions(ii).name];
    run_time(ii) = solutions(ii).runtime;
    n_iter(ii) = solutions(ii).n_iter;
    p_res( ii) = solutions(ii).sol_res_p.abs;
    p_res_rel(ii) = solutions(ii).sol_res_p.rel;
    p_res_L2(ii) = solutions(ii).sol_res_p.L2;
    d_res(ii) = solutions(ii).sol_res_d.abs;
    d_res_rel(ii) = solutions(ii).sol_res_d.rel;
    d_res_L2(ii) = solutions(ii).sol_res_d.L2;
    objval(ii) = solutions(ii).sol_obj_val;
    if(show_time)
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
  end
  if(show_time)
    T = make_table(first_col_title, name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all, ...
                    i_time, s_time, m_time);
  else
    T = make_table(first_col_title, name, run_time, n_iter, p_res, p_res_rel, p_res_L2, ...
                    d_res, d_res_rel, d_res_L2, objval, show_all);
  end
  disp(" ")
  disp("####")
  disp(msg);
  disp(" ");
  disp(T);
end



