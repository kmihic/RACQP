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


function sol = use_osqp(m, max_time, epsilon, max_iter)

  m_os = osqp;

  settings = m_os.default_settings();
  settings.eps_abs = epsilon;
  settings.eps_rel = epsilon;
  settings.polish = false;
  settings.verbose = true;
  settings.max_iter = max_iter;
  if(max_time > 0)
    settings.time_limit = max_time;
  end
  %settings.scaling = false;

  % osqp does not return dual variable for variable 
  % lower/upper bounds. Converting bounds to constraints 
  % to get duals for dual residual to match Gurobi.
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
  l = [m.beq;-inf(length(m.bineq),1)];
  u = [m.beq;m.bineq];

  x0 = zeros(m.size,1);
  y0 = zeros(size(A,1),1);
  m_os.setup(Q, c, A, l, u, settings);
  m_os.warm_start('x', x0, 'y', y0);


  res = m_os.solve();

% prepare the output to be read by run*test scripts
 % if(strcmpi(res.info.status, 'solved'))
    sol.sol_x = res.x;
    sol.sol_y = res.y;
    sol.sol_obj_val = res.info.obj_val;    
    sol.n_iter = res.info.iter;
    
    % get residuals
    % turning off dual for bounds; having them as inequalities
    do_bounds = false;        
    % osqp does Lagrangian funct f(x) + y'A, we do -y'A
    [res_p res_d ] = get_residuals(m, res.x, -res.y, true, false, ...
                          [], do_bounds);
    sol.sol_res_p = res_p;
    sol.sol_res_d = res_d;
%   else
%     sol.sol_x = [];
%     sol.sol_y = [];
%     sol.sol_obj_val = inf;
%     sol.n_iter = -inf;
%     sol.sol_res_p.abs = -inf;
%     sol.sol_res_d.abs = -inf;
%     sol.sol_res_d.L2 = -inf;
%     sol.sol_res_p.rel = -inf;
%     sol.sol_res_d.rel = -inf;
%     sol.sol_res_d.L2 = -inf;
%   end
  sol.runtime = res.info.run_time;

end
