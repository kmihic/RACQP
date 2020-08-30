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


function run_admm_comparison_test_bin()

  addpath('/opt/gurobi/gurobi811/linux64/matlab/');
  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../utils');
 solver = 'racqp';
 
  solutions = [];
  sol_stat = [];
  beta = 50;
  rnd_seed = 123456;
  r_time = Inf;
  epsilon = 0;
  iterations = 200;
  n_size = 1000;  
  n_grp = 50;
  r = round(n_size/2);
  k = 1e-5;
  admm_v = get_admm_type(); 
  model = get_model(n_size,r,k);
  for kk = 1:size(admm_v,1)
    for n_iter = iterations  
      model.group_mode = admm_v(kk,1);
      disp("ADMM variant: "+admm_v(kk,1));
      disp("Number of groups: "+n_grp);
      disp("Max number of iterations: "+n_iter);
      inst_param = get_instance_run_params(model, 'racqp', r_time, false, ...
                        epsilon, n_iter);
      run_param = get_markowitz_run_params(n_grp, beta, ...
                   r_time, epsilon, n_iter, rnd_seed, model.Q);
      inst_param.racqp_run_p = default_mip_run_params(n_size, 1, r_time,run_param);      
      % testing a single loop, no perturb
      inst_param.racqp_run_p.n_perturb_trial = inf;
      inst_param.racqp_run_p.mip_epsilon = 0;  
     inst_param.racqp_run_p.debug = 1;    
      inst_param.racqp_run_p.max_iter = 1;
      sol = solve_instance(inst_param);
      sol.name = admm_v(kk,1)+"_NITER_"+n_iter;      
      solutions = [solutions,sol]; 
    end %for each n iter
  end %for each admm variant
  plot_admm_variants_evol_bin(solutions);

end

function run_params = get_markowitz_run_params(n_blocks, beta, max_time, epsilon,...
                                   max_iter, rnd_seed,Q)

  run_params = default_run_params(n_blocks, beta);
  %change some parameters
  run_params.max_iter = max_iter;
  run_params.rnd_seed = rnd_seed;
  %turn on verbose, need to get plot data
  run_params.debug = true;
  %get density
  density = nnz(Q)/size(Q,1)^2;
  if(density > 0.5)
    run_params.use_sparse = false; %default value is true
  end
  gurobi_params = default_gurobi_parameters;
  run_params.gurobi_params = gurobi_params;
end

function r_val = get_admm_type()

  r_val = ["RAC";
           "RP";
           "CADMM";
           "DADMM"
];

end

function m = get_model(n,r,k)

  rng(1234);
  m.size = n;
  m.Q = construct_matrix(n,n, 10) + k*eye(n,n);
  m.c = rand(n,1);
  m.Aineq = sparse(0,n);
  m.bineq = [];
  m.x0 = round(rand(n,1));
  m.integers = [];
  m.binary = [1:n];
  m.lb = zeros(n,1);
  m.ub = Inf(n,1);
  m.const = 0;
  m.Aeq = ones(1,n);
  m.beq = r;
end


function A = construct_matrix(n_size, m_rank, kappa)


  P = orth(randn(n_size));
  lambda = zeros(n_size,1);
  lambda(1:m_rank) = kappa*rand(m_rank,1);
  D = diag(lambda);
  A = P*D*P';

end




