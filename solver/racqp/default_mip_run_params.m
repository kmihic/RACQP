function run_mip = default_mip_run_params(nsize, lambda, max_time,run_params, p_trial)

  if(nargin <=4)
   p_trial = 1;
  end
  
  %load default gurobi parameters
  gurobi_params = default_gurobi_parameters;
  %change how long gurobi will spend on each subproblem
  gurobi_params.TimeLimit = 2;
  run_params.gurobi_params = gurobi_params;
  
  %set MIP run parameters
  run_mip.permute_type = 'swap';
  run_mip.permute_dist = 'exponential'; 
  run_mip.permute_min = 2; 
  run_mip.permute_max = nsize; 
  run_mip.permute_mu = lambda*nsize; 
  run_mip.max_iter = Inf;
  run_mip.max_nperturb = Inf;
  run_mip.max_rtime = max_time; 
  run_mip.rnd_seed = 123; 
  run_mip.run_sub = run_params;
  run_mip.debug = 0;
  % testing a single loop, no perturb
  run_mip.n_perturb_trial = round(p_trial*nsize);
  run_mip.mip_epsilon = run_params.epsilon_prim;
  
end