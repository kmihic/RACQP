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


function run_admm_compare_test(r_time,  rnd_seed, print_graph)

solver = 'racqp';
 
  solutions = [];
  sol_stat = [];
  dir = "../data/data_rnd/";
  filename = "MARKOWITZ_N3000_SP0_95_E0_5_Z2_C100_RND1";
  beta = 1;
  disp("Solving: "+filename)
  disp("LOADING THE MODEL...")
  load(dir+filename);
  epsilon = 1e-20;
  iterations = [10 50 100];
  admm_v = get_admm_type();

  for kk = 1:size(admm_v,1)
    for n_iter = iterations  
      if(strcmpi(admm_v(kk,1),'DADMM'))
        N = model.size;
        model.local_constraints.lb=model.lb;
        model.local_constraints.ub=model.ub;
        model.local_constraints.Aineq = sparse(0,N);
        model.local_constraints.bineq = [];
        model.local_constraints.Aeq = sparse(0,N);
        model.local_constraints.beq = [];
        model.lb = -Inf(N,1);
        model.ub = Inf(N,1);
      end
      if(strcmpi(admm_v(kk,1),'CADMM2'))
        % single block mode - one block for x, the other for x_hat -- split method
        % to handle bounds
        n_grp = 1;
        model.group_mode = 'CADMM';
      else
        n_grp = 50;
        model.group_mode = admm_v(kk,1);
      end
      disp("ADMM variant: "+admm_v(kk,1));
      disp("Number of groups: "+n_grp);
      disp("Max number of iterations: "+n_iter);
      inst_param = get_instance_run_params(model, 'racqp', r_time, false, ...
                        epsilon, n_iter);
      inst_param.racqp_run_p = get_markowitz_run_params(n_grp, beta, ...
                   r_time, epsilon, n_iter, rnd_seed, model.Q,print_graph);
      if(strcmpi(admm_v(kk,1),'DADMM'))
        inst_param.racqp_run_p.sub_solver_type = 'gurobi';
      end
      sol = solve_instance(inst_param);
      sol.name = admm_v(kk,1)+"_NITER_"+n_iter;      
      solutions = [solutions,sol]; 
    end %for each n iter
  end %for each admm variant
  disp(" ")
  disp("####################")
  disp(" SUMMARY ")
  print_solutions(solutions, false, false, 'ADMM_variant_Niter');
  if(print_graph)
     plot_admm_variants_evol(solutions);
  end

end

function run_params = get_markowitz_run_params(n_blocks, beta, max_time, epsilon,...
                                   max_iter, rnd_seed,Q, debug )

  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  run_params.debug = debug;
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
           "DADMM"];

end


